# -*- coding: utf-8 -*-
"""
Aplicação desktop para processamento geoespacial em Python/Tkinter.

Abas:
1) Grade de pontos a partir de shapefile de polígonos.
2) Buffers a partir de pontos (diâmetro em m ou km).
3) Recorte de rasters por feição de buffers (GeoTIFF UTM), em subpastas por ano.
4) Métricas de paisagem (PyLandStats e/ou landscapemetrics do R via rpy2) nos níveis patch/classe/paisagem (vizinhança 4 ou 8).
5) Junção de buffers com métricas em shapefile.
6) Reprojeção de rasters em lote (sugestão automática de EPSG SIRGAS 2000/UTM).
7) Fundido (em memória): recorte + métricas sem salvar o recorte, paralelo e retomável.
8) Consolidar tabelas: junta CSVs (e/ou Parquets) numa tabela única.

Histórico de versões
--------------------
V6  - Recorte e métricas paralelizados; subpastas por ano; nomes padronizados.
V9  - Abas novas: Reprojeção de rasters e Fundido em memória (recorte+métrica sem disco).
V10 - Fundido com ProcessPoolExecutor (escapa do GIL), 1 tarefa por raster, cache do GDAL,
      opção de pular reprojeção quando a origem já está em metros.
V11 - Barra de progresso por feição e estimativa de tempo (ETA) no Fundido.
V12 - Correção do erro de exceção diferida (NameError 'exc'); retomada inicial.
V13 - Retomada do Fundido baseada nos CSVs já existentes (pasta + subpastas);
      barra/contagem/ETA globais; robustez (um raster com erro não derruba o lote;
      tratamento de BrokenProcessPool); escrita atômica dos CSVs.
V14 - Consolidação movida para a aba 8 (Consolidar); Fundido só grava CSV por feição.
V15 - Consolidação com DuckDB (lê todos os CSVs de uma vez, paralelo): Parquet único,
      Parquet por ano e CSV; fallback em pyarrow/concatenação por bytes.
V16 - Consolidação em lotes mostrando número de arquivos e ETA na leitura.
V17 (atual) -
      * Aba Fundido grava Parquet POR ANO (um arquivo em cada pasta de ano), em fluxo
        (baixa memória), e gera ao final um PARQUET ÚNICO com todos os dados.
      * Retomada do Fundido passa a reconhecer os Parquets por ano já gravados.
      * Mantida a opção de salvar CSV por feição no Fundido.
      * Aba Consolidar ganha a opção de LER vários Parquets e juntá-los num único Parquet
        (além de continuar consolidando CSVs).
      * Aba 9 (Buffer + Métrica) aceita métricas em CSV ou Parquet, pareando pelo campo
        id_pt (id) de ambos para integrar à geometria do shapefile.

V18 -
      * Aba 6 (Métricas) ganha o motor R `landscapemetrics` (via rpy2), além do PyLandStats.
        O usuário pode marcar um motor, o outro, ou ambos. Painel de seleção de métricas por
        nível (patch/classe/paisagem) para o R, com as funções lsm_* do pacote. Campo cellsize
        (options(landscapemetrics_cellsize = 30) por padrão), vizinhança via `directions` (4/8),
        leitura recursiva de subpastas e núcleos paralelos (ProcessPoolExecutor, sessão R por
        processo, com fallback sequencial). Mantido o fluxo de salvamento (CSV por raster em
        subpasta do ANO + CSV geral, agora com prefixo do motor e coluna `engine`), o
        detalhamento dos logs e as barras de progresso.

Execução: python app_grade_buffer_shapefile_v17_pylandstats.py
"""

import math
import os
import queue
import re
import shutil
import threading
import time
import json
import hashlib
import multiprocessing
import unicodedata
import zipfile
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
from rasterio.transform import array_bounds
from rasterio.warp import Resampling, calculate_default_transform, reproject, transform_bounds

# Usa pyogrio como motor de leitura/gravação para evitar dependência direta do Fiona via pip.
gpd.options.io_engine = "pyogrio"
from pyproj import CRS
from shapely.geometry import Point, mapping
from shapely.ops import unary_union
from shapely.prepared import prep

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
except ImportError as exc:
    raise RuntimeError("Tkinter não está disponível nesta instalação do Python.") from exc


# =============================================================================
# Métricas PyLandStats disponíveis na versão 3.x
# =============================================================================

PATCH_METRICS = [
    "area",
    "perimeter",
    "perimeter_area_ratio",
    "shape_index",
    "fractal_dimension",
    "core_area",
    "number_of_core_areas",
    "core_area_index",
    "euclidean_nearest_neighbor",
]

# PyLandStats cria métricas de distribuição a partir das métricas de patch.
_PATCH_METRICS_FOR_DISTR = [m for m in PATCH_METRICS if m != "number_of_core_areas"] + ["disjunct_core_area"]
DISTR_METRICS = [
    f"{metric}_{suffix}"
    for metric in _PATCH_METRICS_FOR_DISTR
    for suffix in ["mn", "am", "md", "ra", "sd", "cv"]
]

CLASS_METRICS = [
    "total_area",
    "proportion_of_landscape",
    "number_of_patches",
    "patch_density",
    "largest_patch_index",
    "total_edge",
    "edge_density",
    "total_core_area",
    "core_area_proportion_of_landscape",
    "number_of_disjunct_core_areas",
    "landscape_shape_index",
    "effective_mesh_size",
] + DISTR_METRICS

ENTROPY_METRICS = [
    "entropy",
    "shannon_diversity_index",
    "joint_entropy",
    "conditional_entropy",
    "mutual_information",
    "relative_mutual_information",
    "contagion",
]

LANDSCAPE_METRICS = [
    "total_area",
    "number_of_patches",
    "patch_density",
    "largest_patch_index",
    "total_edge",
    "edge_density",
    "total_core_area",
    "number_of_disjunct_core_areas",
    "landscape_shape_index",
    "effective_mesh_size",
] + ENTROPY_METRICS + DISTR_METRICS


# =============================================================================
# Métricas do pacote R  landscapemetrics  (código -> função lsm_*)
# Fonte: tabelas fornecidas pelo usuário (níveis patch / class / landscape).
# Cada item é (CÓDIGO, "lsm_x_...") e é passado ao calculate_lsm(what=...).
# =============================================================================

R_PATCH_METRICS = [
    ("AREA", "lsm_p_area"), ("CAI", "lsm_p_cai"), ("CIRCLE", "lsm_p_circle"),
    ("CONTIG", "lsm_p_contig"), ("CORE", "lsm_p_core"), ("ENN", "lsm_p_enn"),
    ("FRAC", "lsm_p_frac"), ("GYRATE", "lsm_p_gyrate"), ("NCORE", "lsm_p_ncore"),
    ("PARA", "lsm_p_para"), ("PERIM", "lsm_p_perim"), ("SHAPE", "lsm_p_shape"),
]

R_CLASS_METRICS = [
    ("AI", "lsm_c_ai"), ("AREA_CV", "lsm_c_area_cv"), ("AREA_MN", "lsm_c_area_mn"),
    ("AREA_SD", "lsm_c_area_sd"), ("CA", "lsm_c_ca"), ("CAI_CV", "lsm_c_cai_cv"),
    ("CAI_MN", "lsm_c_cai_mn"), ("CAI_SD", "lsm_c_cai_sd"), ("CIRCLE_CV", "lsm_c_circle_cv"),
    ("CIRCLE_MN", "lsm_c_circle_mn"), ("CIRCLE_SD", "lsm_c_circle_sd"), ("CLUMPY", "lsm_c_clumpy"),
    ("COHESION", "lsm_c_cohesion"), ("CONTIG_CV", "lsm_c_contig_cv"), ("CONTIG_MN", "lsm_c_contig_mn"),
    ("CONTIG_SD", "lsm_c_contig_sd"), ("CORE_CV", "lsm_c_core_cv"), ("CORE_MN", "lsm_c_core_mn"),
    ("CORE_SD", "lsm_c_core_sd"), ("CPLAND", "lsm_c_cpland"), ("DCAD", "lsm_c_dcad"),
    ("DCORE_CV", "lsm_c_dcore_cv"), ("DCORE_MN", "lsm_c_dcore_mn"), ("DCORE_SD", "lsm_c_dcore_sd"),
    ("DIVISION", "lsm_c_division"), ("ED", "lsm_c_ed"), ("ENN_CV", "lsm_c_enn_cv"),
    ("ENN_MN", "lsm_c_enn_mn"), ("ENN_SD", "lsm_c_enn_sd"), ("FRAC_CV", "lsm_c_frac_cv"),
    ("FRAC_MN", "lsm_c_frac_mn"), ("FRAC_SD", "lsm_c_frac_sd"), ("GYRATE_CV", "lsm_c_gyrate_cv"),
    ("GYRATE_MN", "lsm_c_gyrate_mn"), ("GYRATE_SD", "lsm_c_gyrate_sd"), ("IJI", "lsm_c_iji"),
    ("LPI", "lsm_c_lpi"), ("LSI", "lsm_c_lsi"), ("MESH", "lsm_c_mesh"), ("NDCA", "lsm_c_ndca"),
    ("nLSI", "lsm_c_nlsi"), ("NP", "lsm_c_np"), ("PAFRAC", "lsm_c_pafrac"),
    ("PARA_CV", "lsm_c_para_cv"), ("PARA_MN", "lsm_c_para_mn"), ("PARA_SD", "lsm_c_para_sd"),
    ("PD", "lsm_c_pd"), ("PLADJ", "lsm_c_pladj"), ("PLAND", "lsm_c_pland"),
    ("SHAPE_CV", "lsm_c_shape_cv"), ("SHAPE_MN", "lsm_c_shape_mn"), ("SHAPE_SD", "lsm_c_shape_sd"),
    ("SPLIT", "lsm_c_split"), ("TCA", "lsm_c_tca"), ("TE", "lsm_c_te"),
]

R_LANDSCAPE_METRICS = [
    ("AI", "lsm_l_ai"), ("AREA_CV", "lsm_l_area_cv"), ("AREA_MN", "lsm_l_area_mn"),
    ("AREA_SD", "lsm_l_area_sd"), ("CAI_CV", "lsm_l_cai_cv"), ("CAI_MN", "lsm_l_cai_mn"),
    ("CAI_SD", "lsm_l_cai_sd"), ("CIRCLE_CV", "lsm_l_circle_cv"), ("CIRCLE_MN", "lsm_l_circle_mn"),
    ("CIRCLE_SD", "lsm_l_circle_sd"), ("COHESION", "lsm_l_cohesion"), ("CONDENT", "lsm_l_condent"),
    ("CONTAG", "lsm_l_contag"), ("CONTIG_CV", "lsm_l_contig_cv"), ("CONTIG_MN", "lsm_l_contig_mn"),
    ("CONTIG_SD", "lsm_l_contig_sd"), ("CORE_CV", "lsm_l_core_cv"), ("CORE_MN", "lsm_l_core_mn"),
    ("CORE_SD", "lsm_l_core_sd"), ("DCAD", "lsm_l_dcad"), ("DCORE_CV", "lsm_l_dcore_cv"),
    ("DCORE_MN", "lsm_l_dcore_mn"), ("DCORE_SD", "lsm_l_dcore_sd"), ("DIVISION", "lsm_l_division"),
    ("ED", "lsm_l_ed"), ("ENN_CV", "lsm_l_enn_cv"), ("ENN_MN", "lsm_l_enn_mn"),
    ("ENN_SD", "lsm_l_enn_sd"), ("ENT", "lsm_l_ent"), ("FRAC_CV", "lsm_l_frac_cv"),
    ("FRAC_MN", "lsm_l_frac_mn"), ("FRAC_SD", "lsm_l_frac_sd"), ("GYRATE_CV", "lsm_l_gyrate_cv"),
    ("GYRATE_MN", "lsm_l_gyrate_mn"), ("GYRATE_SD", "lsm_l_gyrate_sd"), ("IJI", "lsm_l_iji"),
    ("JOINENT", "lsm_l_joinent"), ("LPI", "lsm_l_lpi"), ("LSI", "lsm_l_lsi"), ("MESH", "lsm_l_mesh"),
    ("MSIDI", "lsm_l_msidi"), ("MSIEI", "lsm_l_msiei"), ("MUTINF", "lsm_l_mutinf"),
    ("NDCA", "lsm_l_ndca"), ("NP", "lsm_l_np"), ("PAFRAC", "lsm_l_pafrac"),
    ("PARA_CV", "lsm_l_para_cv"), ("PARA_MN", "lsm_l_para_mn"), ("PARA_SD", "lsm_l_para_sd"),
    ("PD", "lsm_l_pd"), ("PLADJ", "lsm_l_pladj"), ("PR", "lsm_l_pr"), ("PRD", "lsm_l_prd"),
    ("RELMUTINF", "lsm_l_relmutinf"), ("RPR/RPD", "lsm_l_rpr"), ("SHAPE_CV", "lsm_l_shape_cv"),
    ("SHAPE_MN", "lsm_l_shape_mn"), ("SHAPE_SD", "lsm_l_shape_sd"), ("SHDI", "lsm_l_shdi"),
    ("SHEI", "lsm_l_shei"), ("SIDI", "lsm_l_sidi"), ("SIEI", "lsm_l_siei"), ("SPLIT", "lsm_l_split"),
    ("TA", "lsm_l_ta"), ("TCA", "lsm_l_tca"), ("TE", "lsm_l_te"),
]

# Rótulos exibidos nas listboxes do R ("CÓDIGO — lsm_x_...").
R_METRICS_BY_LEVEL = {"patch": R_PATCH_METRICS, "class": R_CLASS_METRICS, "landscape": R_LANDSCAPE_METRICS}
R_LEVEL_ARG = {"patch": "patch", "class": "class", "landscape": "landscape"}


# =============================================================================
# Funções utilitárias gerais
# =============================================================================

class ProcessamentoInterrompido(Exception):
    """Exceção usada quando o usuário clica no botão Parar."""


def checar_parada(stop_event: Optional[threading.Event]):
    if stop_event is not None and stop_event.is_set():
        raise ProcessamentoInterrompido("Processamento interrompido pelo usuário.")


def normalizar_nome(texto: str, max_len: int = 90) -> str:
    """
    Normaliza texto para nome seguro de arquivo:
    - remove acentos;
    - troca espaços por _;
    - remove símbolos problemáticos.
    """
    texto = unicodedata.normalize("NFKD", str(texto))
    texto = texto.encode("ascii", "ignore").decode("ascii")
    texto = re.sub(r"[^\w\-.]+", "_", texto)
    texto = re.sub(r"_+", "_", texto).strip("_")
    if not texto:
        texto = "sem_nome"
    return texto[:max_len] if len(texto) > max_len else texto


def rotulo_valor(valor: float, unidade: str) -> str:
    """Cria rótulo curto para distância: 500m, 1p5km etc."""
    if float(valor).is_integer():
        txt = str(int(valor))
    else:
        txt = str(valor).replace(".", "p")
    return f"{txt}{unidade}"


def converter_para_metros(valores: Iterable[float], unidade: str) -> List[float]:
    """Converte valores em m ou km para metros."""
    fator = 1000.0 if unidade.lower() == "km" else 1.0
    return [float(v) * fator for v in valores]


def parse_lista_numerica(texto: str) -> List[float]:
    """
    Lê lista de números digitada pelo usuário.
    Aceita separação por vírgula, ponto e vírgula, espaço ou quebra de linha.
    Exemplos:
        100, 250, 500
        0.5; 1; 2
    """
    if not texto or not texto.strip():
        raise ValueError("Informe pelo menos um valor de distância/diâmetro.")

    partes = re.split(r"[,\n;\t ]+", texto.strip())
    valores = []
    for item in partes:
        if not item:
            continue
        item = item.replace(",", ".") if item.count(",") == 1 and "." not in item else item
        try:
            valor = float(item)
        except ValueError as exc:
            raise ValueError(f"Valor inválido: {item}") from exc
        if valor <= 0:
            raise ValueError("Todos os valores devem ser maiores que zero.")
        valores.append(valor)

    # Remove duplicados preservando a ordem.
    saida = []
    for v in valores:
        if v not in saida:
            saida.append(v)
    return saida


def unir_geometrias(gdf: gpd.GeoDataFrame):
    """Une geometrias de modo compatível com diferentes versões do GeoPandas."""
    try:
        return gdf.geometry.union_all()
    except AttributeError:
        return unary_union(gdf.geometry)


def ler_shapefile(caminho: str) -> gpd.GeoDataFrame:
    """Lê shapefile e faz validações mínimas."""
    caminho = Path(caminho)
    if not caminho.exists() or caminho.suffix.lower() != ".shp":
        raise FileNotFoundError(f"Arquivo .shp inválido: {caminho}")

    gdf = gpd.read_file(caminho, engine="pyogrio")

    if gdf.empty:
        raise ValueError(f"O shapefile está vazio: {caminho.name}")

    if gdf.crs is None:
        raise ValueError(
            f"O shapefile {caminho.name} não possui CRS definido. "
            "Defina a projeção antes de processar."
        )

    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()

    if gdf.empty:
        raise ValueError(f"O shapefile não possui geometrias válidas: {caminho.name}")

    return gdf


def obter_crs_metrico(gdf: gpd.GeoDataFrame) -> CRS:
    """
    Retorna um CRS métrico adequado.
    - Se o CRS de entrada já for projetado em metros, usa o CRS original.
    - Se for geográfico, calcula automaticamente uma zona UTM pelo centroide.
    """
    crs_original = CRS.from_user_input(gdf.crs)

    if crs_original.is_projected:
        try:
            unidade = crs_original.axis_info[0].unit_name.lower()
            if "metre" in unidade or "meter" in unidade:
                return crs_original
        except Exception:
            return crs_original

    gdf_wgs84 = gdf.to_crs(epsg=4326)
    centroide = unir_geometrias(gdf_wgs84).centroid
    lon, lat = centroide.x, centroide.y

    zona = int(math.floor((lon + 180.0) / 6.0) + 1)
    epsg = 32600 + zona if lat >= 0 else 32700 + zona
    return CRS.from_epsg(epsg)


def obter_crs_utm_por_geometria(gdf: gpd.GeoDataFrame) -> CRS:
    """
    Calcula obrigatoriamente uma projeção UTM WGS84 com base no centroide da geometria.
    """
    if gdf.crs is None:
        raise ValueError("Não é possível calcular UTM porque a geometria não possui CRS.")

    gdf_wgs84 = gdf.to_crs(epsg=4326)
    geom = unir_geometrias(gdf_wgs84)
    centroide = geom.centroid
    lon, lat = centroide.x, centroide.y
    zona = int(math.floor((lon + 180.0) / 6.0) + 1)
    epsg = 32600 + zona if lat >= 0 else 32700 + zona
    return CRS.from_epsg(epsg)


def corrigir_geometrias(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Corrige geometrias inválidas quando possível.
    Usa buffer(0), método clássico e compatível com várias versões do Shapely.
    """
    gdf = gdf.copy()
    invalidas = ~gdf.geometry.is_valid
    if invalidas.any():
        gdf.loc[invalidas, "geometry"] = gdf.loc[invalidas, "geometry"].buffer(0)
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()
    return gdf


def compactar_shapefile(shp_path: Path, zip_path: Path) -> Path:
    """Compacta todos os componentes do shapefile em um arquivo .zip."""
    base = shp_path.with_suffix("")
    if zip_path.exists():
        zip_path.unlink()

    componentes = list(base.parent.glob(base.name + ".*"))
    if not componentes:
        raise FileNotFoundError(f"Nenhum componente encontrado para {shp_path}")

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for componente in componentes:
            zf.write(componente, arcname=componente.name)

    return zip_path


def salvar_gdf_como_shapefile_zip(
    gdf: gpd.GeoDataFrame,
    pasta_saida: str,
    nome_base: str,
    encoding: str = "utf-8",
) -> Path:
    """Salva GeoDataFrame em shapefile temporário e retorna apenas o .zip final."""
    pasta_saida = Path(pasta_saida)
    pasta_saida.mkdir(parents=True, exist_ok=True)

    nome_base = normalizar_nome(nome_base)
    pasta_temp = pasta_saida / f"__tmp_{nome_base}"

    if pasta_temp.exists():
        shutil.rmtree(pasta_temp)
    pasta_temp.mkdir(parents=True, exist_ok=True)

    shp_path = pasta_temp / f"{nome_base}.shp"
    zip_path = pasta_saida / f"{nome_base}.zip"

    gdf.to_file(shp_path, driver="ESRI Shapefile", encoding=encoding, index=False, engine="pyogrio")

    compactar_shapefile(shp_path, zip_path)
    shutil.rmtree(pasta_temp, ignore_errors=True)
    return zip_path


def filtrar_colunas_shapefile(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Mantém colunas compatíveis com shapefile."""
    gdf = gdf.copy()
    for col in list(gdf.columns):
        if col == gdf.geometry.name:
            continue
        if gdf[col].dtype == "object":
            gdf[col] = gdf[col].astype(str).str.slice(0, 250)
    return gdf


def extrair_ano_do_nome(caminho: str) -> str:
    """Extrai o primeiro ano com 4 dígitos encontrado no nome do arquivo."""
    nome = Path(caminho).stem
    match = re.search(r"(?<!\d)(19\d{2}|20\d{2}|21\d{2})(?!\d)", nome)
    return match.group(1) if match else "anoNA"


def extrair_ano_de_recorte(caminho: str) -> str:
    """
    Extrai o ANO de um arquivo recortado no padrão clip_id_<ID>_diam_<DIAM>_<ANO>.

    O ano é sempre o ÚLTIMO campo separado por '_' no nome (antes da extensão).
    Isso evita confundir o diâmetro (ex.: 1900, 2000, 2100) com o ano quando
    ambos têm 4 dígitos — algo que a busca genérica de ano faria errado.
    Para nomes fora desse padrão (ex.: rasters de origem), usa a extração genérica.

    Exemplos:
        clip_id_1_diam_1000_1985.tif  -> 1985
        clip_id_1_diam_2000_1985.tif  -> 1985  (não confunde com o diâmetro 2000)
        clip_id_7_diam_1000_anoNA.tif -> anoNA
    """
    stem = Path(caminho).stem
    # Aceita um prefixo antes de "clip_id_" (ex.: "Reclass_", da aba 5).
    if re.match(r"(?i)^(?:.+_)?clip_id_.+_diam_.+_[^_]+$", stem):
        return stem.rsplit("_", 1)[-1]
    return extrair_ano_do_nome(caminho)


def formatar_valor_atributo(valor) -> str:
    """Formata atributo para uso seguro no nome do arquivo."""
    if valor is None:
        return "semValor"
    try:
        if isinstance(valor, float) and float(valor).is_integer():
            valor = int(valor)
    except Exception:
        pass
    return normalizar_nome(str(valor).replace(".", "p"), max_len=60)


def numero_padrao_de_nucleos() -> int:
    """Número padrão de workers paralelos, deixando folga para o sistema."""
    cpus = os.cpu_count() or 4
    return max(1, min(8, cpus - 1)) if cpus > 1 else 1


def formatar_duracao_horas(segundos: float) -> str:
    """Formata a estimativa de tempo restante (minutos abaixo de 1 h; horas acima)."""
    segundos = max(0.0, float(segundos))
    if segundos < 60:
        return "< 1 min"
    h = int(segundos // 3600)
    m = int((segundos % 3600) // 60)
    if segundos < 3600:
        return f"{m} min"
    horas = segundos / 3600.0
    return f"{horas:.1f} h ({h}h{m:02d}m)"


def extrair_ano_raster(caminho: str) -> Optional[str]:
    """
    Extrai o ANO (4 dígitos, faixa 19xx/20xx/21xx) do NOME de um raster.

    Regra documentada: se houver mais de um padrão de 4 dígitos no nome, usa o
    ÚLTIMO encontrado (o ano costuma vir no fim do nome, ex.: C10_1_reclass_1985).
    Retorna None se nenhum ano for encontrado (o chamador deve avisar e pular).
    """
    nome = Path(caminho).stem
    achados = re.findall(r"(?<!\d)(19\d{2}|20\d{2}|21\d{2})(?!\d)", nome)
    return achados[-1] if achados else None


def sugerir_epsg_utm_sirgas2000(lon: float, lat: float) -> int:
    """
    Sugere um EPSG de UTM a partir do centróide (lon/lat em graus).

    Preferência para o contexto brasileiro: SIRGAS 2000 / UTM.
    - Zonas Sul 17S..25S -> EPSG 31977..31985.
    - Zonas Norte 11N..22N -> EPSG 31965..31976.
    Fora dessas faixas (raro no Brasil), cai para WGS84 / UTM (326xx/327xx).
    """
    zona = int(math.floor((lon + 180.0) / 6.0) + 1)
    zona = max(1, min(60, zona))
    if lat >= 0:
        if 11 <= zona <= 22:
            return 31954 + zona  # 11N->31965 ... 22N->31976
        return 32600 + zona       # fallback WGS84/UTM Norte
    else:
        if 17 <= zona <= 25:
            return 31960 + zona  # 17S->31977 ... 25S->31985
        return 32700 + zona       # fallback WGS84/UTM Sul


def sugerir_epsg_de_raster(caminho: str) -> Optional[int]:
    """Abre o raster, leva o extent para 4326 e sugere o EPSG UTM do centróide."""
    try:
        with rasterio.open(caminho) as src:
            if src.crs is None:
                return None
            try:
                oeste, sul, leste, norte = transform_bounds(src.crs, "EPSG:4326", *src.bounds, densify_pts=21)
            except Exception:
                oeste, sul, leste, norte = src.bounds
            lon = (oeste + leste) / 2.0
            lat = (sul + norte) / 2.0
            return sugerir_epsg_utm_sirgas2000(lon, lat)
    except Exception:
        return None


def calcular_utm_por_feature(gdf: gpd.GeoDataFrame) -> Dict[object, CRS]:
    """
    Calcula, de uma só vez, a projeção UTM WGS84 adequada para cada feição.

    Em vez de chamar ``to_crs(4326)`` uma vez por feição (lento quando há
    muitas feições e muitos rasters), reprojetamos o GeoDataFrame inteiro
    uma única vez e derivamos a zona UTM do centróide de cada geometria.
    """
    if gdf.crs is None:
        raise ValueError("Não é possível calcular UTM porque o shapefile não possui CRS.")

    gdf_wgs84 = gdf.to_crs(epsg=4326)
    utm_por_feature: Dict[object, CRS] = {}
    cache_epsg: Dict[int, CRS] = {}
    for idx_feature, geom in gdf_wgs84.geometry.items():
        if geom is None or geom.is_empty:
            continue
        centroide = geom.centroid
        lon, lat = centroide.x, centroide.y
        zona = int(math.floor((lon + 180.0) / 6.0) + 1)
        epsg = (32600 + zona) if lat >= 0 else (32700 + zona)
        crs = cache_epsg.get(epsg)
        if crs is None:
            crs = CRS.from_epsg(epsg)
            cache_epsg[epsg] = crs
        utm_por_feature[idx_feature] = crs
    return utm_por_feature


def reservar_caminho_unico(pasta: Path, base: str, extensao: str, reservados: set, lock: threading.Lock) -> Path:
    """
    Reserva um caminho de arquivo único dentro de ``pasta`` de forma segura
    para uso por múltiplas threads, evitando sobrescrita acidental.
    """
    extensao = extensao if extensao.startswith(".") else f".{extensao}"
    with lock:
        candidato = base
        contador = 2
        caminho = pasta / f"{candidato}{extensao}"
        chave = str(caminho).lower()
        while chave in reservados or caminho.exists():
            candidato = f"{base}_{contador}"
            caminho = pasta / f"{candidato}{extensao}"
            chave = str(caminho).lower()
            contador += 1
        reservados.add(chave)
        return caminho


# =============================================================================
# Processamento 1: grade de pontos
# =============================================================================

def criar_grade_pontos(
    shp_poligono: str,
    pasta_saida: str,
    distancias_usuario: List[float],
    unidade: str,
    limite_estimado_pontos: int = 2_000_000,
    log_callback=None,
    load_callback=None,
    progress_callback=None,
    stop_event: Optional[threading.Event] = None,
) -> List[Path]:
    """Cria uma grade de pontos uniforme dentro da geometria de polígonos."""
    def log(msg: str):
        if log_callback:
            log_callback(msg)

    def set_load(acao: str, mensagem: str = ""):
        if load_callback:
            load_callback(acao, mensagem)

    def set_progress(valor: float, maximo: float, mensagem: str = ""):
        if progress_callback:
            progress_callback(valor, maximo, mensagem)

    checar_parada(stop_event)
    set_load("start", f"Carregando shapefile de polígonos: {Path(shp_poligono).name}")
    try:
        gdf_poly = ler_shapefile(shp_poligono)
    finally:
        set_load("stop", "Carregamento finalizado")
    checar_parada(stop_event)

    gdf_poly = corrigir_geometrias(gdf_poly)
    tipos = set(gdf_poly.geom_type.unique())
    if not any("Polygon" in t for t in tipos):
        raise ValueError("A aba Grade de Pontos exige shapefile com geometrias de polígonos.")

    crs_original = gdf_poly.crs
    crs_metrico = obter_crs_metrico(gdf_poly)

    log(f"CRS original: {crs_original}")
    log(f"CRS métrico usado no processamento: {crs_metrico.to_string()}")

    gdf_m = gdf_poly.to_crs(crs_metrico)
    geometria_base = unir_geometrias(gdf_m)
    geometria_preparada = prep(geometria_base)

    minx, miny, maxx, maxy = geometria_base.bounds
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    distancias_m = converter_para_metros(distancias_usuario, unidade)
    arquivos_saida = []
    total_progresso = max(1, len(distancias_m) * 100)
    set_progress(0, total_progresso, "Aguardando início do processamento da grade")

    for idx_dist, (valor_user, distancia_m) in enumerate(zip(distancias_usuario, distancias_m)):
        checar_parada(stop_event)
        rotulo = rotulo_valor(valor_user, unidade)
        base_prog = idx_dist * 100
        log("-" * 70)
        log(f"Gerando grade com espaçamento de {valor_user} {unidade} ({distancia_m:.3f} m).")
        set_progress(base_prog, total_progresso, f"Preparando grade {idx_dist + 1}/{len(distancias_m)}: {valor_user} {unidade}")

        xs = np.arange(minx, maxx + distancia_m, distancia_m)
        ys = np.arange(miny, maxy + distancia_m, distancia_m)
        estimativa = len(xs) * len(ys)
        log(f"Estimativa inicial de pontos no retângulo envolvente: {estimativa:,}".replace(",", "."))

        if estimativa > limite_estimado_pontos:
            raise MemoryError(
                f"A distância {valor_user} {unidade} pode gerar {estimativa:,} pontos antes do recorte. "
                f"Aumente o espaçamento ou eleve o limite no script."
            )

        pontos = []
        ids = []
        contador = 1
        total_linhas = max(1, len(ys))
        passo_log = max(1, total_linhas // 50)

        for idx_y, y in enumerate(ys, start=1):
            checar_parada(stop_event)
            for x in xs:
                p = Point(float(x), float(y))
                try:
                    dentro = geometria_preparada.covers(p)
                except AttributeError:
                    dentro = geometria_base.covers(p)
                if dentro:
                    pontos.append(p)
                    ids.append(contador)
                    contador += 1

            if idx_y % passo_log == 0 or idx_y == total_linhas:
                parcial = min(85.0, (idx_y / total_linhas) * 85.0)
                set_progress(
                    base_prog + parcial,
                    total_progresso,
                    f"Grade {idx_dist + 1}/{len(distancias_m)}: linha {idx_y}/{total_linhas}; pontos válidos: {len(pontos):,}".replace(",", "."),
                )

        if not pontos:
            log(f"Nenhum ponto gerado para {valor_user} {unidade}.")
            set_progress(base_prog + 100, total_progresso, f"Grade {idx_dist + 1}/{len(distancias_m)} concluída sem pontos")
            continue

        gdf_pts = gpd.GeoDataFrame(
            {
                "id_pt": ids,
                "dist_m": [float(distancia_m)] * len(pontos),
                "dist_orig": [str(valor_user)] * len(pontos),
                "unit": [unidade] * len(pontos),
            },
            geometry=pontos,
            crs=crs_metrico,
        )

        gdf_pts = gdf_pts.to_crs(crs_original)
        gdf_pts = filtrar_colunas_shapefile(gdf_pts)

        nome_poly = normalizar_nome(Path(shp_poligono).stem)
        nome_base = f"grade_pontos_{nome_poly}_{rotulo}_{timestamp}"
        set_progress(base_prog + 92, total_progresso, f"Salvando shapefile e ZIP da grade {idx_dist + 1}/{len(distancias_m)}")
        checar_parada(stop_event)
        zip_path = salvar_gdf_como_shapefile_zip(gdf_pts, pasta_saida, nome_base)
        set_progress(base_prog + 100, total_progresso, f"Grade {idx_dist + 1}/{len(distancias_m)} concluída")

        log(f"Pontos gerados: {len(gdf_pts):,}".replace(",", "."))
        log(f"Arquivo salvo: {zip_path}")
        arquivos_saida.append(zip_path)

    return arquivos_saida


# =============================================================================
# Processamento 2: buffers em pontos
# =============================================================================

def gerar_buffers_pontos(
    shp_pontos_lista: List[str],
    pasta_saida: str,
    diametros_usuario: List[float],
    unidade: str,
    log_callback=None,
    load_callback=None,
    progress_callback=None,
    stop_event: Optional[threading.Event] = None,
) -> List[Path]:
    """Gera buffers a partir de um ou mais shapefiles de pontos."""
    def log(msg: str):
        if log_callback:
            log_callback(msg)

    def set_load(acao: str, mensagem: str = ""):
        if load_callback:
            load_callback(acao, mensagem)

    def set_progress(valor: float, maximo: float, mensagem: str = ""):
        if progress_callback:
            progress_callback(valor, maximo, mensagem)

    if not shp_pontos_lista:
        raise ValueError("Selecione pelo menos um shapefile de pontos.")

    diametros_m = converter_para_metros(diametros_usuario, unidade)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    arquivos_saida = []
    total_tarefas = max(1, len(shp_pontos_lista) * len(diametros_m))
    total_progresso = total_tarefas * 100
    tarefa_idx = 0
    set_progress(0, total_progresso, "Aguardando início do processamento dos buffers")

    for idx_arquivo, shp_pontos in enumerate(shp_pontos_lista, start=1):
        checar_parada(stop_event)
        log("-" * 70)
        log(f"Lendo shapefile de pontos: {shp_pontos}")

        set_load("start", f"Carregando shapefile de pontos {idx_arquivo}/{len(shp_pontos_lista)}: {Path(shp_pontos).name}")
        try:
            gdf_pts = ler_shapefile(shp_pontos)
        finally:
            set_load("stop", f"Carregamento finalizado: {Path(shp_pontos).name}")
        checar_parada(stop_event)

        crs_original = gdf_pts.crs
        crs_metrico = obter_crs_metrico(gdf_pts)

        tipos = set(gdf_pts.geom_type.unique())
        if not any("Point" in t for t in tipos):
            raise ValueError(
                f"O arquivo {Path(shp_pontos).name} não parece ser um shapefile de pontos. "
                f"Tipos encontrados: {', '.join(sorted(tipos))}"
            )

        log(f"CRS original: {crs_original}")
        log(f"CRS métrico usado no processamento: {crs_metrico.to_string()}")
        gdf_m = gdf_pts.to_crs(crs_metrico)

        for idx_diam, (valor_user, diametro_m) in enumerate(zip(diametros_usuario, diametros_m), start=1):
            checar_parada(stop_event)
            raio_m = diametro_m / 2.0
            rotulo = rotulo_valor(valor_user, unidade)
            base_prog = tarefa_idx * 100

            log(f"Gerando buffer com DIÂMETRO de {valor_user} {unidade}; raio usado = {raio_m:.3f} m.")
            set_progress(
                base_prog + 10,
                total_progresso,
                f"Arquivo {idx_arquivo}/{len(shp_pontos_lista)} | buffer {idx_diam}/{len(diametros_m)}: preparando {valor_user} {unidade}",
            )

            gdf_buf = gdf_m.copy()
            gdf_buf["diam_m"] = float(diametro_m)
            gdf_buf["radius_m"] = float(raio_m)
            gdf_buf["diam_orig"] = str(valor_user)
            gdf_buf["unit"] = unidade

            set_progress(base_prog + 35, total_progresso, f"Calculando geometrias de buffer: {valor_user} {unidade}")
            gdf_buf["geometry"] = gdf_buf.geometry.buffer(raio_m)
            gdf_buf = gdf_buf.set_geometry("geometry")
            checar_parada(stop_event)

            set_progress(base_prog + 60, total_progresso, "Reprojetando resultado para o CRS original")
            gdf_buf = gdf_buf.to_crs(crs_original)
            gdf_buf = filtrar_colunas_shapefile(gdf_buf)

            nome_pts = normalizar_nome(Path(shp_pontos).stem)
            nome_base = f"buffer_{nome_pts}_diam_{rotulo}_{timestamp}"
            set_progress(base_prog + 85, total_progresso, f"Salvando shapefile e ZIP: {nome_base}")
            checar_parada(stop_event)
            zip_path = salvar_gdf_como_shapefile_zip(gdf_buf, pasta_saida, nome_base)
            set_progress(base_prog + 100, total_progresso, f"Buffer concluído: {valor_user} {unidade}")
            tarefa_idx += 1

            log(f"Buffers gerados: {len(gdf_buf):,}".replace(",", "."))
            log(f"Arquivo salvo: {zip_path}")
            arquivos_saida.append(zip_path)

    return arquivos_saida


# =============================================================================
# Processamento 3: recorte de rasters por shapefile de buffers
# =============================================================================

def salvar_array_reprojetado_utm(
    array_clipado: np.ndarray,
    transform_clipado,
    src_dataset,
    crs_destino: CRS,
    caminho_saida: Path,
    log_callback=None,
    stop_event: Optional[threading.Event] = None,
    num_threads: int = 1,
):
    """Reprojeta o raster já recortado para UTM e salva em GeoTIFF."""
    def log(msg: str):
        if log_callback:
            log_callback(msg)

    altura = array_clipado.shape[1]
    largura = array_clipado.shape[2]
    bounds = array_bounds(altura, largura, transform_clipado)

    dst_transform, dst_width, dst_height = calculate_default_transform(
        src_dataset.crs,
        crs_destino,
        largura,
        altura,
        *bounds,
    )

    nodata_saida = src_dataset.nodata
    if nodata_saida is None:
        nodata_saida = 0

    destino = np.full(
        (src_dataset.count, dst_height, dst_width),
        nodata_saida,
        dtype=array_clipado.dtype,
    )

    num_threads = max(1, int(num_threads))
    for banda in range(src_dataset.count):
        checar_parada(stop_event)
        reproject(
            source=array_clipado[banda],
            destination=destino[banda],
            src_transform=transform_clipado,
            src_crs=src_dataset.crs,
            src_nodata=nodata_saida,
            dst_transform=dst_transform,
            dst_crs=crs_destino,
            dst_nodata=nodata_saida,
            resampling=Resampling.nearest,
            num_threads=num_threads,
        )

    perfil = src_dataset.profile.copy()
    perfil.update(
        driver="GTiff",
        height=dst_height,
        width=dst_width,
        transform=dst_transform,
        crs=crs_destino,
        nodata=nodata_saida,
        compress="lzw",
        tiled=True,
        BIGTIFF="IF_SAFER",
    )

    caminho_saida.parent.mkdir(parents=True, exist_ok=True)
    checar_parada(stop_event)
    with rasterio.open(caminho_saida, "w", **perfil) as dst:
        dst.write(destino)
    remover_auxiliares_raster(caminho_saida)

    log(f"Raster salvo em UTM: {caminho_saida}")


def recortar_rasters_por_buffers(
    shp_buffers: str,
    rasters_lista: List[str],
    pasta_saida: str,
    campo_id: str = "id_pt",
    campo_diametro: str = "diam_m",
    all_touched: bool = True,
    max_workers: Optional[int] = None,
    log_callback=None,
    load_callback=None,
    progress_callback=None,
    stop_event: Optional[threading.Event] = None,
) -> List[Path]:
    """
    Recorta cada raster por cada feição do shapefile de buffers e salva em UTM.

    Otimizações da V6:
    - As projeções UTM de cada feição são calculadas uma única vez (não por raster).
    - As geometrias do buffer reprojetadas para o CRS do raster são reaproveitadas
      entre rasters que compartilham o mesmo CRS.
    - O recorte de cada feição roda em paralelo (ThreadPoolExecutor), pois o GDAL
      libera o GIL durante recorte/reprojeção.
    - Cada raster é salvo em uma subpasta com o ANO extraído do seu nome.
    - Nome padronizado: clip_id_<ID>_diam_<DIAMETRO>_<ANO>.tif
    """
    def log(msg: str):
        if log_callback:
            log_callback(msg)

    def set_load(acao: str, mensagem: str = ""):
        if load_callback:
            load_callback(acao, mensagem)

    def set_progress(valor: float, maximo: float, mensagem: str = ""):
        if progress_callback:
            progress_callback(valor, maximo, mensagem)

    if not shp_buffers:
        raise ValueError("Selecione o shapefile de buffers.")
    if not rasters_lista:
        raise ValueError("Selecione pelo menos um arquivo raster.")
    if not pasta_saida:
        raise ValueError("Selecione a pasta de saída dos rasters recortados.")

    if not max_workers or max_workers < 1:
        max_workers = numero_padrao_de_nucleos()

    checar_parada(stop_event)
    set_load("start", f"Carregando shapefile de buffers: {Path(shp_buffers).name}")
    try:
        gdf_buffers = ler_shapefile(shp_buffers)
    finally:
        set_load("stop", "Carregamento do shapefile de buffers finalizado")
    checar_parada(stop_event)

    gdf_buffers = corrigir_geometrias(gdf_buffers)
    tipos = set(gdf_buffers.geom_type.unique())
    if not any("Polygon" in t for t in tipos):
        raise ValueError(
            "A aba de recorte raster exige shapefile com geometrias de polígonos/buffers. "
            f"Tipos encontrados: {', '.join(sorted(tipos))}"
        )

    if campo_id not in gdf_buffers.columns:
        raise ValueError(f"O campo '{campo_id}' não existe na tabela de atributos do shapefile de buffers.")
    if campo_diametro not in gdf_buffers.columns:
        raise ValueError(f"O campo '{campo_diametro}' não existe na tabela de atributos do shapefile de buffers.")

    log(f"Feições de buffer carregadas: {len(gdf_buffers):,}".replace(",", "."))
    log(f"Campo de identificação usado: {campo_id}")
    log(f"Campo de diâmetro usado: {campo_diametro}")
    log(f"CRS do shapefile de buffers: {gdf_buffers.crs}")
    log(f"Núcleos paralelos usados no recorte: {max_workers}")

    # Pré-calcula, de uma só vez, a projeção UTM e os rótulos id/diâmetro de cada feição.
    set_load("start", "Calculando projeções UTM e rótulos das feições")
    try:
        utm_por_feature = calcular_utm_por_feature(gdf_buffers)
        rotulos_por_feature: Dict[object, tuple] = {}
        for idx_feature in gdf_buffers.index:
            row = gdf_buffers.loc[idx_feature]
            rotulos_por_feature[idx_feature] = (
                formatar_valor_atributo(row[campo_id]),
                formatar_valor_atributo(row[campo_diametro]),
            )
    finally:
        set_load("stop", "Projeções UTM e rótulos calculados")
    checar_parada(stop_event)

    rasters_validos = []
    for idx, raster_path in enumerate(rasters_lista, start=1):
        checar_parada(stop_event)
        set_load("start", f"Verificando raster {idx}/{len(rasters_lista)}: {Path(raster_path).name}")
        try:
            with rasterio.open(raster_path) as src:
                if src.crs is None:
                    raise ValueError(f"O raster não possui CRS definido: {Path(raster_path).name}")
                log(
                    f"Raster {idx}/{len(rasters_lista)}: {Path(raster_path).name} | "
                    f"bandas={src.count}, tamanho={src.width}x{src.height}, CRS={src.crs}"
                )
                rasters_validos.append((raster_path, str(src.crs)))
        finally:
            set_load("stop", f"Raster verificado: {Path(raster_path).name}")

    total_tarefas = max(1, len(rasters_validos) * len(gdf_buffers))
    total_progresso = total_tarefas * 100
    pasta_saida = Path(pasta_saida)
    arquivos_saida: List[Path] = []
    set_progress(0, total_progresso, "Aguardando início dos recortes raster")

    # Distribui as threads do GDAL (warp) entre os workers para não saturar a CPU.
    cpus = os.cpu_count() or 4
    reproject_threads = max(1, cpus // max(1, max_workers))

    # Estado compartilhado entre threads.
    contador_lock = threading.Lock()
    concluidas = {"n": 0}
    nomes_reservados: set = set()
    nomes_lock = threading.Lock()

    # Cache de geometrias do buffer já reprojetadas, por CRS do raster.
    geoms_por_crs: Dict[str, Dict[object, object]] = {}

    def obter_geoms_para_crs(crs_str: str) -> Dict[object, object]:
        if crs_str not in geoms_por_crs:
            gdf_no_crs = gdf_buffers.to_crs(crs_str)
            geoms_por_crs[crs_str] = {
                idx: gdf_no_crs.geometry.loc[idx] for idx in gdf_buffers.index
            }
        return geoms_por_crs[crs_str]

    def avancar(mensagem: str):
        with contador_lock:
            concluidas["n"] += 1
            n = concluidas["n"]
        set_progress(n * 100, total_progresso, mensagem)

    for idx_raster, (raster_path, crs_str) in enumerate(rasters_validos, start=1):
        checar_parada(stop_event)
        ano = extrair_ano_do_nome(raster_path)
        nome_raster = Path(raster_path).name
        log("-" * 70)
        log(f"Processando raster {idx_raster}/{len(rasters_validos)}: {nome_raster}")
        log(f"Ano extraído do nome do raster: {ano}")

        # Subpasta com o ano do raster.
        pasta_ano = pasta_saida / normalizar_nome(ano, max_len=20)
        pasta_ano.mkdir(parents=True, exist_ok=True)

        geoms_feature = obter_geoms_para_crs(crs_str)

        def processar_feature(idx_feature):
            checar_parada(stop_event)
            id_val, diam_val = rotulos_por_feature[idx_feature]
            geom_raster = geoms_feature.get(idx_feature)

            if geom_raster is None or geom_raster.is_empty:
                avancar(f"Feição ignorada (vazia): id={id_val}, diam={diam_val}")
                return ("vazia", id_val, diam_val, None)

            try:
                with rasterio.open(raster_path) as src:
                    array_clipado, transform_clipado = mask(
                        src,
                        [mapping(geom_raster)],
                        crop=True,
                        all_touched=all_touched,
                        filled=True,
                    )
                    if array_clipado.size == 0:
                        avancar(f"Recorte vazio: id={id_val}, diam={diam_val}")
                        return ("recorte_vazio", id_val, diam_val, None)

                    crs_utm = utm_por_feature.get(idx_feature)
                    if crs_utm is None:
                        crs_utm = obter_crs_utm_por_geometria(gdf_buffers.loc[[idx_feature]])
                    epsg_utm = crs_utm.to_epsg() or "utm"

                    base_nome = normalizar_nome(
                        f"clip_id_{id_val}_diam_{diam_val}_{ano}", max_len=150
                    )
                    caminho_saida = reservar_caminho_unico(
                        pasta_ano, base_nome, ".tif", nomes_reservados, nomes_lock
                    )

                    salvar_array_reprojetado_utm(
                        array_clipado=array_clipado,
                        transform_clipado=transform_clipado,
                        src_dataset=src,
                        crs_destino=crs_utm,
                        caminho_saida=caminho_saida,
                        log_callback=None,
                        stop_event=stop_event,
                        num_threads=reproject_threads,
                    )
            except ProcessamentoInterrompido:
                raise
            except ValueError:
                avancar(f"Sem sobreposição: id={id_val}, diam={diam_val}")
                return ("sem_sobreposicao", id_val, diam_val, None)

            avancar(f"Recorte salvo: {caminho_saida.name} (EPSG:{epsg_utm})")
            return ("ok", id_val, diam_val, caminho_saida)

        indices = list(gdf_buffers.index)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futuros = {executor.submit(processar_feature, idx): idx for idx in indices}
            try:
                for futuro in as_completed(futuros):
                    if stop_event is not None and stop_event.is_set():
                        for f in futuros:
                            f.cancel()
                        raise ProcessamentoInterrompido("Processamento interrompido pelo usuário.")
                    status, id_val, diam_val, caminho = futuro.result()
                    if status == "ok" and caminho is not None:
                        arquivos_saida.append(caminho)
                        log(f"  -> [{ano}] id={id_val} diam={diam_val} salvo: {caminho.name}")
                    elif status == "sem_sobreposicao":
                        log(f"  -> Sem sobreposição: raster={nome_raster}, id={id_val}, diam={diam_val}")
                    elif status in ("vazia", "recorte_vazio"):
                        log(f"  -> Recorte ignorado ({status}): id={id_val}, diam={diam_val}")
            except ProcessamentoInterrompido:
                for f in futuros:
                    f.cancel()
                raise

        log(f"Raster {nome_raster} concluído. Recortes salvos em: {pasta_ano}")

    return arquivos_saida


# =============================================================================
# Processamento 4: métricas de paisagem com PyLandStats
# =============================================================================

def listar_rasters_pasta(pasta: str, recursivo: bool = False) -> List[Path]:
    pasta = Path(pasta)
    if not pasta.exists() or not pasta.is_dir():
        raise ValueError("Selecione uma pasta válida contendo os rasters recortados.")
    padroes = ["*.tif", "*.tiff", "*.TIF", "*.TIFF"]
    arquivos: List[Path] = []
    for padrao in padroes:
        if recursivo:
            arquivos.extend(pasta.rglob(padrao))
        else:
            arquivos.extend(pasta.glob(padrao))
    # Remove duplicados e ordena.
    arquivos = sorted(set(arquivos), key=lambda p: str(p).lower())
    # Ignora a subpasta em que a aba 6 grava os rasters reprojetados para metros
    # (quando ela fica dentro da pasta escolhida), para não ler o mesmo raster duas vezes.
    arquivos = [p for p in arquivos if PASTA_REPROJ_METRICAS not in p.relative_to(pasta).parts[:-1]]
    return arquivos


# Subpasta (dentro da pasta de saída dos CSVs) onde a aba 6 grava os rasters
# reprojetados para metros antes do cálculo das métricas.
PASTA_REPROJ_METRICAS = "rasters_reprojetados_metros"


def _info_projecao_raster(caminho: str):
    """Retorna (situacao, descricao) do CRS de um raster.

    situacao: "metros" (projetado em metros), "nao_metros" (graus ou outra
    unidade), "sem_crs" (sem CRS definido) ou "erro" (não foi possível abrir).
    """
    try:
        with rasterio.open(caminho) as src:
            crs = src.crs
            if crs is None:
                return "sem_crs", "sem CRS definido"
            try:
                desc = crs.to_string()
            except Exception:
                desc = str(crs)
            if len(desc) > 60:
                desc = desc[:57] + "..."
            if _crs_e_metrico(crs):
                return "metros", desc
            try:
                unidade = "graus" if crs.is_geographic else (crs.linear_units or "desconhecida")
            except Exception:
                unidade = "desconhecida"
            return "nao_metros", f"{desc}, unidade: {unidade}"
    except Exception as exc:
        return "erro", str(exc)


def verificar_projecao_rasters(rasters, progress_callback=None, stop_event=None) -> Dict[str, list]:
    """Verifica, raster a raster, se o CRS é projetado em metros.

    Retorna um dicionário com as chaves "metros", "nao_metros", "sem_crs" e
    "erro"; cada valor é uma lista de (caminho, descrição do CRS ou do erro).
    """
    caminhos = [str(r) for r in rasters]
    resultado: Dict[str, list] = {"metros": [], "nao_metros": [], "sem_crs": [], "erro": []}
    total = len(caminhos)
    executor = ThreadPoolExecutor(max_workers=8)
    try:
        for i, (caminho, (situacao, desc)) in enumerate(
                zip(caminhos, executor.map(_info_projecao_raster, caminhos)), start=1):
            checar_parada(stop_event)
            resultado[situacao].append((caminho, desc))
            if progress_callback and (i % 25 == 0 or i == total):
                progress_callback(i, max(1, total), f"Projeção verificada: {i}/{total}")
    finally:
        executor.shutdown(wait=False, cancel_futures=True)
    return resultado


def registrar_verificacao_projecao(resultado: Dict[str, list], pasta_base: str, log) -> str:
    """Escreve no log o resultado de verificar_projecao_rasters (um aviso por
    raster fora de metros) e devolve o resumo em uma linha."""
    def rel(caminho: str) -> str:
        try:
            return os.path.relpath(caminho, pasta_base)
        except ValueError:
            return caminho

    for chave, rotulo in (("nao_metros", "NÃO está em metros"), ("sem_crs", "SEM CRS definido"),
                          ("erro", "ERRO ao abrir")):
        linhas = [f"{rotulo}: {rel(c)} ({d})" for c, d in resultado[chave]]
        for i in range(0, len(linhas), 200):  # em blocos, para não travar a janela com milhares de linhas
            log("\n".join(linhas[i:i + 200]))
    total = sum(len(v) for v in resultado.values())
    resumo = (f"Rasters verificados: {total} | em metros: {len(resultado['metros'])} | "
              f"NÃO em metros: {len(resultado['nao_metros'])} | sem CRS: {len(resultado['sem_crs'])} | "
              f"erros de leitura: {len(resultado['erro'])}")
    log(resumo)
    return resumo


def preparar_rasters_em_metros(
    pasta_rasters: str,
    pasta_saida: str,
    recursivo: bool = True,
    max_workers: Optional[int] = None,
    log_callback=None,
    load_callback=None,
    progress_callback=None,
    stop_event: Optional[threading.Event] = None,
) -> List[str]:
    """
    Garante que os rasters usados nas métricas estejam em metros.

    Verifica o CRS de cada raster da pasta; os que não estiverem em metros são
    reprojetados (vizinho mais próximo) para o SIRGAS 2000 / UTM da zona do seu
    centróide e gravados em <pasta_saida>/rasters_reprojetados_metros, com o
    mesmo nome e a mesma estrutura de subpastas. Devolve a lista de rasters a
    usar no cálculo: os originais que já estavam em metros e as cópias
    reprojetadas dos demais. Rasters sem CRS, ilegíveis ou cuja reprojeção
    falhou ficam fora da lista e são informados no log.
    """
    def log(msg: str):
        if log_callback:
            log_callback(msg)

    if load_callback:
        load_callback("start", f"Listando rasters em: {pasta_rasters}")
    try:
        rasters = listar_rasters_pasta(pasta_rasters, recursivo=recursivo)
    finally:
        if load_callback:
            load_callback("stop", "Listagem de rasters concluída")
    if not rasters:
        raise ValueError("Nenhum raster .tif/.tiff foi encontrado na pasta selecionada.")

    log(f"Verificando se cada raster está em metros ({len(rasters)} raster(s))...")
    resultado = verificar_projecao_rasters(rasters, progress_callback=progress_callback, stop_event=stop_event)
    registrar_verificacao_projecao(resultado, pasta_rasters, log)

    finais: Dict[str, str] = {c: c for c, _ in resultado["metros"]}
    nao_metros = [c for c, _ in resultado["nao_metros"]]
    if nao_metros:
        destino = Path(pasta_saida) / PASTA_REPROJ_METRICAS
        log(f"Reprojetando {len(nao_metros)} raster(s) para metros (UTM) em: {destino}")
        grupos: Dict[int, List[str]] = {}
        for caminho in nao_metros:
            checar_parada(stop_event)
            epsg = sugerir_epsg_de_raster(caminho)
            if epsg is None:
                log(f"ERRO: não foi possível definir o EPSG UTM de {Path(caminho).name} (raster ignorado no cálculo).")
                continue
            grupos.setdefault(epsg, []).append(caminho)
        for epsg, lista in sorted(grupos.items()):
            res = reprojetar_rasters_lote(
                rasters_lista=lista,
                pasta_saida=str(destino),
                epsg_destino=epsg,
                sufixo_epsg=False,
                resampling=Resampling.nearest,
                max_workers=max_workers,
                log_callback=log_callback,
                progress_callback=progress_callback,
                stop_event=stop_event,
                pasta_base=pasta_rasters,
                normalizar_nomes=False,
            )
            finais.update(res["mapa"])

    usados = [finais[str(r)] for r in rasters if str(r) in finais]
    reprojetados = len(usados) - len(resultado["metros"])
    log(f"Rasters prontos para o cálculo: {len(usados)} | já em metros: {len(resultado['metros'])} | "
        f"reprojetados: {reprojetados} | ignorados: {len(rasters) - len(usados)}")
    if not usados:
        raise ValueError("Nenhum raster em metros disponível para o cálculo (veja o log).")
    return usados


def preparar_df_metricas(df: pd.DataFrame, nivel: str, raster_path: Path) -> pd.DataFrame:
    """Padroniza DataFrame de métricas para exportação."""
    if df is None:
        return pd.DataFrame()
    if isinstance(df, pd.Series):
        df = df.to_frame().T
    df = df.copy()
    df = df.reset_index()

    # Renomeia coluna index genérica quando existir.
    if "index" in df.columns:
        if nivel == "patch":
            df = df.rename(columns={"index": "patch_id"})
        elif nivel == "class":
            df = df.rename(columns={"index": "class_val"})
        else:
            df = df.rename(columns={"index": "row_id"})

    # Garante metadados exigidos pelo usuário e úteis para rastreabilidade.
    df.insert(0, "level", nivel)
    df.insert(0, "year", extrair_ano_de_recorte(str(raster_path)))
    df.insert(0, "source_stem", raster_path.stem)
    df.insert(0, "source_file", raster_path.name)
    df.insert(0, "source_path", str(raster_path))
    return df


def importar_pylandstats():
    try:
        import pylandstats as pls
    except ImportError as exc:
        raise ImportError(
            "O pacote PyLandStats não está instalado neste ambiente. Instale com: pip install pylandstats"
        ) from exc
    return pls


def calcular_metricas_paisagem_pylandstats(
    pasta_rasters: str,
    pasta_saida: str,
    niveis: Dict[str, bool],
    metricas_por_nivel: Dict[str, Sequence[str]],
    vizinhanca: str = "8",
    recursivo: bool = True,
    max_workers: Optional[int] = None,
    engine_label: str = "pylandstats",
    log_callback=None,
    load_callback=None,
    progress_callback=None,
    stop_event: Optional[threading.Event] = None,
    rasters_lista: Optional[Sequence[str]] = None,
) -> Dict[str, object]:
    """
    Calcula métricas de paisagem para todos os rasters de uma pasta (e subpastas).
    Salva um CSV individual por raster (em subpasta do ANO) e um CSV geral consolidado.

    Otimizações da V6:
    - Cada raster é processado em paralelo (ThreadPoolExecutor); a leitura do raster
      e a rotulagem de manchas do PyLandStats (scipy.ndimage) liberam o GIL.
    - Validação de rasters mais leve (sem ler a banda inteira).
    """
    def log(msg: str):
        if log_callback:
            log_callback(msg)

    def set_load(acao: str, mensagem: str = ""):
        if load_callback:
            load_callback(acao, mensagem)

    def set_progress(valor: float, maximo: float, mensagem: str = ""):
        if progress_callback:
            progress_callback(valor, maximo, mensagem)

    if vizinhanca not in {"4", "8"}:
        raise ValueError("A vizinhança deve ser '4' ou '8'.")
    if not any(niveis.values()):
        raise ValueError("Selecione pelo menos um nível escalar: patch, classe ou paisagem.")

    pasta_saida_path = Path(pasta_saida)
    if not pasta_saida_path:
        raise ValueError("Selecione a pasta de saída dos CSVs.")
    pasta_saida_path.mkdir(parents=True, exist_ok=True)

    checar_parada(stop_event)
    set_load("start", "Importando PyLandStats")
    pls = importar_pylandstats()
    set_load("stop", "PyLandStats importado")
    checar_parada(stop_event)

    set_load("start", f"Listando rasters em: {pasta_rasters}")
    if rasters_lista is not None:  # lista já preparada (ex.: rasters reprojetados para metros)
        rasters = [Path(p) for p in rasters_lista]
    else:
        rasters = listar_rasters_pasta(pasta_rasters, recursivo=recursivo)
    set_load("stop", f"Rasters encontrados: {len(rasters)}")

    if not rasters:
        raise ValueError("Nenhum raster .tif/.tiff foi encontrado na pasta selecionada.")

    if not max_workers or max_workers < 1:
        max_workers = numero_padrao_de_nucleos()

    # Validação rápida com rasterio (sem ler a banda inteira, para ganhar velocidade).
    rasters_validos: List[Path] = []
    for idx, raster_path in enumerate(rasters, start=1):
        checar_parada(stop_event)
        set_load("start", f"Verificando raster {idx}/{len(rasters)}: {raster_path.name}")
        try:
            with rasterio.open(raster_path) as src:
                if src.count < 1:
                    raise ValueError(f"Raster sem bandas: {raster_path.name}")
                if src.crs is None:
                    log(f"Aviso: raster sem CRS definido, mas as métricas categóricas ainda podem ser calculadas: {raster_path.name}")
                log(f"Raster validado: {raster_path.name} | {src.width}x{src.height} | bandas={src.count} | dtype={src.dtypes[0]}")
                rasters_validos.append(raster_path)
        except Exception as exc:
            log(f"Raster ignorado (inválido): {raster_path.name} -> {exc}")
        finally:
            set_load("stop", f"Raster verificado: {raster_path.name}")

    if not rasters_validos:
        raise ValueError("Nenhum raster válido foi encontrado para o cálculo de métricas.")

    niveis_ativos = [nivel for nivel, ativo in niveis.items() if ativo]
    total_progresso = max(1, len(rasters_validos)) * 100
    set_progress(0, total_progresso, "Aguardando início do cálculo das métricas")
    log(f"Núcleos paralelos usados nas métricas: {max_workers}")

    arquivos_individuais: List[Path] = []
    dfs_gerais: List[pd.DataFrame] = []
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    contador_lock = threading.Lock()
    concluidas = {"n": 0}

    def computar_raster(raster_path: Path):
        """Calcula todos os níveis selecionados de um raster e devolve o DataFrame."""
        checar_parada(stop_event)
        mensagens: List[str] = []
        try:
            ls = pls.Landscape(str(raster_path), neighborhood_rule=str(vizinhanca))
        except Exception as exc:
            return raster_path, None, [f"ERRO ao criar Landscape para {raster_path.name}: {exc}"]

        dfs_raster: List[pd.DataFrame] = []
        for nivel in niveis_ativos:
            checar_parada(stop_event)
            metricas = list(metricas_por_nivel.get(nivel, []))
            metricas_arg = metricas if metricas else None
            try:
                if nivel == "patch":
                    df = ls.compute_patch_metrics_df(metrics=metricas_arg)
                elif nivel == "class":
                    df = ls.compute_class_metrics_df(metrics=metricas_arg)
                elif nivel == "landscape":
                    df = ls.compute_landscape_metrics_df(metrics=metricas_arg)
                else:
                    raise ValueError(f"Nível desconhecido: {nivel}")

                df = preparar_df_metricas(df, nivel=nivel, raster_path=raster_path)
                df["neighborhood_rule"] = str(vizinhanca)
                df["engine"] = engine_label
                dfs_raster.append(df)
                mensagens.append(f"{raster_path.name} | nível {nivel}: {df.shape[0]} linhas x {df.shape[1]} colunas")
            except Exception as exc:
                mensagens.append(f"ERRO no nível {nivel} para {raster_path.name}: {exc}")

        if dfs_raster:
            df_final = pd.concat(dfs_raster, ignore_index=True, sort=False)
        else:
            df_final = None
        return raster_path, df_final, mensagens

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futuros = {executor.submit(computar_raster, rp): rp for rp in rasters_validos}
        try:
            for futuro in as_completed(futuros):
                if stop_event is not None and stop_event.is_set():
                    for f in futuros:
                        f.cancel()
                    raise ProcessamentoInterrompido("Processamento interrompido pelo usuário.")

                raster_path, df_final, mensagens = futuro.result()
                for m in mensagens:
                    log(m)

                if df_final is not None and not df_final.empty:
                    # Subpasta do ANO específico, lido do PADRÃO do nome do recorte
                    # (clip_id_<ID>_diam_<DIAM>_<ANO>): o ano é o último campo do nome.
                    ano = extrair_ano_de_recorte(str(raster_path))
                    pasta_ano = pasta_saida_path / normalizar_nome(ano, max_len=20)
                    pasta_ano.mkdir(parents=True, exist_ok=True)

                    nome_csv = normalizar_nome(f"metricas_paisagem_{engine_label}_{raster_path.stem}", max_len=160) + ".csv"
                    caminho_csv = pasta_ano / nome_csv
                    df_final.to_csv(caminho_csv, index=False, encoding="utf-8-sig")
                    arquivos_individuais.append(caminho_csv)
                    dfs_gerais.append(df_final)
                    log(f"CSV individual salvo: {caminho_csv}")
                else:
                    log(f"Nenhuma métrica foi gerada para: {raster_path.name}")

                with contador_lock:
                    concluidas["n"] += 1
                    n = concluidas["n"]
                set_progress(n * 100, total_progresso, f"Rasters concluídos: {n}/{len(rasters_validos)}")
        except ProcessamentoInterrompido:
            for f in futuros:
                f.cancel()
            raise

    if dfs_gerais:
        df_geral = pd.concat(dfs_gerais, ignore_index=True, sort=False)
        caminho_geral = pasta_saida_path / f"metricas_paisagem_GERAL_{engine_label}_{timestamp}.csv"
        df_geral.to_csv(caminho_geral, index=False, encoding="utf-8-sig")
        log(f"CSV geral consolidado salvo: {caminho_geral}")
    else:
        caminho_geral = None

    return {
        "rasters_processados": len(rasters_validos),
        "csvs_individuais": arquivos_individuais,
        "csv_geral": caminho_geral,
    }


# =============================================================================
# Processamento 4-R: métricas de paisagem com o pacote R `landscapemetrics`
# (chamado a partir do Python via rpy2). Mantém o mesmo fluxo de salvamento:
# CSV individual por raster em subpasta do ANO + CSV geral consolidado.
# =============================================================================

# Estado por processo (para o ProcessPoolExecutor reaproveitar a sessão R).
_R_STATE: Dict[str, object] = {}


def _lsm_r_init(cellsize: float = 30.0):
    """Inicializa a sessão R (rpy2) num processo/worker e importa os pacotes."""
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri, StrVector
    from rpy2.robjects.conversion import localconverter

    terra = importr("terra")
    lsm = importr("landscapemetrics")
    try:
        # Opção pedida pelo usuário (default 30). Harmless se a versão ignorar.
        ro.r("options(landscapemetrics_cellsize = %s)" % repr(float(cellsize)))
    except Exception:
        pass
    _R_STATE.update(dict(ro=ro, importr=importr, pandas2ri=pandas2ri,
                         StrVector=StrVector, localconverter=localconverter,
                         lsm=lsm, terra=terra, cellsize=float(cellsize)))


def _lsm_pivot_meta(pdf, nivel: str, raster_path: Path, vizinhanca: str):
    """Converte a saída longa do landscapemetrics em formato largo + metadados,
    no mesmo esquema do PyLandStats (source_*, year, level, ...)."""
    if pdf is None or len(pdf) == 0:
        return None
    idx = [c for c in ["level"] if c in pdf.columns]
    if "class" in pdf.columns:
        idx.append("class")
    if nivel == "patch" and "id" in pdf.columns:
        idx.append("id")
    if not idx:
        idx = None
    try:
        wide = pdf.pivot_table(index=idx, columns="metric", values="value", aggfunc="first").reset_index()
    except Exception:
        wide = pdf.copy()
    try:
        wide.columns.name = None
    except Exception:
        pass
    ren = {}
    if "class" in wide.columns:
        ren["class"] = "class_val"
    if "id" in wide.columns:
        ren["id"] = "patch_id"
    wide = wide.rename(columns=ren)
    if "level" in wide.columns:
        wide = wide.drop(columns=["level"])
    wide.insert(0, "level", nivel)
    wide.insert(0, "year", extrair_ano_de_recorte(str(raster_path)))
    wide.insert(0, "source_stem", raster_path.stem)
    wide.insert(0, "source_file", raster_path.name)
    wide.insert(0, "source_path", str(raster_path))
    wide["neighborhood_rule"] = str(vizinhanca)
    wide["engine"] = "landscapemetrics"
    return wide


def _lsm_r_compute(args):
    """Calcula, no R, todos os níveis selecionados de UM raster. Executa em
    processo separado (ProcessPoolExecutor) ou em sequência."""
    raster_path_str, niveis_ativos, metricas_func_por_nivel, vizinhanca, cellsize = args
    raster_path = Path(raster_path_str)
    mensagens: List[str] = []
    if not _R_STATE:
        try:
            _lsm_r_init(cellsize)
        except Exception as exc:
            return raster_path_str, None, ["ERRO ao inicializar R/landscapemetrics: %s" % exc]

    ro = _R_STATE["ro"]; lsm = _R_STATE["lsm"]; terra = _R_STATE["terra"]
    StrVector = _R_STATE["StrVector"]; pandas2ri = _R_STATE["pandas2ri"]
    localconverter = _R_STATE["localconverter"]

    try:
        r_rast = terra.rast(str(raster_path))
    except Exception as exc:
        return raster_path_str, None, ["ERRO ao abrir raster no R (%s): %s" % (raster_path.name, exc)]

    dfs_raster: List[pd.DataFrame] = []
    for nivel in niveis_ativos:
        funcs = list(metricas_func_por_nivel.get(nivel, []) or [])
        if not funcs:
            continue
        try:
            res = lsm.calculate_lsm(r_rast, what=StrVector(funcs),
                                    directions=int(vizinhanca), verbose=False)
            with localconverter(ro.default_converter + pandas2ri.converter):
                pdf = ro.conversion.rpy2py(res)
            wide = _lsm_pivot_meta(pdf, nivel, raster_path, vizinhanca)
            if wide is not None and len(wide):
                dfs_raster.append(wide)
                mensagens.append("%s | nível %s (R): %d linha(s) x %d coluna(s)"
                                 % (raster_path.name, nivel, wide.shape[0], wide.shape[1]))
            else:
                mensagens.append("%s | nível %s (R): sem resultados" % (raster_path.name, nivel))
        except Exception as exc:
            mensagens.append("ERRO no nível %s (R) para %s: %s" % (nivel, raster_path.name, exc))

    df_final = pd.concat(dfs_raster, ignore_index=True, sort=False) if dfs_raster else None
    return raster_path_str, df_final, mensagens


def calcular_metricas_paisagem_r(
    pasta_rasters: str,
    pasta_saida: str,
    niveis: Dict[str, bool],
    metricas_func_por_nivel: Dict[str, Sequence[str]],
    vizinhanca: str = "8",
    cellsize: float = 30.0,
    recursivo: bool = True,
    max_workers: Optional[int] = None,
    log_callback=None,
    load_callback=None,
    progress_callback=None,
    stop_event: Optional[threading.Event] = None,
    rasters_lista: Optional[Sequence[str]] = None,
) -> Dict[str, object]:
    """Calcula métricas com o pacote R landscapemetrics (via rpy2), para todos
    os rasters da pasta (e subpastas). Salva CSV por raster (subpasta do ANO) e
    um CSV geral. Paraleliza por raster com ProcessPoolExecutor (cada processo
    tem sua própria sessão R); cai para modo sequencial se necessário."""
    def log(msg: str):
        if log_callback:
            log_callback(msg)

    def set_load(acao: str, mensagem: str = ""):
        if load_callback:
            load_callback(acao, mensagem)

    def set_progress(valor: float, maximo: float, mensagem: str = ""):
        if progress_callback:
            progress_callback(valor, maximo, mensagem)

    if vizinhanca not in {"4", "8"}:
        raise ValueError("A vizinhança (directions) deve ser '4' ou '8'.")
    if not any(niveis.values()):
        raise ValueError("Selecione pelo menos um nível escalar: patch, classe ou paisagem.")

    pasta_saida_path = Path(pasta_saida)
    pasta_saida_path.mkdir(parents=True, exist_ok=True)

    # Validação do ambiente R/rpy2 (mensagem clara se faltar).
    checar_parada(stop_event)
    set_load("start", "Verificando rpy2 + R (landscapemetrics/terra)")
    import importlib.util
    if importlib.util.find_spec("rpy2") is None:
        raise ImportError(
            "O pacote Python 'rpy2' não está instalado. Para usar o motor R:\n"
            "  1) Instale o R (https://cran.r-project.org) e, no R, os pacotes:\n"
            "     install.packages(c('landscapemetrics','terra'))\n"
            "  2) No Python: pip install rpy2\n"
            "  3) Garanta que a variável R_HOME aponte para a instalação do R.")
    try:
        _lsm_r_init(cellsize)  # valida R + pacotes no processo principal
        log("R/landscapemetrics/terra carregados. options(landscapemetrics_cellsize=%s)." % repr(float(cellsize)))
    except Exception as exc:
        raise ImportError(
            "Não foi possível iniciar o R/landscapemetrics via rpy2: %s\n"
            "Verifique se o R está instalado, se R_HOME está definido e se os pacotes "
            "R 'landscapemetrics' e 'terra' foram instalados." % exc)
    set_load("stop", "Ambiente R pronto")

    set_load("start", "Listando rasters em: %s" % pasta_rasters)
    if rasters_lista is not None:  # lista já preparada (ex.: rasters reprojetados para metros)
        rasters = [Path(p) for p in rasters_lista]
    else:
        rasters = listar_rasters_pasta(pasta_rasters, recursivo=recursivo)
    set_load("stop", "Rasters encontrados: %d" % len(rasters))
    if not rasters:
        raise ValueError("Nenhum raster .tif/.tiff foi encontrado na pasta selecionada.")

    if not max_workers or max_workers < 1:
        max_workers = numero_padrao_de_nucleos()

    niveis_ativos = [n for n, ativo in niveis.items() if ativo]
    total_progresso = max(1, len(rasters)) * 100
    set_progress(0, total_progresso, "Aguardando início do cálculo (R)")
    log("Núcleos paralelos usados nas métricas (R): %d" % max_workers)
    log("Níveis (R): %s | vizinhança(directions)=%s | cellsize=%s"
        % (", ".join(niveis_ativos), vizinhanca, repr(float(cellsize))))

    arquivos_individuais: List[Path] = []
    dfs_gerais: List[pd.DataFrame] = []
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    concluidas = {"n": 0}

    def salvar_resultado(raster_path_str, df_final, mensagens):
        raster_path = Path(raster_path_str)
        for m in mensagens:
            log(m)
        if df_final is not None and not df_final.empty:
            ano = extrair_ano_de_recorte(str(raster_path))
            pasta_ano = pasta_saida_path / normalizar_nome(ano, max_len=20)
            pasta_ano.mkdir(parents=True, exist_ok=True)
            nome_csv = normalizar_nome("metricas_paisagem_landscapemetrics_%s" % raster_path.stem, max_len=160) + ".csv"
            caminho_csv = pasta_ano / nome_csv
            df_final.to_csv(caminho_csv, index=False, encoding="utf-8-sig")
            arquivos_individuais.append(caminho_csv)
            dfs_gerais.append(df_final)
            log("CSV individual salvo (R): %s" % caminho_csv)
        else:
            log("Nenhuma métrica (R) foi gerada para: %s" % raster_path.name)
        concluidas["n"] += 1
        set_progress(concluidas["n"] * 100, total_progresso,
                     "Rasters concluídos (R): %d/%d" % (concluidas["n"], len(rasters)))

    tarefas = [(str(rp), niveis_ativos, {k: list(v) for k, v in metricas_func_por_nivel.items()},
                str(vizinhanca), float(cellsize)) for rp in rasters]

    usar_paralelo = max_workers and max_workers > 1
    if usar_paralelo:
        try:
            import multiprocessing as mp
            ctx = mp.get_context("spawn")  # sessão R limpa por processo
            log("Executando em paralelo (processos) — cada processo abre sua própria sessão R.")
            with ProcessPoolExecutor(max_workers=max_workers, mp_context=ctx,
                                     initializer=_lsm_r_init, initargs=(float(cellsize),)) as executor:
                futuros = {executor.submit(_lsm_r_compute, t): t[0] for t in tarefas}
                try:
                    for futuro in as_completed(futuros):
                        checar_parada(stop_event)
                        raster_path_str, df_final, mensagens = futuro.result()
                        salvar_resultado(raster_path_str, df_final, mensagens)
                except ProcessamentoInterrompido:
                    for f in futuros:
                        f.cancel()
                    raise
        except (BrokenProcessPool, OSError) as exc:
            log("Falha no modo paralelo (%s). Recalculando em modo sequencial…" % exc)
            usar_paralelo = False
            concluidas["n"] = 0
            arquivos_individuais.clear()
            dfs_gerais.clear()

    if not usar_paralelo:
        log("Executando em modo sequencial (1 sessão R).")
        for t in tarefas:
            checar_parada(stop_event)
            salvar_resultado(*_lsm_r_compute(t))

    if dfs_gerais:
        df_geral = pd.concat(dfs_gerais, ignore_index=True, sort=False)
        caminho_geral = pasta_saida_path / ("metricas_paisagem_GERAL_landscapemetrics_%s.csv" % timestamp)
        df_geral.to_csv(caminho_geral, index=False, encoding="utf-8-sig")
        log("CSV geral consolidado salvo (R): %s" % caminho_geral)
    else:
        caminho_geral = None

    return {
        "rasters_processados": len(rasters),
        "csvs_individuais": arquivos_individuais,
        "csv_geral": caminho_geral,
    }


# =============================================================================
# Processamento 5: juntar buffer + métricas e exportar shapefiles (.shp)
# =============================================================================

# Limite de tamanho de nome de campo no formato DBF/Shapefile.
LIMITE_NOME_CAMPO_SHP = 10


def normalizar_chave_id(valor) -> str:
    """
    Normaliza um valor para servir de chave textual de pareamento.

    A comparação de IDs é feita como texto (str(valor).strip()), porém um número
    armazenado como float (ex.: 1.0) precisa virar "1" para casar com o ID
    extraído de ``source_file`` (ex.: "1"). Por isso, floats inteiros são
    convertidos para int antes de virar string.
    """
    if valor is None:
        return ""
    if isinstance(valor, float):
        if math.isnan(valor):
            return ""
        if float(valor).is_integer():
            valor = int(valor)
    return str(valor).strip()


def rotulo_ano(valor) -> str:
    """Rótulo seguro de ano para nome de arquivo (1985.0 -> 1985)."""
    if valor is None:
        return "anoNA"
    if isinstance(valor, float):
        if math.isnan(valor):
            return "anoNA"
        if float(valor).is_integer():
            valor = int(valor)
    return normalizar_nome(str(valor), max_len=20)


def gerar_nomes_curtos_campos(
    colunas: Sequence[str],
    limite: int = LIMITE_NOME_CAMPO_SHP,
) -> Dict[str, str]:
    """
    Gera nomes curtos (<= ``limite`` caracteres) e únicos para os campos de um
    shapefile, resolvendo automaticamente colisões causadas pela truncagem.

    Regra:
    - Nome com <= limite caracteres e ainda não usado: mantido como está.
    - Caso contrário, trunca para ``limite`` caracteres; se o truncado for único,
      ele é usado (mais legível).
    - Em caso de colisão, anexa um sufixo numérico incremental (_01, _02, ...),
      reduzindo o prefixo o quanto necessário para caber no limite
      (ex.: perimet_01, perimet_02, ...).

    Retorna um dicionário {nome_original: nome_curto}. A ordem de entrada é
    preservada, de modo que o resultado é determinístico.
    """
    usados: set = set()
    original_para_curto: Dict[str, str] = {}

    for nome in colunas:
        nome = str(nome)
        if len(nome) <= limite and nome not in usados:
            curto = nome
        else:
            truncado = nome[:limite]
            if truncado not in usados:
                curto = truncado
            else:
                curto = None
                contador = 1
                while curto is None:
                    sufixo = f"_{contador:02d}"
                    prefixo = nome[: max(0, limite - len(sufixo))]
                    candidato = f"{prefixo}{sufixo}"
                    if candidato not in usados:
                        curto = candidato
                    contador += 1
        usados.add(curto)
        original_para_curto[nome] = curto

    return original_para_curto


def gerar_shapefiles_metricas_fundido(
    shp_buffers_lista: List[str],
    pasta_saida: str,
    campo_id: str = "id_pt",
    incluir_consolidado: bool = False,
    log_callback=None,
    load_callback=None,
    progress_callback=None,
    stop_event: Optional[threading.Event] = None,
) -> Dict[str, object]:
    """
    Gera shapefiles de métricas da aba Fundido, no MESMO formato da aba 9:
    formato longo (geometria replicada por linha de métrica), nomes de campo
    encurtados para 10 caracteres, um .shp por ano e (opcional) um consolidado,
    além do CSV de equivalência de nomes de campo.

    A junção é feita pela coluna 'id' das métricas (= campo de ID do buffer),
    lendo as métricas já gravadas pela aba Fundido (Parquet por ano, ou CSVs por
    feição) em cada subpasta de ano.
    """
    def log(msg: str):
        if log_callback:
            log_callback(msg)

    def set_load(acao: str, mensagem: str = ""):
        if load_callback:
            load_callback(acao, mensagem)

    def set_progress(valor: float, maximo: float, mensagem: str = ""):
        if progress_callback:
            progress_callback(valor, maximo, mensagem)

    campo_id = (campo_id or "id_pt").strip()
    pasta_saida_path = Path(pasta_saida)
    saida_shp = pasta_saida_path / "shapefiles_metricas"
    saida_shp.mkdir(parents=True, exist_ok=True)

    # ---- 1) Carrega e combina os buffers (CRS do primeiro como referência) ----
    set_load("start", "Carregando shapefile(s) de buffer")
    try:
        gdfs = []
        crs_ref = None
        for shp in shp_buffers_lista:
            g = corrigir_geometrias(ler_shapefile(shp))
            if campo_id not in g.columns:
                raise ValueError(f"O campo '{campo_id}' não existe em {Path(shp).name}.")
            if crs_ref is None:
                crs_ref = g.crs
            elif g.crs != crs_ref:
                g = g.to_crs(crs_ref)
            gdfs.append(g)
        gdf_buffer = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True), crs=crs_ref)
    finally:
        set_load("stop", "Buffers carregados")

    geom_col = gdf_buffer.geometry.name
    buf_attrs = pd.DataFrame(gdf_buffer.drop(columns=[geom_col]))
    buf_attrs["__buf_idx__"] = np.arange(len(gdf_buffer))
    buf_attrs["__join_key__"] = [normalizar_chave_id(v) for v in gdf_buffer[campo_id].tolist()]
    geometrias_buffer = list(gdf_buffer.geometry.values)
    log(f"Feições no buffer combinado: {len(gdf_buffer)} | CRS: {crs_ref}")

    # ---- 2) Descobre as pastas de ano com métricas ----
    pastas_ano = []
    for d in sorted(pasta_saida_path.iterdir()):
        if not d.is_dir():
            continue
        if d.name.startswith("__") or d.name.lower().startswith("shapefiles"):
            continue
        if list(d.glob("*.parquet")) or list(d.glob("*.csv")):
            pastas_ano.append(d)
    if not pastas_ano:
        raise ValueError("Nenhuma pasta de ano com Parquet/CSV encontrada. "
                         "Rode o Fundido com 'Parquet por ano' ou 'CSV por feição' ligado.")

    tem_pyarrow = _modulo_disponivel("pyarrow")

    def ler_metricas_do_ano(pasta_ano: Path) -> "pd.DataFrame":
        parquets = sorted(pasta_ano.glob("*.parquet"))
        if parquets and tem_pyarrow:
            return pd.concat([pd.read_parquet(p) for p in parquets], ignore_index=True)
        csvs = sorted(pasta_ano.glob("*.csv"))
        if csvs:
            return pd.concat([pd.read_csv(c, encoding="utf-8-sig", low_memory=False) for c in csvs], ignore_index=True)
        if parquets and not tem_pyarrow:
            raise ValueError("Há Parquets por ano, mas 'pyarrow' não está instalado para lê-los. "
                             "Instale 'pip install pyarrow' ou use os CSVs por feição.")
        return pd.DataFrame()

    mapa_original_para_curto = None
    year_curto = "year"
    arquivos_saida: List[Path] = []
    acumulado = [] if incluir_consolidado else None
    total = len(pastas_ano)
    feicoes_total = 0

    for idx, pasta_ano in enumerate(pastas_ano, start=1):
        checar_parada(stop_event)
        set_progress(idx - 1, total, f"Gerando shapefile {idx}/{total} (ano {pasta_ano.name})")
        df_ano = ler_metricas_do_ano(pasta_ano)
        if df_ano.empty:
            log(f"Ano {pasta_ano.name}: sem dados — pulado.")
            continue
        df_ano.columns = [str(c).strip() for c in df_ano.columns]
        if "id" not in df_ano.columns:
            log(f"Ano {pasta_ano.name}: métricas sem coluna 'id' — pulado.")
            continue
        df_ano["__join_key__"] = [normalizar_chave_id(v) for v in df_ano["id"].tolist()]

        # Resolve colisões de nome entre atributos do buffer e métricas.
        colisoes = (set(buf_attrs.columns) & set(df_ano.columns)) - {"__join_key__"}
        if colisoes:
            df_ano = df_ano.rename(columns={c: f"{c}_m" for c in colisoes})

        merged = buf_attrs.merge(df_ano, on="__join_key__", how="inner")
        if merged.empty:
            log(f"Ano {pasta_ano.name}: nenhuma correspondência id↔buffer — pulado.")
            continue
        geom_saida = [geometrias_buffer[i] for i in merged["__buf_idx__"].to_numpy()]
        merged = merged.drop(columns=["__buf_idx__", "__join_key__"])
        gdf_out = gpd.GeoDataFrame(merged, geometry=geom_saida, crs=crs_ref)

        # Nomes curtos: calculados uma vez (mesmas colunas em todos os anos).
        if mapa_original_para_curto is None:
            colunas_atributo = [c for c in gdf_out.columns if c != gdf_out.geometry.name]
            mapa_original_para_curto = gerar_nomes_curtos_campos(colunas_atributo)
            year_curto = mapa_original_para_curto.get("year", "year")
            mapa_df = pd.DataFrame(
                [(curto, orig) for orig, curto in mapa_original_para_curto.items()],
                columns=["nome_curto", "nome_original"],
            )
            mapa_df.to_csv(saida_shp / "buffer_metricas_mapa_campos.csv", index=False, encoding="utf-8-sig")

        gdf_curto = filtrar_colunas_shapefile(gdf_out.rename(columns=mapa_original_para_curto))
        nome_base = f"buffer_metricas_{rotulo_ano(pasta_ano.name)}"
        shp_path = saida_shp / f"{nome_base}.shp"
        gdf_curto.to_file(shp_path, driver="ESRI Shapefile", encoding="utf-8", index=False, engine="pyogrio")
        arquivos_saida.append(shp_path)
        feicoes_total += len(gdf_curto)
        log(f"Ano {pasta_ano.name}: {len(gdf_curto)} feições -> {shp_path.name}")
        if acumulado is not None:
            acumulado.append(gdf_curto)

    consolidado = None
    if incluir_consolidado and acumulado:
        checar_parada(stop_event)
        set_progress(total, total, "Gravando shapefile consolidado (todos os anos)")
        try:
            gdf_todos = gpd.GeoDataFrame(pd.concat(acumulado, ignore_index=True), crs=crs_ref)
            consolidado = saida_shp / "buffer_metricas_TODOS_ANOS.shp"
            gdf_todos.to_file(consolidado, driver="ESRI Shapefile", encoding="utf-8", index=False, engine="pyogrio")
            log(f"Consolidado: {len(gdf_todos)} feições -> {consolidado.name}")
        except Exception as exc:
            log(f"AVISO: não foi possível gerar o consolidado (provável limite de tamanho do shapefile): {exc}")
            consolidado = None

    set_progress(total, total, "Shapefiles de métricas concluídos.")
    return {
        "shapefiles": arquivos_saida,
        "consolidado": consolidado,
        "feicoes_total": feicoes_total,
        "pasta": saida_shp,
    }


def gerar_shapefiles_buffer_metricas(
    shp_buffer: str,
    csv_metricas: str,
    pasta_saida: str,
    campo_id: str = "id_pt",
    tipo_entrada: str = "csv",
    log_callback=None,
    load_callback=None,
    progress_callback=None,
    stop_event: Optional[threading.Event] = None,
) -> Dict[str, object]:
    """
    Enriquece cada feição de um shapefile de buffer com as métricas de paisagem
    (PyLandStats) exportadas em CSV, gerando shapefiles (.shp) em formato longo:
    uma feição replicada por linha de métrica (várias classes x vários anos).

    Saídas (somente .shp):
    - Um shapefile por ano: buffer_metricas_<ano>.shp
    - Um consolidado com todos os anos: buffer_metricas_TODOS_ANOS.shp
    - Um CSV de equivalência nome_curto -> nome_original (mapa de campos)
    """
    def log(msg: str):
        if log_callback:
            log_callback(msg)

    def set_load(acao: str, mensagem: str = ""):
        if load_callback:
            load_callback(acao, mensagem)

    def set_progress(valor: float, maximo: float, mensagem: str = ""):
        if progress_callback:
            progress_callback(valor, maximo, mensagem)

    campo_id = (campo_id or "id_pt").strip()
    pasta_saida = Path(pasta_saida)
    pasta_saida.mkdir(parents=True, exist_ok=True)

    set_progress(0, 100, "Aguardando início do pareamento")

    # ------------------------------------------------------------------
    # 1) Leitura do shapefile de buffer
    # ------------------------------------------------------------------
    checar_parada(stop_event)
    set_load("start", f"Carregando shapefile de buffer: {Path(shp_buffer).name}")
    try:
        gdf_buffer = ler_shapefile(shp_buffer)
    finally:
        set_load("stop", "Carregamento do buffer finalizado")
    checar_parada(stop_event)

    crs_buffer = gdf_buffer.crs
    geom_col = gdf_buffer.geometry.name
    log(f"Feições no buffer: {len(gdf_buffer):,}".replace(",", "."))
    log(f"CRS do buffer (mantido na saída): {crs_buffer}")

    if campo_id not in gdf_buffer.columns:
        raise ValueError(
            f"O campo de ID '{campo_id}' não existe na tabela do buffer. "
            f"Campos disponíveis: {', '.join([c for c in gdf_buffer.columns if c != geom_col])}"
        )

    # ------------------------------------------------------------------
    # 2) Leitura das métricas (CSV ou Parquet) e definição da chave de junção
    # ------------------------------------------------------------------
    tipo_entrada = (tipo_entrada or "csv").lower()
    checar_parada(stop_event)

    if tipo_entrada == "parquet":
        if not _modulo_disponivel("pyarrow"):
            raise ValueError("Leitura de Parquet requer 'pyarrow'. Instale com: pip install pyarrow")
        fonte = Path(csv_metricas)
        set_load("start", f"Carregando Parquet de métricas: {fonte.name}")
        try:
            if fonte.is_dir():
                arqs = sorted(fonte.rglob("*.parquet"))
                if not arqs:
                    raise ValueError("Nenhum arquivo .parquet encontrado na pasta selecionada.")
                df_csv = pd.concat([pd.read_parquet(a) for a in arqs], ignore_index=True)
                log(f"Parquets lidos da pasta: {len(arqs)}")
            else:
                df_csv = pd.read_parquet(fonte)
        finally:
            set_load("stop", "Carregamento do Parquet finalizado")
        df_csv.columns = [str(c).strip() for c in df_csv.columns]
        log(f"Linhas no Parquet de métricas: {len(df_csv):,}".replace(",", "."))
        log(f"Colunas no Parquet de métricas: {len(df_csv.columns)}")

        # Pareamento direto pela coluna de ID (id_pt de ambos; aceita 'id' como alternativa).
        col_id_metricas = None
        for cand in (campo_id, "id_pt", "id"):
            if cand in df_csv.columns:
                col_id_metricas = cand
                break
        if col_id_metricas is None:
            raise ValueError(f"O Parquet não tem coluna de ID ('{campo_id}', 'id_pt' ou 'id') para parear.")
        log(f"Coluna de ID usada no Parquet: '{col_id_metricas}' (pareada com '{campo_id}' do buffer)")
        if "year" not in df_csv.columns:
            df_csv["year"] = "anoNA"
        df_csv["__join_key__"] = [normalizar_chave_id(v) for v in df_csv[col_id_metricas].tolist()]
        df_csv = df_csv[df_csv["__join_key__"] != ""].copy()
    else:
        set_load("start", f"Carregando CSV de métricas: {Path(csv_metricas).name}")
        try:
            # utf-8-sig remove eventual BOM no cabeçalho.
            df_csv = pd.read_csv(
                csv_metricas,
                sep=",",
                decimal=".",
                encoding="utf-8-sig",
                low_memory=False,
            )
        finally:
            set_load("stop", "Carregamento do CSV finalizado")
        checar_parada(stop_event)

        df_csv.columns = [str(c).strip() for c in df_csv.columns]
        log(f"Linhas no CSV de métricas: {len(df_csv):,}".replace(",", "."))
        log(f"Colunas no CSV de métricas: {len(df_csv.columns)}")

        if "year" not in df_csv.columns:
            raise ValueError("O CSV de métricas não possui a coluna obrigatória 'year'.")

        # Garante que colunas numéricas sejam interpretadas como número (decimal '.').
        colunas_texto_protegidas = {"source_path", "source_file", "source_stem", "level", "neighborhood_rule"}
        convertidas = 0
        for coluna in df_csv.columns:
            if coluna in colunas_texto_protegidas:
                continue
            if df_csv[coluna].dtype == object:
                convertido = pd.to_numeric(df_csv[coluna], errors="coerce")
                if convertido.notna().sum() >= df_csv[coluna].notna().sum():
                    df_csv[coluna] = convertido
                    convertidas += 1
        if convertidas:
            log(f"Colunas convertidas para numérico (decimal '.'): {convertidas}")

        # Chave de junção: 'id'/'id_pt' direto se existir; senão extrai de source_file (id_<N>).
        col_id_direta = next((c for c in (campo_id, "id_pt", "id") if c in df_csv.columns), None)
        if col_id_direta is not None:
            log(f"Pareamento pela coluna de ID do CSV: '{col_id_direta}'")
            df_csv["__join_key__"] = [normalizar_chave_id(v) for v in df_csv[col_id_direta].tolist()]
        elif "source_file" in df_csv.columns:
            log("Pareamento pelo ID extraído de 'source_file' (padrão id_<N>)")
            df_csv["__join_key__"] = (
                df_csv["source_file"].astype(str).str.extract(r"id_(\d+)", expand=False).str.strip()
            )
        else:
            raise ValueError("O CSV não tem coluna de ID ('id_pt'/'id') nem 'source_file' para parear.")
        sem_id = int(df_csv["__join_key__"].isna().sum())
        if sem_id:
            log(f"Aviso: {sem_id} linha(s) sem ID para pareamento (ignoradas).")
        df_csv = df_csv[df_csv["__join_key__"].notna()].copy()
        df_csv = df_csv[df_csv["__join_key__"].astype(str) != ""].copy()

    # ------------------------------------------------------------------
    # 4) Chave textual do buffer (str().strip(), robusto a número/float)
    # ------------------------------------------------------------------
    buf_attrs = pd.DataFrame(gdf_buffer.drop(columns=[geom_col]))
    buf_attrs["__buf_idx__"] = np.arange(len(gdf_buffer))
    buf_attrs["__join_key__"] = [normalizar_chave_id(v) for v in gdf_buffer[campo_id].tolist()]
    geometrias_buffer = list(gdf_buffer.geometry.values)

    # ------------------------------------------------------------------
    # 5) Diagnóstico de pareamento (registra, mas não interrompe)
    # ------------------------------------------------------------------
    chaves_buffer = set(buf_attrs["__join_key__"]) - {""}
    chaves_csv = set(df_csv["__join_key__"]) - {""}
    csv_sem_buffer = sorted(chaves_csv - chaves_buffer, key=lambda s: (len(s), s))
    buffer_sem_csv = sorted(chaves_buffer - chaves_csv, key=lambda s: (len(s), s))

    log("-" * 70)
    log(f"IDs distintos no buffer: {len(chaves_buffer)} | IDs distintos no CSV: {len(chaves_csv)}")
    if csv_sem_buffer:
        amostra = ", ".join(csv_sem_buffer[:20])
        log(f"IDs no CSV SEM feição correspondente no buffer: {len(csv_sem_buffer)} (ex.: {amostra}{' ...' if len(csv_sem_buffer) > 20 else ''})")
    if buffer_sem_csv:
        amostra = ", ".join(buffer_sem_csv[:20])
        log(f"IDs no buffer SEM linha correspondente no CSV: {len(buffer_sem_csv)} (ex.: {amostra}{' ...' if len(buffer_sem_csv) > 20 else ''})")

    # ------------------------------------------------------------------
    # 6) Join 1-para-N (replicação da geometria por linha de métrica)
    # ------------------------------------------------------------------
    checar_parada(stop_event)
    set_progress(20, 100, "Pareando buffer com métricas (join 1-para-N)")

    # Resolve colisão de nomes entre atributos do buffer e colunas do CSV,
    # preservando ambos (a coluna do CSV recebe sufixo '_m').
    colisoes = (set(buf_attrs.columns) & set(df_csv.columns)) - {"__join_key__"}
    if colisoes:
        renomeio_csv = {c: f"{c}_m" for c in colisoes}
        df_csv = df_csv.rename(columns=renomeio_csv)
        log(f"Colunas do CSV renomeadas por colisão com o buffer: {', '.join(f'{k}->{v}' for k, v in renomeio_csv.items())}")

    merged = buf_attrs.merge(df_csv, on="__join_key__", how="inner")
    if merged.empty:
        raise ValueError(
            "Nenhuma correspondência entre o buffer e o CSV. "
            "Verifique o 'Campo ID' selecionado e o padrão de source_file."
        )

    # Reconstrói a geometria replicada a partir do índice do buffer.
    geom_saida = [geometrias_buffer[i] for i in merged["__buf_idx__"].to_numpy()]
    merged = merged.drop(columns=["__buf_idx__", "__join_key__"])

    gdf_out = gpd.GeoDataFrame(merged, geometry=geom_saida, crs=crs_buffer)
    log("-" * 70)
    log(f"Feições geradas (formato longo): {len(gdf_out):,}".replace(",", "."))

    # ------------------------------------------------------------------
    # 7) Nomes curtos (<=10) únicos para os campos do shapefile
    # ------------------------------------------------------------------
    colunas_atributo = [c for c in gdf_out.columns if c != gdf_out.geometry.name]
    mapa_original_para_curto = gerar_nomes_curtos_campos(colunas_atributo)
    renomeados = {orig: curto for orig, curto in mapa_original_para_curto.items() if orig != curto}
    if renomeados:
        log(f"Campos renomeados para caber em 10 caracteres (.shp): {len(renomeados)}")
        for orig, curto in renomeados.items():
            log(f"   {orig} -> {curto}")
    else:
        log("Nenhum campo precisou ser encurtado para o formato .shp.")

    # Salva o CSV de equivalência (nome_curto -> nome_original completo).
    mapa_df = pd.DataFrame(
        [(curto, orig) for orig, curto in mapa_original_para_curto.items()],
        columns=["nome_curto", "nome_original"],
    )
    caminho_mapa = pasta_saida / "buffer_metricas_mapa_campos.csv"
    mapa_df.to_csv(caminho_mapa, index=False, encoding="utf-8-sig")
    log(f"Mapa de equivalência de campos salvo: {caminho_mapa}")

    gdf_out_curto = gdf_out.rename(columns=mapa_original_para_curto)
    gdf_out_curto = filtrar_colunas_shapefile(gdf_out_curto)

    # Nome curto correspondente à coluna 'year' (para separar por ano).
    year_curto = mapa_original_para_curto.get("year", "year")

    # ------------------------------------------------------------------
    # 8) Exportação: um .shp por ano + um consolidado
    # ------------------------------------------------------------------
    anos = list(pd.unique(gdf_out_curto[year_curto]))
    # Ordena anos de forma estável (numérica quando possível).
    def _ordenar_ano(v):
        try:
            return (0, float(v))
        except (TypeError, ValueError):
            return (1, str(v))
    anos = sorted(anos, key=_ordenar_ano)

    total_etapas = max(1, len(anos) + 1)
    arquivos_saida: List[Path] = []

    for idx_ano, ano in enumerate(anos, start=1):
        checar_parada(stop_event)
        gdf_ano = gdf_out_curto[gdf_out_curto[year_curto] == ano].copy()
        nome_base = f"buffer_metricas_{rotulo_ano(ano)}"
        shp_path = pasta_saida / f"{nome_base}.shp"
        gdf_ano.to_file(shp_path, driver="ESRI Shapefile", encoding="utf-8", index=False, engine="pyogrio")
        arquivos_saida.append(shp_path)
        log(f"[{rotulo_ano(ano)}] feições: {len(gdf_ano):,}".replace(",", ".") + f" -> {shp_path.name}")
        set_progress(20 + (idx_ano / total_etapas) * 75, 100, f"Salvando shapefile do ano {rotulo_ano(ano)}")

    # Consolidado com todos os anos.
    checar_parada(stop_event)
    shp_consolidado = pasta_saida / "buffer_metricas_TODOS_ANOS.shp"
    gdf_out_curto.to_file(shp_consolidado, driver="ESRI Shapefile", encoding="utf-8", index=False, engine="pyogrio")
    arquivos_saida.append(shp_consolidado)
    set_progress(100, 100, "Concluído")
    log(f"[TODOS_ANOS] feições: {len(gdf_out_curto):,}".replace(",", ".") + f" -> {shp_consolidado.name}")

    return {
        "shapefiles": arquivos_saida,
        "consolidado": shp_consolidado,
        "mapa_campos": caminho_mapa,
        "anos": [rotulo_ano(a) for a in anos],
        "feicoes_total": len(gdf_out_curto),
        "csv_sem_buffer": csv_sem_buffer,
        "buffer_sem_csv": buffer_sem_csv,
        "campos_renomeados": renomeados,
    }


# =============================================================================
# Processamento 6: reprojeção de rasters em lote
# =============================================================================

RESAMPLING_REPROJECAO = {
    "nearest (categórico)": Resampling.nearest,
    "bilinear": Resampling.bilinear,
    "cubic": Resampling.cubic,
    "mode": Resampling.mode,
}


def reprojetar_raster_arquivo(
    caminho_entrada: str,
    caminho_saida: Path,
    epsg_destino: int,
    resampling: Resampling = Resampling.nearest,
    num_threads: int = 1,
    stop_event: Optional[threading.Event] = None,
) -> Path:
    """Reprojeta um raster inteiro para o EPSG informado, gravando em GeoTIFF."""
    with rasterio.open(caminho_entrada) as src:
        if src.crs is None:
            raise ValueError(f"O raster não possui CRS definido: {Path(caminho_entrada).name}")

        crs_destino = CRS.from_epsg(int(epsg_destino))
        dst_transform, dst_width, dst_height = calculate_default_transform(
            src.crs, crs_destino, src.width, src.height, *src.bounds
        )

        perfil = src.profile.copy()
        perfil.update(
            driver="GTiff",
            crs=crs_destino,
            transform=dst_transform,
            width=dst_width,
            height=dst_height,
            compress="lzw",
            tiled=True,
            BIGTIFF="IF_SAFER",
        )

        caminho_saida.parent.mkdir(parents=True, exist_ok=True)
        num_threads = max(1, int(num_threads))
        remover_auxiliares_raster(caminho_saida)
        with rasterio.open(caminho_saida, "w", **perfil) as dst:
            for banda in range(1, src.count + 1):
                checar_parada(stop_event)
                reproject(
                    source=rasterio.band(src, banda),
                    destination=rasterio.band(dst, banda),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    src_nodata=src.nodata,
                    dst_transform=dst_transform,
                    dst_crs=crs_destino,
                    dst_nodata=src.nodata,
                    resampling=resampling,
                    num_threads=num_threads,
                )
    return caminho_saida


def reprojetar_rasters_lote(
    rasters_lista: List[str],
    pasta_saida: str,
    epsg_destino: int,
    sufixo_epsg: bool = True,
    resampling: Resampling = Resampling.nearest,
    max_workers: Optional[int] = None,
    log_callback=None,
    load_callback=None,
    progress_callback=None,
    stop_event: Optional[threading.Event] = None,
    pasta_base: Optional[str] = None,
    normalizar_nomes: bool = True,
) -> Dict[str, object]:
    """Reprojeta vários rasters em paralelo. Erros por arquivo não abortam o lote.

    Com `pasta_base`, a estrutura de subpastas de cada raster em relação a essa
    pasta é replicada dentro de `pasta_saida`. O resultado traz, além das
    saídas e dos erros, o dicionário "mapa" {raster de entrada: raster de saída}.
    """
    def log(msg: str):
        if log_callback:
            log_callback(msg)

    def set_load(acao: str, mensagem: str = ""):
        if load_callback:
            load_callback(acao, mensagem)

    def set_progress(valor: float, maximo: float, mensagem: str = ""):
        if progress_callback:
            progress_callback(valor, maximo, mensagem)

    if not rasters_lista:
        raise ValueError("Selecione pelo menos um raster para reprojetar.")
    if not pasta_saida:
        raise ValueError("Selecione a pasta de saída dos rasters reprojetados.")
    try:
        epsg_destino = int(epsg_destino)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"EPSG de destino inválido: {epsg_destino}") from exc

    if not max_workers or max_workers < 1:
        max_workers = numero_padrao_de_nucleos()

    pasta_saida_path = Path(pasta_saida)
    pasta_saida_path.mkdir(parents=True, exist_ok=True)

    cpus = os.cpu_count() or 4
    reproject_threads = max(1, cpus // max(1, max_workers))

    total_progresso = max(1, len(rasters_lista)) * 100
    set_progress(0, total_progresso, "Aguardando início da reprojeção")
    log(f"EPSG de destino: {epsg_destino}")
    log(f"Núcleos paralelos: {max_workers}")
    if pasta_base:
        log(f"Estrutura de subpastas replicada a partir de: {pasta_base}")

    contador_lock = threading.Lock()
    concluidas = {"n": 0}
    saidas: List[Path] = []
    erros: List[str] = []
    mapa: Dict[str, str] = {}

    def tarefa(caminho_entrada: str):
        checar_parada(stop_event)
        nome = Path(caminho_entrada).name
        try:
            stem = Path(caminho_entrada).stem
            if normalizar_nomes:
                stem = normalizar_nome(stem, max_len=120)
            sufixo = f"_{epsg_destino}" if sufixo_epsg else ""
            caminho_saida = pasta_saida_path / _subpasta_relativa(caminho_entrada, pasta_base) / f"{stem}{sufixo}.tif"
            if os.path.normcase(os.path.abspath(caminho_saida)) == os.path.normcase(os.path.abspath(caminho_entrada)):
                raise ValueError("o arquivo de saída coincide com o de entrada; escolha outra pasta de saída")
            reprojetar_raster_arquivo(
                caminho_entrada=caminho_entrada,
                caminho_saida=caminho_saida,
                epsg_destino=epsg_destino,
                resampling=resampling,
                num_threads=reproject_threads,
                stop_event=stop_event,
            )
            return ("ok", nome, caminho_saida, None)
        except ProcessamentoInterrompido:
            raise
        except Exception as exc:
            return ("erro", nome, None, str(exc))

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futuros = {executor.submit(tarefa, rp): rp for rp in rasters_lista}
        try:
            for futuro in as_completed(futuros):
                if stop_event is not None and stop_event.is_set():
                    for f in futuros:
                        f.cancel()
                    raise ProcessamentoInterrompido("Processamento interrompido pelo usuário.")
                status, nome, caminho, erro = futuro.result()
                if status == "ok":
                    saidas.append(caminho)
                    mapa[str(futuros[futuro])] = str(caminho)
                    log(f"Reprojetado: {nome} -> {os.path.relpath(caminho, pasta_saida_path)}")
                else:
                    erros.append(f"{nome}: {erro}")
                    log(f"ERRO ao reprojetar {nome}: {erro}")
                with contador_lock:
                    concluidas["n"] += 1
                    n = concluidas["n"]
                set_progress(n * 100, total_progresso, f"Rasters concluídos: {n}/{len(rasters_lista)}")
        except ProcessamentoInterrompido:
            for f in futuros:
                f.cancel()
            raise

    return {"saidas": saidas, "erros": erros, "total": len(rasters_lista), "mapa": mapa}


def _subpasta_relativa(caminho_arquivo: str, pasta_base: Optional[str]) -> str:
    """Subpasta do arquivo em relação a pasta_base ("" na raiz ou fora dela)."""
    if not pasta_base:
        return ""
    try:
        rel = os.path.relpath(os.path.dirname(os.path.abspath(caminho_arquivo)), os.path.abspath(pasta_base))
    except ValueError:  # unidades de disco diferentes
        return ""
    return "" if rel == "." or rel.startswith("..") else rel


# =============================================================================
# Processamento 9: reclassificação de rasters por regras (Val_Min, Val_Max, Novo_Val)
# Portado da aba "Reclassificar" do GeoRasterTools (equivalente ao terra::classify).
# =============================================================================

# Regras padrão (nomenclatura MapBiomas): cada tupla é (Val_Min, Val_Max, Novo_Val).
# Arquivos que ArcGIS/QGIS/GDAL criam ao lado de um raster: tabela de atributos (.vat.dbf),
# estatísticas e histograma (.aux.xml) e pirâmides (.ovr). Se o raster é sobrescrito e eles
# ficam, o ArcGIS continua mostrando as classes ANTIGAS (ex.: valor 1 que não existe mais).
SUFIXOS_AUXILIARES_RASTER = (".aux.xml", ".aux", ".ovr", ".vat.dbf", ".vat.cpg")


def remover_auxiliares_raster(caminho) -> List[str]:
    """Apaga os arquivos auxiliares de `caminho`; devolve os que não puderam ser
    apagados (ex.: raster aberto no ArcGIS)."""
    presos: List[str] = []
    for sufixo in SUFIXOS_AUXILIARES_RASTER:
        aux = str(caminho) + sufixo
        if os.path.exists(aux):
            try:
                os.remove(aux)
            except OSError:
                presos.append(os.path.basename(aux))
    return presos


REGRAS_RECLASS_PADRAO = [
    (1, 1, 1), (3, 3, 1), (4, 4, 1), (5, 5, 1),
    (6, 6, 1), (49, 49, 1), (10, 10, 1), (11, 11, 1),
    (12, 12, 1), (32, 32, 1), (29, 29, 1), (50, 50, 1),
    (14, 14, 14), (15, 15, 14), (18, 18, 14), (19, 19, 14),
    (39, 39, 14), (20, 20, 14), (40, 40, 14), (62, 62, 14),
    (41, 41, 14), (36, 36, 14), (46, 46, 14), (47, 47, 14),
    (35, 35, 14), (48, 48, 14), (9, 9, 14), (21, 21, 14),
    (22, 22, 22), (23, 23, 22), (75, 75, 22), (25, 25, 22),
    (24, 24, 24),
    (30, 30, 22),
    (26, 26, 26), (33, 33, 26), (31, 31, 26),
]
COLUNAS_REGRA_RECLASS = ("Val_Min", "Val_Max", "Novo_Val")
NOME_REGRA_PADRAO = "Padrao"

# Tipos de dado de saída: código no estilo terra/GDAL -> (dtype numpy, bits).
RECLASS_DTYPES = {
    "INT1U": ("uint8", 8), "INT2U": ("uint16", 16), "INT2S": ("int16", 16), "INT4S": ("int32", 32),
    "FLT4S": ("float32", 32), "FLT8S": ("float64", 64), "BYTE": ("uint8", 8),
}
# Nodata gravado quando os pixels não mapeados ficam "NA" (sem valor).
RECLASS_NODATA_NA = {"uint8": 255, "int16": -9999, "uint16": 0, "int32": -9999}
RECLASS_NAO_MAPEADOS = {"NA (sem valor)": "NA", "Zero (0)": "0"}
RECLASS_INTERVALOS = {"Ambos os lados (NA)": None, "Direita (TRUE)": True, "Esquerda (FALSE)": False}


def rotulo_dtype_reclass(codigo: str) -> str:
    """'INT2U' -> 'INT2U - 16 bits - uint16' (o número do código conta BYTES, não bits)."""
    np_dtype, bits = RECLASS_DTYPES[codigo]
    return f"{codigo} - {bits} bits - {np_dtype}"


def regras_reclass_padrao() -> List[tuple]:
    """Cópia das regras padrão."""
    return [tuple(map(float, r)) for r in REGRAS_RECLASS_PADRAO]


def ler_regras_reclass_txt(caminho: str) -> List[tuple]:
    """Lê um .txt com linhas 'Val_Min,Val_Max,Novo_Val' (separador: vírgula, ponto e
    vírgula, tab ou espaço), ignorando comentários (#) e linhas em branco."""
    regras: List[tuple] = []
    with open(caminho, "r", encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if re.match(r"^\s*#|^\s*$", linha):
                continue
            nums: List[float] = []
            for parte in re.split(r"[,;\t ]+", linha.strip()):
                try:
                    nums.append(float(parte))
                except ValueError:
                    continue
            if len(nums) >= 3:
                regras.append((nums[0], nums[1], nums[2]))
    return regras


def salvar_regras_reclass_txt(regras, caminho: str, nome_regra: str = "Regra") -> None:
    """Grava as regras no mesmo formato de texto usado pelo GeoRasterTools."""
    linhas = [
        "# GeoLand - Regras de Reclassificacao",
        f"# Regra: {nome_regra}",
        "# Formato: Val_Min,Val_Max,Novo_Val",
        "# ----------------------------------------",
    ] + [f"{r[0]:g},{r[1]:g},{r[2]:g}" for r in regras]
    with open(caminho, "w", encoding="utf-8") as fh:
        fh.write("\n".join(linhas) + "\n")


def reclassificar_array(dados: np.ndarray, regras, outros: Optional[float], direita: Optional[bool] = None,
                        incluir_minimo: bool = True, validos: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Aplica as regras a um array (resultado em float64).

    `direita` segue o terra::classify: True -> (min, max]; False -> [min, max);
    None -> [min, max]. Com direita=True e incluir_minimo=True, o menor Val_Min
    de todas as regras também é incluído. As regras valem na ordem dada: a
    primeira que casar define o pixel. `outros` é o valor dos pixels que não
    casam com nenhuma regra (None = NaN) e dos pixels fora de `validos`.
    """
    vazio = np.nan if outros is None else float(outros)
    saida = np.full(dados.shape, vazio, dtype=np.float64)
    atribuido = np.zeros(dados.shape, dtype=bool)
    menor = min(r[0] for r in regras) if regras else None
    for vmin, vmax, novo in regras:
        if direita is True:
            m = (dados > vmin) & (dados <= vmax)
            if incluir_minimo and vmin == menor:
                m |= dados == vmin
        elif direita is False:
            m = (dados >= vmin) & (dados < vmax)
        else:
            m = (dados >= vmin) & (dados <= vmax)
        m &= ~atribuido
        saida[m] = novo
        atribuido |= m
    if validos is not None:
        saida[~validos] = vazio
    return saida


def reclassificar_raster_arquivo(
    caminho_entrada: str,
    regras,
    caminho_saida: str,
    dtype: str = "INT1U",
    outros: str = "NA",
    direita: Optional[bool] = None,
    incluir_minimo: bool = True,
    sobrescrever: bool = True,
    stop_event: Optional[threading.Event] = None,
) -> List[str]:
    """
    Reclassifica todas as bandas de um raster e grava em GeoTIFF (LZW).
    Apaga a tabela de atributos/estatísticas antigas da saída (.vat.dbf, .aux.xml...)
    e devolve os auxiliares que não puderam ser apagados.

    Processa em faixas de linhas, para a memória não crescer com o tamanho do
    raster. Em rasters inteiros de até 16 bits usa uma tabela com o resultado
    de cada valor possível, muito mais rápida que testar regra a regra em cada
    pixel (o resultado é o mesmo).
    """
    from rasterio.windows import Window

    if not sobrescrever and os.path.exists(caminho_saida):
        raise FileExistsError(f"arquivo já existe: {caminho_saida}")
    if os.path.normcase(os.path.abspath(caminho_saida)) == os.path.normcase(os.path.abspath(caminho_entrada)):
        raise ValueError("o arquivo de saída coincide com o de entrada; mude a pasta de saída ou o prefixo")

    outros_val = None if outros == "NA" else 0.0
    vazio = np.nan if outros_val is None else outros_val
    np_dtype = RECLASS_DTYPES[dtype][0]
    nodata_saida = RECLASS_NODATA_NA.get(np_dtype, np.nan) if outros_val is None else 0

    with rasterio.open(caminho_entrada) as src:
        perfil = src.meta.copy()
        perfil.update(driver="GTiff", dtype=np_dtype, nodata=nodata_saida, compress="LZW", tiled=True,
                      BIGTIFF="IF_SAFER")
        tabela, minimo = None, 0
        tipo_entrada = np.dtype(src.dtypes[0])
        if len(set(src.dtypes)) == 1 and np.issubdtype(tipo_entrada, np.integer) and tipo_entrada.itemsize <= 2:
            info = np.iinfo(tipo_entrada)
            minimo = int(info.min)
            tabela = reclassificar_array(np.arange(info.min, info.max + 1, dtype=np.float64),
                                         regras, outros_val, direita, incluir_minimo)
        linhas_por_faixa = max(1, 4_000_000 // max(1, src.width * src.count))
        Path(caminho_saida).parent.mkdir(parents=True, exist_ok=True)
        remover_auxiliares_raster(caminho_saida)
        try:
            with rasterio.open(caminho_saida, "w", **perfil) as dst:
                for linha in range(0, src.height, linhas_por_faixa):
                    checar_parada(stop_event)
                    janela = Window(0, linha, src.width, min(linhas_por_faixa, src.height - linha))
                    bruto = src.read(window=janela)
                    if tabela is not None:
                        resultado = tabela[bruto.astype(np.int64) - minimo]
                        if src.nodata is not None:
                            resultado[bruto == src.nodata] = vazio
                    else:
                        dados = bruto.astype(np.float64)
                        validos = np.isfinite(dados)
                        if src.nodata is not None:
                            validos &= dados != src.nodata
                        resultado = reclassificar_array(dados, regras, outros_val, direita, incluir_minimo, validos)
                    if outros_val is None:
                        resultado = np.where(np.isnan(resultado), nodata_saida, resultado)
                    dst.write(resultado.astype(np_dtype), window=janela)
        except BaseException:
            try:  # não deixa um raster pela metade na pasta de saída
                os.remove(caminho_saida)
            except OSError:
                pass
            raise
    return remover_auxiliares_raster(caminho_saida)


def reclassificar_rasters_lote(
    arquivos: Sequence[tuple],
    regras_por_nome: Dict[str, list],
    pasta_saida: str,
    dtype: str = "INT1U",
    outros: str = "NA",
    direita: Optional[bool] = None,
    incluir_minimo: bool = True,
    sobrescrever: bool = True,
    prefixo: str = "Reclass_",
    pasta_base: Optional[str] = None,
    max_workers: Optional[int] = None,
    log_callback=None,
    progress_callback=None,
    stop_event: Optional[threading.Event] = None,
) -> Dict[str, object]:
    """
    Reclassifica vários rasters em paralelo. `arquivos` é uma lista de
    (caminho do raster, nome da regra) e `regras_por_nome` traz as regras de
    cada nome. A saída é <pasta_saida>/<subpasta>/<prefixo><nome do arquivo>;
    com `pasta_base`, as subpastas da entrada são replicadas. Erros por
    arquivo não abortam o lote.
    """
    def log(msg: str):
        if log_callback:
            log_callback(msg)

    if not arquivos:
        raise ValueError("Nenhum raster para reclassificar.")
    if not pasta_saida:
        raise ValueError("Selecione a pasta de saída.")
    if dtype not in RECLASS_DTYPES:
        raise ValueError(f"Tipo de dado de saída inválido: {dtype}")
    if not max_workers or max_workers < 1:
        max_workers = numero_padrao_de_nucleos()
    pasta_saida_path = Path(pasta_saida)
    pasta_saida_path.mkdir(parents=True, exist_ok=True)

    # Avisos sobre valores que o tipo de dado de saída não consegue guardar.
    tipo_saida = np.dtype(RECLASS_DTYPES[dtype][0])
    nodata_saida = RECLASS_NODATA_NA.get(tipo_saida.name, np.nan) if outros == "NA" else 0
    for nome_regra in sorted({r for _, r in arquivos}):
        novos = sorted({r[2] for r in regras_por_nome.get(nome_regra) or []})
        if np.issubdtype(tipo_saida, np.integer):
            info = np.iinfo(tipo_saida)
            fora = [v for v in novos if v < info.min or v > info.max or v != int(v)]
            if fora:
                log(f"AVISO: a regra '{nome_regra}' tem Novo_Val que não cabe em {tipo_saida.name} "
                    f"({', '.join(f'{v:g}' for v in fora[:8])}): escolha outro tipo de dado de saída.")
        if nodata_saida in novos:
            log(f"AVISO: na regra '{nome_regra}', Novo_Val = {nodata_saida:g} é também o valor de 'sem dado' "
                f"da saída; esses pixels serão lidos como sem dado.")

    total = len(arquivos)
    total_progresso = total * 100
    if progress_callback:
        progress_callback(0, total_progresso, "Aguardando início da reclassificação")

    def tarefa(item):
        caminho, nome_regra = item
        checar_parada(stop_event)
        regras = regras_por_nome.get(nome_regra)
        if not regras:
            return ("erro", None, f"regra '{nome_regra}' não encontrada ou vazia")
        try:
            saida = pasta_saida_path / _subpasta_relativa(caminho, pasta_base) / f"{prefixo}{Path(caminho).name}"
            presos = reclassificar_raster_arquivo(
                caminho, regras, str(saida), dtype=dtype, outros=outros, direita=direita,
                incluir_minimo=incluir_minimo, sobrescrever=sobrescrever, stop_event=stop_event)
            return ("ok", saida, presos)
        except ProcessamentoInterrompido:
            raise
        except Exception as exc:
            return ("erro", None, str(exc))

    saidas: List[Path] = []
    erros: List[str] = []
    avisos_aux: List[str] = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futuros = {executor.submit(tarefa, item): item for item in arquivos}
        try:
            for n, futuro in enumerate(as_completed(futuros), start=1):
                if stop_event is not None and stop_event.is_set():
                    raise ProcessamentoInterrompido("Processamento interrompido pelo usuário.")
                caminho, nome_regra = futuros[futuro]
                status, saida, erro = futuro.result()
                nome = Path(caminho).name
                if status == "ok":
                    saidas.append(saida)
                    log(f"[{n}/{total}] {nome} (regra: {nome_regra}) -> OK: {os.path.relpath(saida, pasta_saida_path)}")
                    if erro:
                        avisos_aux.append(nome)
                        log(f"  AVISO: não foi possível apagar a tabela de atributos antiga ({', '.join(erro)}); "
                            "o ArcGIS pode mostrar classes antigas. Feche o raster no ArcGIS e apague esses arquivos.")
                else:
                    erros.append(f"{nome}: {erro}")
                    log(f"[{n}/{total}] {nome} (regra: {nome_regra}) -> ERRO: {erro}")
                if progress_callback:
                    progress_callback(n * 100, total_progresso, f"Rasters concluídos: {n}/{total}")
        except ProcessamentoInterrompido:
            for f in futuros:
                f.cancel()
            raise
    return {"saidas": saidas, "erros": erros, "total": total, "auxiliares_presos": avisos_aux}


# =============================================================================
# Processamento 7: recorte + métricas FUNDIDO em memória (sem salvar recorte)
# =============================================================================

def _recorte_band1(src, geom, all_touched: bool):
    """Recorta a banda 1 pelo polígono. Retorna (array2d, transform) ou None."""
    try:
        arr, tr = mask(src, [mapping(geom)], crop=True, all_touched=all_touched, filled=True)
    except ValueError:
        return None
    if arr.size == 0:
        return None
    return arr[0], tr


def _reprojetar_band1_para_utm(banda1, transform_clipado, src_crs, src_nodata, crs_utm, num_threads: int = 1):
    """Reprojeta a banda 1 recortada para UTM. Retorna (array2d, res_x, res_y, nodata)."""
    altura, largura = banda1.shape
    bounds = array_bounds(altura, largura, transform_clipado)
    dst_transform, dst_w, dst_h = calculate_default_transform(src_crs, crs_utm, largura, altura, *bounds)
    nodata_val = src_nodata if src_nodata is not None else 0
    destino = np.full((dst_h, dst_w), nodata_val, dtype=banda1.dtype)
    reproject(
        source=banda1,
        destination=destino,
        src_transform=transform_clipado,
        src_crs=src_crs,
        src_nodata=nodata_val,
        dst_transform=dst_transform,
        dst_crs=crs_utm,
        dst_nodata=nodata_val,
        resampling=Resampling.nearest,
        num_threads=max(1, int(num_threads)),
    )
    return destino, abs(dst_transform.a), abs(dst_transform.e), nodata_val


def _caminho_unico_fs(pasta: Path, base: str, extensao: str) -> Path:
    """Nome de arquivo único baseado no sistema de arquivos (seguro entre processos)."""
    extensao = extensao if extensao.startswith(".") else f".{extensao}"
    caminho = pasta / f"{base}{extensao}"
    if not caminho.exists():
        return caminho
    contador = 2
    while True:
        candidato = pasta / f"{base}_{contador}{extensao}"
        if not candidato.exists():
            return candidato
        contador += 1


def _crs_e_metrico(src_crs) -> bool:
    """Indica se o CRS é projetado em metros (para permitir pular a reprojeção)."""
    try:
        if not src_crs.is_projected:
            return False
        unidade = (src_crs.linear_units or "").lower()
        return ("metre" in unidade) or ("meter" in unidade) or (unidade == "m")
    except Exception:
        return False


def _nome_partial_raster(raster_path: str) -> str:
    """
    Nome determinístico do CSV parcial de um raster, usado pelo recurso de retomada.

    Combina o nome do arquivo com um hash do caminho completo, para que dois rasters
    com o mesmo nome em pastas diferentes não colidam no cache de retomada.
    """
    stem = normalizar_nome(Path(raster_path).stem, max_len=80)
    h = hashlib.md5(str(Path(raster_path).resolve()).encode("utf-8")).hexdigest()[:8]
    return f"{stem}__{h}.csv"


def _assinatura_config_fundido(
    shp_buffers_lista, niveis_ativos, metricas_norm, vizinhanca,
    campo_id, campo_diametro, all_touched, pular_reproj_metrico,
) -> str:
    """
    Assinatura (hash) das opções de processamento. Se a configuração mudar entre
    execuções, o cache de retomada é descartado para não misturar resultados.
    """
    config = {
        "buffers": sorted(Path(s).name for s in shp_buffers_lista),
        "niveis": sorted(niveis_ativos),
        "metricas": {k: sorted(v) for k, v in sorted(metricas_norm.items())},
        "vizinhanca": str(vizinhanca),
        "campo_id": campo_id,
        "campo_diametro": campo_diametro,
        "all_touched": bool(all_touched),
        "pular_reproj_metrico": bool(pular_reproj_metrico),
    }
    bruto = json.dumps(config, sort_keys=True, ensure_ascii=False)
    return hashlib.md5(bruto.encode("utf-8")).hexdigest()


def _escanear_stems_csv_existentes(pasta_saida_path: Path) -> set:
    """
    Varre a pasta de saída e TODAS as subpastas, devolvendo o conjunto de nomes
    (sem extensão, minúsculos) dos CSVs já gerados. Usado para retomar a partir
    dos arquivos realmente existentes no disco (independe da versão que os criou).
    Ignora o CSV consolidado e arquivos auxiliares.
    """
    stems = set()
    try:
        for p in pasta_saida_path.rglob("*.csv"):
            nome = p.stem.lower()
            if nome.startswith("metricas_fundido_geral"):
                continue
            stems.add(nome)
    except Exception:
        pass
    return stems


def _worker_fundido_raster(tarefa: dict, stop_flag) -> dict:
    """
    Processa UM raster contra as feições recebidas (uma tarefa por raster).

    Executa em PROCESSO separado (ProcessPoolExecutor) ou em thread (fallback).
    Quando "salvar CSV por feição" está ligado, grava cada CSV de forma ATÔMICA
    (escreve .tmp e renomeia) e NÃO acumula DataFrames em memória — o consolidado
    é montado depois a partir dos CSVs no disco. Isso reduz muito o uso de memória
    (menos risco de o processo ser encerrado por falta de RAM) e permite retomar.
    """
    import pylandstats as pls

    raster_path = tarefa["raster_path"]
    ano = tarefa["ano"]
    niveis_ativos = tarefa["niveis_ativos"]
    metricas_por_nivel = tarefa["metricas_por_nivel"]
    vizinhanca = str(tarefa["vizinhanca"])
    all_touched = tarefa["all_touched"]
    salvar_csv_por_feicao = tarefa["salvar_csv_por_feicao"]
    salvar_parquet_por_ano = tarefa.get("salvar_parquet_por_ano", False)
    parquet_target = Path(tarefa["parquet_target"]) if tarefa.get("parquet_target") else None
    pular_reproj_metrico = tarefa["pular_reproj_metrico"]
    reproject_threads = tarefa.get("reproject_threads", 1)
    pasta_saida = Path(tarefa["pasta_saida"])
    partial_target = Path(tarefa["partial_target"])
    progress_q = tarefa.get("progress_q")

    mensagens: List[str] = []
    n_ok = 0
    n_pulados = 0
    n_csv = 0
    dfs_raster: List["pd.DataFrame"] = []  # usado SÓ quando não salva por feição
    pq_writer = None
    pq_schema = None
    pq_cols = None
    pq_tmp = None
    if salvar_parquet_por_ano and parquet_target is not None:
        import pyarrow as pa  # noqa
        import pyarrow.parquet as pq  # noqa

    def parar() -> bool:
        return stop_flag is not None and stop_flag.is_set()

    try:
        with rasterio.Env(GDAL_CACHEMAX=128):
            with rasterio.open(raster_path) as src:
                src_crs = src.crs
                src_nodata = src.nodata
                usar_nativo = bool(pular_reproj_metrico and _crs_e_metrico(src_crs))

                for feat in tarefa["features"]:
                    if parar():
                        return {"status": "interrompido", "raster": Path(raster_path).name}

                    # Reporta o avanço (1 feição) para a barra de progresso do app.
                    if progress_q is not None:
                        try:
                            progress_q.put(1)
                        except Exception:
                            pass

                    geom = feat["geom"]
                    if geom is None or getattr(geom, "is_empty", False):
                        continue

                    rec = _recorte_band1(src, geom, all_touched)
                    if rec is None:
                        n_pulados += 1
                        continue
                    banda1, transform_clipado = rec

                    if usar_nativo:
                        array2d = banda1
                        res_x = abs(transform_clipado.a)
                        res_y = abs(transform_clipado.e)
                        nodata_val = src_nodata if src_nodata is not None else 0
                    else:
                        try:
                            crs_utm = CRS.from_epsg(int(feat["utm_epsg"]))
                        except Exception:
                            crs_utm = CRS.from_epsg(32700)
                        array2d, res_x, res_y, nodata_val = _reprojetar_band1_para_utm(
                            banda1, transform_clipado, src_crs, src_nodata, crs_utm, num_threads=reproject_threads
                        )

                    try:
                        ls = pls.Landscape(array2d, res=(res_x, res_y), nodata=nodata_val, neighborhood_rule=vizinhanca)
                    except Exception as exc:
                        mensagens.append(f"Landscape id={feat['id_token']} ano={ano}: {exc}")
                        continue

                    dfs_niveis = []
                    for nivel in niveis_ativos:
                        metricas = metricas_por_nivel.get(nivel) or None
                        try:
                            if nivel == "patch":
                                df = ls.compute_patch_metrics_df(metrics=metricas)
                            elif nivel == "class":
                                df = ls.compute_class_metrics_df(metrics=metricas)
                            else:
                                df = ls.compute_landscape_metrics_df(metrics=metricas)
                            df = preparar_df_metricas(df, nivel=nivel, raster_path=Path(raster_path))
                            df["year"] = ano
                            df["neighborhood_rule"] = vizinhanca
                            df.insert(0, "diametro", feat["diam_raw"])
                            df.insert(0, "id", feat["id_raw"])
                            dfs_niveis.append(df)
                        except Exception as exc:
                            mensagens.append(f"{nivel} id={feat['id_token']} ano={ano}: {exc}")

                    if not dfs_niveis:
                        continue

                    df_clip = pd.concat(dfs_niveis, ignore_index=True, sort=False)
                    n_ok += 1

                    if salvar_csv_por_feicao:
                        pasta_ano = pasta_saida / normalizar_nome(ano, max_len=20)
                        pasta_ano.mkdir(parents=True, exist_ok=True)
                        base = normalizar_nome(f"id_{feat['id_token']}_diam_{feat['diam_token']}_{ano}", max_len=150)
                        caminho_csv = pasta_ano / f"{base}.csv"
                        # Escrita atômica: o CSV final só aparece completo (retomada confiável).
                        tmp = caminho_csv.with_suffix(".csv.tmp")
                        df_clip.to_csv(tmp, index=False, encoding="utf-8-sig")
                        os.replace(tmp, caminho_csv)
                        n_csv += 1

                    if salvar_parquet_por_ano and parquet_target is not None:
                        import pyarrow as pa
                        import pyarrow.parquet as pq
                        if pq_writer is None:
                            pq_cols = list(df_clip.columns)
                            tabela = pa.Table.from_pandas(df_clip, preserve_index=False)
                            pq_schema = tabela.schema
                            parquet_target.parent.mkdir(parents=True, exist_ok=True)
                            pq_tmp = parquet_target.with_suffix(".parquet.tmp")
                            pq_writer = pq.ParquetWriter(str(pq_tmp), pq_schema, compression="zstd")
                            pq_writer.write_table(tabela)
                        else:
                            df_alinhado = df_clip.reindex(columns=pq_cols)
                            tabela = pa.Table.from_pandas(df_alinhado, preserve_index=False)
                            if tabela.schema != pq_schema:
                                tabela = tabela.cast(pq_schema, safe=False)
                            pq_writer.write_table(tabela)

                    if (not salvar_csv_por_feicao) and (not salvar_parquet_por_ano):
                        dfs_raster.append(df_clip)
    except Exception as exc:
        if pq_writer is not None:
            try:
                pq_writer.close()
            except Exception:
                pass
        return {"status": "erro_raster", "raster": Path(raster_path).name, "erro": str(exc), "mensagens": mensagens}

    parquet_path = None
    if salvar_parquet_por_ano and parquet_target is not None and pq_writer is not None:
        pq_writer.close()
        os.replace(pq_tmp, parquet_target)
        parquet_path = str(parquet_target)

    partial_path = None
    if (not salvar_csv_por_feicao) and (not salvar_parquet_por_ano) and dfs_raster:
        df_raster = pd.concat(dfs_raster, ignore_index=True, sort=False)
        partial_target.parent.mkdir(parents=True, exist_ok=True)
        tmp = partial_target.with_suffix(".csv.tmp")
        df_raster.to_csv(tmp, index=False, encoding="utf-8-sig")
        os.replace(tmp, partial_target)
        partial_path = str(partial_target)

    return {
        "status": "ok",
        "raster": Path(raster_path).name,
        "ano": ano,
        "n_ok": n_ok,
        "n_pulados": n_pulados,
        "n_csv": n_csv,
        "partial_path": partial_path,
        "parquet_path": parquet_path,
        "mensagens": mensagens,
    }


def processar_buffers_rasters_em_memoria(
    shp_buffers_lista: List[str],
    rasters_lista: List[str],
    pasta_saida: str,
    niveis: Dict[str, bool],
    metricas_por_nivel: Dict[str, Sequence[str]],
    vizinhanca: str = "8",
    campo_id: str = "id_pt",
    campo_diametro: str = "diam_m",
    all_touched: bool = True,
    max_workers: Optional[int] = None,
    usar_processos: bool = True,
    salvar_csv_por_feicao: bool = True,
    salvar_parquet_por_ano: bool = False,
    gerar_parquet_unico_final: bool = False,
    pular_reproj_metrico: bool = False,
    retomar: bool = True,
    log_callback=None,
    load_callback=None,
    progress_callback=None,
    stop_event: Optional[threading.Event] = None,
) -> Dict[str, object]:
    """
    Recorte + métricas EM MEMÓRIA (sem salvar recortes), com paralelismo e RETOMADA.

    Retomada (corrigida): baseada nos CSVs realmente existentes na pasta de saída e
    em TODAS as subpastas. Cada feição cujo CSV (id_<ID>_diam_<DIAM>_<ANO>.csv) já
    existe é pulada. A barra de progresso, a contagem de arquivos e a estimativa de
    tempo refletem o estado GLOBAL (já existentes + novos).
    """
    def log(msg: str):
        if log_callback:
            log_callback(msg)

    def set_load(acao: str, mensagem: str = ""):
        if load_callback:
            load_callback(acao, mensagem)

    def set_progress(valor: float, maximo: float, mensagem: str = ""):
        if progress_callback:
            progress_callback(valor, maximo, mensagem)

    if not shp_buffers_lista:
        raise ValueError("Selecione pelo menos um shapefile de buffers.")
    if not rasters_lista:
        raise ValueError("Selecione pelo menos um raster.")
    if not pasta_saida:
        raise ValueError("Selecione a pasta de saída dos CSVs.")
    if vizinhanca not in {"4", "8"}:
        raise ValueError("A vizinhança deve ser '4' ou '8'.")
    if not any(niveis.values()):
        raise ValueError("Selecione pelo menos um nível escalar: patch, classe ou paisagem.")
    if not (salvar_csv_por_feicao or salvar_parquet_por_ano):
        raise ValueError("Selecione pelo menos uma saída no Fundido: CSV por feição e/ou Parquet por ano.")

    # Parquet exige pyarrow; o Parquet único final é montado depois (DuckDB/pyarrow).
    if salvar_parquet_por_ano and not _modulo_disponivel("pyarrow"):
        if salvar_csv_por_feicao:
            log("AVISO: 'pyarrow' não instalado — Parquet por ano desativado; mantendo apenas CSV por feição. "
                "Instale com: pip install pyarrow")
            salvar_parquet_por_ano = False
            gerar_parquet_unico_final = False
        else:
            raise ValueError("Saída Parquet selecionada, mas 'pyarrow' não está instalado. "
                             "Instale 'pip install pyarrow' ou marque 'Salvar CSV por feição'.")
    if gerar_parquet_unico_final and not salvar_parquet_por_ano:
        gerar_parquet_unico_final = False

    if not max_workers or max_workers < 1:
        max_workers = numero_padrao_de_nucleos()
    cpus = os.cpu_count() or 4
    if usar_processos:
        max_workers = max(1, min(max_workers, cpus))

    pasta_saida_path = Path(pasta_saida)
    pasta_saida_path.mkdir(parents=True, exist_ok=True)

    # Limpeza segura de temporários ÓRFÃOS da própria aplicação (não toca em %TEMP%).
    for antigo in pasta_saida_path.glob("__tmp_fundido_*"):
        if antigo.is_dir():
            shutil.rmtree(antigo, ignore_errors=True)
    # Pasta temporária só para parciais quando NÃO se salva CSV por feição.
    pasta_temp = pasta_saida_path / f"__tmp_fundido_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    pasta_temp.mkdir(parents=True, exist_ok=True)

    set_load("start", "Importando PyLandStats")
    importar_pylandstats()
    set_load("stop", "PyLandStats importado")

    niveis_ativos = [n for n, ativo in niveis.items() if ativo]
    metricas_norm = {nivel: list(metricas_por_nivel.get(nivel, []) or []) for nivel in niveis_ativos}

    # ---- Carrega e valida os shapefiles de buffers ----
    buffers_info = []
    total_features = 0
    for idx_shp, shp in enumerate(shp_buffers_lista, start=1):
        checar_parada(stop_event)
        set_load("start", f"Carregando buffers {idx_shp}/{len(shp_buffers_lista)}: {Path(shp).name}")
        try:
            gdf = corrigir_geometrias(ler_shapefile(shp))
        finally:
            set_load("stop", f"Buffers carregados: {Path(shp).name}")

        if not any("Polygon" in t for t in set(gdf.geom_type.unique())):
            raise ValueError(f"O shapefile {Path(shp).name} não contém polígonos/buffers.")
        if campo_id not in gdf.columns:
            raise ValueError(f"O campo '{campo_id}' não existe em {Path(shp).name}.")
        if campo_diametro not in gdf.columns:
            raise ValueError(f"O campo '{campo_diametro}' não existe em {Path(shp).name}.")

        utm = calcular_utm_por_feature(gdf)
        rotulos = {}
        for idx_feature in gdf.index:
            row = gdf.loc[idx_feature]
            crs_utm = utm.get(idx_feature)
            if crs_utm is None:
                crs_utm = obter_crs_utm_por_geometria(gdf.loc[[idx_feature]])
            rotulos[idx_feature] = {
                "id_raw": row[campo_id],
                "diam_raw": row[campo_diametro],
                "id_token": formatar_valor_atributo(row[campo_id]),
                "diam_token": formatar_valor_atributo(row[campo_diametro]),
                "utm_epsg": int(crs_utm.to_epsg() or 32700),
            }
        buffers_info.append({"gdf": gdf, "rotulos": rotulos})
        total_features += len(gdf)
        log(f"Buffers '{Path(shp).name}': {len(gdf)} feições (campos {campo_id}, {campo_diametro} OK).")

    # ---- Valida rasters e extrai o ano (obrigatório) ----
    rasters_validos = []
    for idx, raster_path in enumerate(rasters_lista, start=1):
        checar_parada(stop_event)
        set_load("start", f"Verificando raster {idx}/{len(rasters_lista)}: {Path(raster_path).name}")
        try:
            ano = extrair_ano_raster(raster_path)
            if ano is None:
                log(f"AVISO: ano não encontrado em '{Path(raster_path).name}'. Raster ignorado.")
                continue
            with rasterio.open(raster_path) as src:
                if src.crs is None:
                    raise ValueError("raster sem CRS definido")
                rasters_validos.append((raster_path, str(src.crs), ano))
                log(f"Raster validado: {Path(raster_path).name} | ano={ano} | CRS={src.crs}")
        except Exception as exc:
            log(f"Raster ignorado (inválido): {Path(raster_path).name} -> {exc}")
        finally:
            set_load("stop", f"Raster verificado: {Path(raster_path).name}")

    if not rasters_validos:
        raise ValueError("Nenhum raster válido (com ano no nome e CRS) foi encontrado.")

    # ---- Controle de retomada (assinatura de configuração) ----
    resume_efetivo = bool(retomar and (salvar_csv_por_feicao or salvar_parquet_por_ano))
    if retomar and not (salvar_csv_por_feicao or salvar_parquet_por_ano):
        log("Aviso: a retomada requer CSV por feição ou Parquet por ano. Será ignorada nesta execução.")

    arquivo_config = pasta_saida_path / "__fundido_config.json"
    assinatura_atual = _assinatura_config_fundido(
        shp_buffers_lista, niveis_ativos, metricas_norm, vizinhanca,
        campo_id, campo_diametro, all_touched, pular_reproj_metrico,
    )
    if resume_efetivo and arquivo_config.exists():
        try:
            cfg_antiga = json.loads(arquivo_config.read_text(encoding="utf-8"))
        except Exception:
            cfg_antiga = {}
        if cfg_antiga.get("assinatura") != assinatura_atual:
            resume_efetivo = False
            log("Configuração mudou desde a última execução (níveis/métricas/vizinhança/campos/buffers). "
                "A retomada foi DESATIVADA para não misturar resultados: tudo será recalculado e sobrescrito.")
    try:
        arquivo_config.write_text(
            json.dumps({"assinatura": assinatura_atual, "atualizado": datetime.now().isoformat(timespec="seconds")},
                       ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
    except Exception:
        pass

    # Retomada por feição (modo CSV puro) usa os CSVs existentes; modo Parquet usa
    # a existência do Parquet do ano (por raster).
    resume_por_feicao = resume_efetivo and salvar_csv_por_feicao and not salvar_parquet_por_ano
    if resume_por_feicao:
        set_load("start", "Lendo CSVs já existentes (retomada)")
        stems_existentes = _escanear_stems_csv_existentes(pasta_saida_path)
        set_load("stop", "Leitura de CSVs existentes concluída")
        log(f"Retomada (CSV): {len(stems_existentes)} CSV(s) já existentes encontrados.")
    else:
        stems_existentes = set()

    # ---- Cache de geometrias reprojetadas por CRS de raster, por shapefile ----
    geoms_cache: Dict[tuple, Dict[object, object]] = {}

    def geoms_para(idx_shp: int, crs_str: str) -> Dict[object, object]:
        chave = (idx_shp, crs_str)
        if chave not in geoms_cache:
            gdf = buffers_info[idx_shp]["gdf"]
            gdf_crs = gdf.to_crs(crs_str)
            geoms_cache[chave] = {idx: gdf_crs.geometry.loc[idx] for idx in gdf.index}
        return geoms_cache[chave]

    def parquet_alvo(raster_path: str, ano: str) -> Path:
        return (pasta_saida_path / normalizar_nome(ano, max_len=20)
                / (normalizar_nome(Path(raster_path).stem, max_len=120) + ".parquet"))

    # ---- Monta uma tarefa por raster, pulando o que já está pronto ----
    reproject_threads = 1
    tarefas = []
    unidades_ja_feitas = 0
    rasters_totalmente_prontos = 0
    for (raster_path, crs_str, ano) in rasters_validos:
        pq_target = parquet_alvo(raster_path, ano)
        # Retomada por raster (modo Parquet): pula raster cujo Parquet do ano já existe.
        if resume_efetivo and salvar_parquet_por_ano and pq_target.exists():
            rasters_totalmente_prontos += 1
            unidades_ja_feitas += total_features
            log(f"Retomada (Parquet): '{Path(raster_path).name}' já tem Parquet do ano — pulando.")
            continue
        features = []
        for idx_shp, info in enumerate(buffers_info):
            geoms = geoms_para(idx_shp, crs_str)
            for idx_feature in info["gdf"].index:
                rot = info["rotulos"][idx_feature]
                if resume_por_feicao:
                    base = normalizar_nome(f"id_{rot['id_token']}_diam_{rot['diam_token']}_{ano}", max_len=150).lower()
                    if base in stems_existentes:
                        unidades_ja_feitas += 1
                        continue
                features.append({
                    "geom": geoms.get(idx_feature),
                    "id_raw": rot["id_raw"],
                    "diam_raw": rot["diam_raw"],
                    "id_token": rot["id_token"],
                    "diam_token": rot["diam_token"],
                    "utm_epsg": rot["utm_epsg"],
                })
        if not features:
            rasters_totalmente_prontos += 1
            continue
        tarefas.append({
            "raster_path": raster_path,
            "ano": ano,
            "features": features,
            "niveis_ativos": niveis_ativos,
            "metricas_por_nivel": metricas_norm,
            "vizinhanca": vizinhanca,
            "all_touched": all_touched,
            "salvar_csv_por_feicao": salvar_csv_por_feicao,
            "salvar_parquet_por_ano": salvar_parquet_por_ano,
            "parquet_target": str(pq_target),
            "pular_reproj_metrico": pular_reproj_metrico,
            "reproject_threads": reproject_threads,
            "pasta_saida": str(pasta_saida_path),
            "partial_target": str(pasta_temp / _nome_partial_raster(raster_path)),
        })

    total_global = max(1, len(rasters_validos) * total_features)
    if resume_efetivo:
        log(f"Retomada: {unidades_ja_feitas} feição(ões) já processada(s) | "
            f"{rasters_totalmente_prontos} raster(s) já prontos serão pulados.")
    if not tarefas:
        log("Nada novo a processar — tudo já existe. Será (re)montado o Parquet único, se solicitado.")

    modo = "processos" if usar_processos else "threads"
    log(f"Rasters a processar agora: {len(tarefas)} | feições já feitas: {unidades_ja_feitas} | "
        f"total global: {total_global} | núcleos: {max_workers} | modo: {modo}")
    log(f"Saídas: CSV por feição={'sim' if salvar_csv_por_feicao else 'não'} | "
        f"Parquet por ano={'sim' if salvar_parquet_por_ano else 'não'} | "
        f"Parquet único final={'sim' if gerar_parquet_unico_final else 'não'} | "
        f"retomada={'sim' if resume_efetivo else 'não'}")

    partial_paths: List[str] = []
    parquet_paths: List[str] = []
    anos_set = set()
    n_csv_total = 0
    n_parquet_ano = 0

    def consumir(resultado: dict):
        nonlocal n_csv_total, n_parquet_ano
        for m in resultado.get("mensagens", []):
            log(f"  -> {m}")
        if resultado.get("status") == "erro_raster":
            log(f"ERRO no raster {resultado.get('raster')}: {resultado.get('erro')}")
            return
        if resultado.get("status") != "ok":
            return
        if resultado.get("partial_path"):
            partial_paths.append(resultado["partial_path"])
        if resultado.get("parquet_path"):
            parquet_paths.append(resultado["parquet_path"])
            n_parquet_ano += 1
        anos_set.add(resultado["ano"])
        n_csv_total += resultado.get("n_csv", 0)
        log(f"Raster {resultado['raster']} (ano {resultado['ano']}): feições OK={resultado['n_ok']}, "
            f"sem sobreposição={resultado['n_pulados']}, CSV={resultado['n_csv']}, "
            f"Parquet ano={'sim' if resultado.get('parquet_path') else 'não'}")

    # ---- Progresso GLOBAL (já feitas + novas) + estimativa de tempo ----
    set_progress(unidades_ja_feitas, total_global,
                 f"Feições {unidades_ja_feitas}/{total_global} (já existiam: {unidades_ja_feitas}) | preparando...")
    inicio_proc = time.time()
    feitas = {"n": unidades_ja_feitas}
    base_ja = unidades_ja_feitas
    parar_poll = threading.Event()

    def iniciar_poller(fila):
        def poller():
            while not parar_poll.is_set():
                drenadas = 0
                try:
                    while True:
                        fila.get_nowait()
                        drenadas += 1
                except queue.Empty:
                    pass
                except Exception:
                    pass
                if drenadas:
                    feitas["n"] += drenadas
                n = min(feitas["n"], total_global)
                novas = max(0, n - base_ja)
                decorrido = time.time() - inicio_proc
                if novas > 0 and decorrido > 0:
                    taxa = novas / decorrido
                    restante = max(0, total_global - n)
                    eta = restante / taxa if taxa > 0 else 0.0
                    msg = (f"Feições {n}/{total_global} (já existiam: {base_ja}) | "
                           f"tempo restante estimado ~ {formatar_duracao_horas(eta)}")
                else:
                    msg = f"Feições {n}/{total_global} (já existiam: {base_ja}) | estimando tempo restante..."
                set_progress(n, total_global, msg)
                time.sleep(0.3)
        th = threading.Thread(target=poller, daemon=True)
        th.start()
        return th

    def tratar_futuro(futuro, futuros):
        """Lê o resultado de um futuro com robustez (um raster com erro não derruba o lote)."""
        try:
            return futuro.result()
        except ProcessamentoInterrompido:
            raise
        except BrokenProcessPool:
            raise RuntimeError(
                "Um processo trabalhador foi encerrado inesperadamente (provável falta de memória "
                "ou falha do GDAL/PyLandStats). Os CSVs já gerados foram preservados — rode novamente "
                "com 'Retomar de onde parou' para continuar. Dica: reduza o número de núcleos ou "
                "desmarque 'Usar processos'."
            )
        except Exception as exc:
            nome = Path(futuros[futuro]["raster_path"]).name
            log(f"ERRO no raster {nome}: {exc}. Continuando com os demais rasters.")
            return None

    # ---- Execução paralela ----
    if usar_processos:
        manager = multiprocessing.Manager()
        mp_stop = manager.Event()
        progress_q = manager.Queue()
        for t in tarefas:
            t["progress_q"] = progress_q
        parar_watcher = threading.Event()

        def watcher():
            while not parar_watcher.is_set():
                if stop_event is not None and stop_event.is_set():
                    mp_stop.set()
                    return
                time.sleep(0.2)

        thread_watcher = threading.Thread(target=watcher, daemon=True)
        thread_watcher.start()
        poller_thread = iniciar_poller(progress_q)
        try:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futuros = {executor.submit(_worker_fundido_raster, t, mp_stop): t for t in tarefas}
                try:
                    for futuro in as_completed(futuros):
                        if stop_event is not None and stop_event.is_set():
                            mp_stop.set()
                            for f in futuros:
                                f.cancel()
                            raise ProcessamentoInterrompido("Processamento interrompido pelo usuário.")
                        resultado = tratar_futuro(futuro, futuros)
                        if resultado is None:
                            continue
                        if resultado.get("status") == "interrompido":
                            raise ProcessamentoInterrompido("Processamento interrompido pelo usuário.")
                        consumir(resultado)
                except ProcessamentoInterrompido:
                    mp_stop.set()
                    for f in futuros:
                        f.cancel()
                    raise
        finally:
            parar_poll.set()
            parar_watcher.set()
            manager.shutdown()
    else:
        progress_q = queue.Queue()
        for t in tarefas:
            t["progress_q"] = progress_q
        poller_thread = iniciar_poller(progress_q)
        try:
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futuros = {executor.submit(_worker_fundido_raster, t, stop_event): t for t in tarefas}
                try:
                    for futuro in as_completed(futuros):
                        if stop_event is not None and stop_event.is_set():
                            for f in futuros:
                                f.cancel()
                            raise ProcessamentoInterrompido("Processamento interrompido pelo usuário.")
                        resultado = tratar_futuro(futuro, futuros)
                        if resultado is None:
                            continue
                        if resultado.get("status") == "interrompido":
                            raise ProcessamentoInterrompido("Processamento interrompido pelo usuário.")
                        consumir(resultado)
                except ProcessamentoInterrompido:
                    for f in futuros:
                        f.cancel()
                    raise
        finally:
            parar_poll.set()

    set_progress(total_global, total_global, "Processamento concluído.")

    anos_finais = set(anos_set)

    # ---- Parquet ÚNICO final (junta os Parquets por ano) ----
    parquet_unico_final = None
    if gerar_parquet_unico_final:
        # Junta TODOS os Parquets por ano existentes na saída (inclui os retomados),
        # exceto Parquets "únicos" anteriores.
        todos_pq = [p for p in pasta_saida_path.rglob("*.parquet")
                    if not p.name.lower().startswith("metricas_fundido_unico")]
        if todos_pq:
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            parquet_unico_final = pasta_saida_path / f"metricas_fundido_UNICO_{ts}.parquet"
            set_load("start", "Gerando Parquet único final")
            try:
                _juntar_parquets_em_unico([str(p) for p in todos_pq], parquet_unico_final, "zstd", log)
                log(f"Parquet único final: {parquet_unico_final} (de {len(todos_pq)} Parquet(s) por ano)")
            except Exception as exc:
                log(f"AVISO: falha ao gerar o Parquet único final: {exc}")
                parquet_unico_final = None
            finally:
                set_load("stop", "Parquet único final concluído")

    # ---- Limpeza dos temporários da própria aplicação ----
    shutil.rmtree(pasta_temp, ignore_errors=True)

    feicoes_concluidas = min(feitas["n"], total_global)
    return {
        "n_csvs_por_feicao": n_csv_total,
        "csv_geral": None,
        "n_parquet_ano": n_parquet_ano,
        "parquet_unico": str(parquet_unico_final) if parquet_unico_final else None,
        "anos": sorted(anos_finais),
        "rasters_processados": len(tarefas),
        "rasters_pulados": rasters_totalmente_prontos,
        "rasters_total": len(rasters_validos),
        "feicoes_ja_existentes": unidades_ja_feitas,
        "feicoes_concluidas": feicoes_concluidas,
        "total_global": total_global,
        "tarefas": len(tarefas),
    }


# =============================================================================
# Processamento 8: consolidação rápida das tabelas (CSV / Parquet)
# =============================================================================

def _juntar_parquets_em_unico(parquets: List[str], destino: Path, compressao: str = "zstd", log=None) -> None:
    """
    Junta vários arquivos Parquet num único Parquet. Usa DuckDB se disponível
    (rápido, baixa memória); senão usa pyarrow (lê e regrava em fluxo).
    """
    def _log(m):
        if log:
            log(m)

    destino = Path(destino)
    destino.parent.mkdir(parents=True, exist_ok=True)
    if not parquets:
        raise ValueError("Nenhum Parquet de entrada para juntar.")

    if _modulo_disponivel("duckdb"):
        import duckdb
        con = duckdb.connect()
        try:
            lista = "[" + ",".join("'" + p.replace("'", "''") + "'" for p in parquets) + "]"
            con.execute(
                f"COPY (SELECT * FROM read_parquet({lista}, union_by_name=true)) "
                f"TO '{str(destino.as_posix()).replace(chr(39), chr(39)*2)}' "
                f"(FORMAT PARQUET, COMPRESSION '{compressao}')"
            )
        finally:
            con.close()
        return

    # Alternativa com pyarrow.
    import pyarrow.parquet as pq
    writer = None
    schema = None
    tmp = destino.with_suffix(".parquet.tmp")
    try:
        for caminho in parquets:
            tabela = pq.read_table(caminho)
            if writer is None:
                schema = tabela.schema
                writer = pq.ParquetWriter(str(tmp), schema, compression=compressao)
                writer.write_table(tabela)
            else:
                if tabela.schema != schema:
                    tabela = tabela.cast(schema, safe=False)
                writer.write_table(tabela)
    finally:
        if writer is not None:
            writer.close()
    os.replace(tmp, destino)


def _modulo_disponivel(nome: str) -> bool:
    import importlib.util
    return importlib.util.find_spec(nome) is not None


def _ler_header_e_corpo(caminho: Path):
    """Lê um CSV em bytes e separa (header_bytes_incl_quebra, corpo_bytes)."""
    dados = caminho.read_bytes()
    pos = dados.find(b"\n")
    if pos == -1:
        return dados, b""  # arquivo só com cabeçalho, sem quebra
    return dados[: pos + 1], dados[pos + 1:]


def _header_normalizado(header_bytes: bytes) -> str:
    """Header como texto, sem BOM e sem quebras, para comparação de esquema."""
    txt = header_bytes.decode("utf-8-sig", errors="replace")
    return txt.replace("\r", "").replace("\n", "").strip()


def _consolidar_com_duckdb(
    pasta_entrada_path: Path,
    pasta_saida_path: Path,
    nome_base: str,
    gerar_csv: bool,
    gerar_parquet_unico: bool,
    gerar_parquet_por_ano: bool,
    coluna_particao: str,
    compressao: str,
    max_workers: int,
    log,
    set_load,
    set_progress,
    stop_event,
    entrada_parquet: bool = False,
) -> Dict[str, object]:
    """
    Consolida com DuckDB mostrando NÚMERO DE ARQUIVOS e ESTIMATIVA DE TEMPO.

    Estratégia: lista os arquivos (glob do DuckDB, rápido), guarda a lista numa
    tabela em disco e faz a leitura em LOTES — assim a barra avança de verdade,
    com 'X/Total arquivos' e ETA. Depois grava CSV/Parquet único/Parquet por ano.
    """
    import duckdb

    def q(caminho: Path) -> str:
        return str(caminho.as_posix()).replace("'", "''")

    ext = "parquet" if entrada_parquet else "csv"
    leitor = "read_parquet" if entrada_parquet else "read_csv_auto"
    glob = q(pasta_entrada_path / "**" / f"*.{ext}")
    csv_out = pasta_saida_path / f"{nome_base}.csv"
    pq_out = pasta_saida_path / f"{nome_base}.parquet"
    pq_part_dir = pasta_saida_path / f"{nome_base}_por_ano"
    temp_db = pasta_saida_path / f"__consol_{datetime.now().strftime('%Y%m%d_%H%M%S')}.duckdb"

    con = duckdb.connect(str(temp_db))
    con.execute(f"PRAGMA threads={max(1, int(max_workers))}")

    saidas = {}
    total_linhas = 0
    try:
        # ---- 1) Listar arquivos (rápido, em C++) e remover saídas anteriores ----
        set_load("start", "Listando arquivos (pode demorar com muitos arquivos)")
        try:
            con.execute(f"CREATE TABLE _files AS SELECT file FROM glob('{glob}')")
            con.execute(
                "DELETE FROM _files WHERE "
                "lower(file) LIKE '%metricas_fundido_geral%' OR "
                "lower(file) LIKE '%metricas_fundido_unico%' OR "
                "lower(file) LIKE '%tabela_consolidada%' OR "
                "lower(file) LIKE '%__tmp_fundido%' OR "
                f"lower(file) LIKE '%{nome_base.lower()}%'"
            )
            total = con.execute("SELECT count(*) FROM _files").fetchone()[0]
        finally:
            set_load("stop", "Listagem concluída")

        if total == 0:
            raise ValueError(f"Nenhum {ext.upper()} encontrado na pasta de entrada (e subpastas).")
        log(f"Arquivos a consolidar: {total}")
        set_progress(0, total, f"0/{total} arquivos lidos | estimando tempo...")

        # ---- 2) Ler em LOTES para mostrar progresso + ETA ----
        inicio = time.time()
        lote = 5000
        primeiro = True
        lidos = 0
        offset = 0
        while offset < total:
            checar_parada(stop_event)
            grupo = [r[0] for r in con.execute(
                f"SELECT file FROM _files LIMIT {lote} OFFSET {offset}"
            ).fetchall()]
            if not grupo:
                break
            lista_sql = "[" + ",".join("'" + p.replace("'", "''") + "'" for p in grupo) + "]"
            fonte = f"{leitor}({lista_sql}, union_by_name=true)"
            try:
                if primeiro:
                    con.execute(f"CREATE TABLE _consol AS SELECT * FROM {fonte}")
                    primeiro = False
                else:
                    con.execute(f"INSERT INTO _consol BY NAME SELECT * FROM {fonte}")
            except Exception as exc:
                log(f"Aviso: falha ao ler um lote de {len(grupo)} arquivo(s): {exc}")

            offset += lote
            lidos = min(offset, total)
            decorrido = time.time() - inicio
            if lidos > 0 and decorrido > 0:
                taxa = lidos / decorrido
                eta = (total - lidos) / taxa if taxa > 0 else 0.0
                set_progress(lidos, total, f"Lendo {lidos}/{total} arquivos | tempo restante ~ {formatar_duracao_horas(eta)}")
            else:
                set_progress(lidos, total, f"Lendo {lidos}/{total} arquivos...")

        if primeiro:
            raise ValueError("Não foi possível ler nenhum arquivo (verifique os CSVs).")

        try:
            total_linhas = con.execute("SELECT count(*) FROM _consol").fetchone()[0]
        except Exception:
            total_linhas = 0
        log(f"Leitura concluída: {total} arquivos, {total_linhas} linhas.")
        set_progress(total, total, f"{total}/{total} arquivos lidos. Gravando saídas...")

        colunas = [r[0] for r in con.execute("DESCRIBE _consol").fetchall()]
        tem_coluna_part = coluna_particao in colunas

        # ---- 3) Gravar as saídas ----
        if gerar_csv:
            checar_parada(stop_event)
            set_load("start", "Gravando CSV consolidado")
            try:
                con.execute(f"COPY _consol TO '{q(csv_out)}' (FORMAT CSV, HEADER true, DELIMITER ',')")
                saidas["csv"] = csv_out
                log(f"CSV consolidado: {csv_out}")
            finally:
                set_load("stop", "CSV gravado")

        if gerar_parquet_unico:
            checar_parada(stop_event)
            set_load("start", "Gravando Parquet único")
            try:
                con.execute(f"COPY _consol TO '{q(pq_out)}' (FORMAT PARQUET, COMPRESSION '{compressao}')")
                saidas["parquet"] = pq_out
                log(f"Parquet único: {pq_out}")
            finally:
                set_load("stop", "Parquet único gravado")

        if gerar_parquet_por_ano:
            checar_parada(stop_event)
            if not tem_coluna_part:
                log(f"AVISO: coluna de partição '{coluna_particao}' não existe — Parquet por ano ignorado.")
            else:
                set_load("start", "Gravando Parquet por ano (particionado)")
                try:
                    con.execute(
                        f"COPY _consol TO '{q(pq_part_dir)}' "
                        f"(FORMAT PARQUET, PARTITION_BY ({coluna_particao}), COMPRESSION '{compressao}', OVERWRITE_OR_IGNORE true)"
                    )
                    saidas["parquet_por_ano"] = pq_part_dir
                    log(f"Parquet por ano: {pq_part_dir}")
                finally:
                    set_load("stop", "Parquet por ano gravado")
    finally:
        try:
            con.close()
        except Exception:
            pass
        # Remove o banco temporário do DuckDB.
        try:
            if temp_db.exists():
                temp_db.unlink()
            wal = temp_db.with_suffix(".duckdb.wal")
            if wal.exists():
                wal.unlink()
        except Exception:
            pass

    set_progress(max(1, total), max(1, total), "Consolidação (DuckDB) concluída.")
    return {"arquivos_lidos": total, "linhas": total_linhas, "saidas": saidas, "motor": "duckdb"}


def consolidar_tabelas(
    pasta_entrada: str,
    pasta_saida: str,
    nome_base: str = "tabela_consolidada",
    gerar_csv: bool = True,
    gerar_parquet_unico: bool = True,
    gerar_parquet_por_ano: bool = False,
    coluna_particao: str = "year",
    compressao: str = "zstd",
    max_workers: Optional[int] = None,
    entrada_parquet: bool = False,
    log_callback=None,
    load_callback=None,
    progress_callback=None,
    stop_event: Optional[threading.Event] = None,
) -> Dict[str, object]:
    """
    Consolida todos os CSVs de uma pasta (e subpastas) em tabela única.

    Motor preferencial: DuckDB (lê todos os CSVs de uma vez, em paralelo, e grava
    Parquet único + Parquet por ano + CSV). Se o DuckDB não estiver instalado,
    cai para um caminho alternativo: CSV por concatenação de bytes e Parquet via
    pyarrow (mais lento, mas funcional).
    """
    def log(msg: str):
        if log_callback:
            log_callback(msg)

    def set_load(acao: str, mensagem: str = ""):
        if load_callback:
            load_callback(acao, mensagem)

    def set_progress(valor: float, maximo: float, mensagem: str = ""):
        if progress_callback:
            progress_callback(valor, maximo, mensagem)

    if not pasta_entrada:
        raise ValueError("Selecione a pasta de entrada (com os subdiretórios de CSVs).")
    if not pasta_saida:
        raise ValueError("Selecione a pasta de saída da tabela consolidada.")
    if not (gerar_csv or gerar_parquet_unico or gerar_parquet_por_ano):
        raise ValueError("Selecione pelo menos um formato de saída.")

    pasta_entrada_path = Path(pasta_entrada)
    pasta_saida_path = Path(pasta_saida)
    pasta_saida_path.mkdir(parents=True, exist_ok=True)
    nome_base = normalizar_nome(nome_base or "tabela_consolidada", max_len=120)
    if not max_workers or max_workers < 1:
        max_workers = numero_padrao_de_nucleos()
    compressao = (compressao or "zstd").lower()

    quer_parquet = gerar_parquet_unico or gerar_parquet_por_ano
    tem_duckdb = _modulo_disponivel("duckdb")
    tem_pyarrow = _modulo_disponivel("pyarrow")

    # ---- Caminho preferencial: DuckDB ----
    if tem_duckdb:
        log(f"Motor: DuckDB | entrada: {'Parquet' if entrada_parquet else 'CSV'} | "
            f"threads: {max_workers} | compressão: {compressao}")
        return _consolidar_com_duckdb(
            pasta_entrada_path, pasta_saida_path, nome_base,
            gerar_csv, gerar_parquet_unico, gerar_parquet_por_ano,
            coluna_particao, compressao, max_workers,
            log, set_load, set_progress, stop_event,
            entrada_parquet=entrada_parquet,
        )

    # ---- Entrada em Parquet, sem DuckDB: junta com pyarrow num único Parquet ----
    if entrada_parquet:
        if not tem_pyarrow:
            raise ValueError("Juntar Parquets requer DuckDB ou pyarrow. Instale 'pip install duckdb' (recomendado).")
        set_load("start", "Procurando Parquets na pasta e subpastas")
        try:
            arquivos = [str(p) for p in sorted(pasta_entrada_path.rglob("*.parquet"))
                        if not (p.stem.lower().startswith(nome_base.lower())
                                or "metricas_fundido_unico" in p.stem.lower())]
        finally:
            set_load("stop", "Busca de Parquets concluída")
        if not arquivos:
            raise ValueError("Nenhum Parquet encontrado na pasta de entrada (e subpastas).")
        log(f"Parquets encontrados: {len(arquivos)}")
        saida_parquet = pasta_saida_path / f"{nome_base}.parquet"
        set_load("start", "Juntando Parquets (pyarrow)")
        try:
            _juntar_parquets_em_unico(arquivos, saida_parquet, compressao, log)
        finally:
            set_load("stop", "Parquets juntados")
        if gerar_parquet_por_ano:
            log("AVISO: Parquet por ano a partir de Parquets requer DuckDB; ignorado nesta execução.")
        set_progress(1, 1, "Consolidação concluída.")
        return {"arquivos_lidos": len(arquivos), "linhas": None,
                "saidas": {"parquet": saida_parquet}, "motor": "pyarrow"}

    # ---- Alternativa sem DuckDB (entrada CSV) ----
    log("DuckDB não está instalado. Usando o caminho alternativo (mais lento). "
        "Para máxima velocidade, instale: pip install duckdb")
    if quer_parquet and not tem_pyarrow:
        log("AVISO: 'pyarrow' também não está instalado — a saída Parquet será ignorada. "
            "Instale com: pip install pyarrow (apenas o CSV será gerado).")
        gerar_parquet_unico = False
        gerar_parquet_por_ano = False
        quer_parquet = False
    if not (gerar_csv or quer_parquet):
        raise ValueError("Sem DuckDB/pyarrow não há formato a gerar. Instale duckdb ou pyarrow.")

    # Coleta os CSVs (recursivo), ignorando saídas anteriores.
    set_load("start", "Procurando CSVs na pasta e subpastas")
    arquivos = []
    try:
        for p in pasta_entrada_path.rglob("*.csv"):
            nome = p.stem.lower()
            if (nome.startswith("metricas_fundido_geral") or nome.startswith("tabela_consolidada")
                    or nome.startswith(nome_base.lower()) or "__tmp_fundido_" in str(p)):
                continue
            arquivos.append(p)
    finally:
        set_load("stop", "Busca de CSVs concluída")
    arquivos = sorted(arquivos)
    if not arquivos:
        raise ValueError("Nenhum CSV encontrado na pasta de entrada (e subpastas).")
    total = len(arquivos)
    log(f"CSVs encontrados: {total}")

    saida_csv = pasta_saida_path / f"{nome_base}.csv"
    saida_parquet = pasta_saida_path / f"{nome_base}.parquet"
    saida_parquet_part = pasta_saida_path / f"{nome_base}_por_ano"
    total_progresso = max(1, total)
    set_progress(0, total_progresso, "Iniciando consolidação...")
    inicio = time.time()

    if gerar_parquet_por_ano and tem_pyarrow:
        saida_parquet_part.mkdir(parents=True, exist_ok=True)

    header_ref_bytes = None
    header_ref_txt = None
    f_csv = None
    writer = None
    schema_ref = None
    cols_parquet_ref = None
    divergencias_header = 0

    try:
        if gerar_csv:
            f_csv = open(saida_csv, "wb")
        for i, caminho in enumerate(arquivos, start=1):
            checar_parada(stop_event)
            try:
                header_bytes, corpo_bytes = _ler_header_e_corpo(caminho)
            except Exception as exc:
                log(f"Aviso: falha ao ler {caminho.name}: {exc}")
                continue
            if header_ref_bytes is None:
                header_ref_bytes = header_bytes
                header_ref_txt = _header_normalizado(header_bytes)
            mesmo_header = _header_normalizado(header_bytes) == header_ref_txt

            if gerar_csv:
                if i == 1:
                    f_csv.write(header_ref_bytes)
                if mesmo_header:
                    f_csv.write(corpo_bytes)
                else:
                    divergencias_header += 1
                    try:
                        df_tmp = pd.read_csv(caminho)
                        df_tmp = df_tmp.reindex(columns=header_ref_txt.split(","))
                        f_csv.write(df_tmp.to_csv(index=False, header=False).encode("utf-8"))
                    except Exception as exc:
                        log(f"Aviso: esquema divergente e falha ao realinhar {caminho.name}: {exc}")

            if gerar_parquet_unico and tem_pyarrow:
                import pyarrow as pa
                import pyarrow.parquet as pq
                try:
                    df_tmp = pd.read_csv(caminho)
                    if not df_tmp.empty:
                        if writer is None:
                            cols_parquet_ref = list(df_tmp.columns)
                            tabela = pa.Table.from_pandas(df_tmp, preserve_index=False)
                            schema_ref = tabela.schema
                            writer = pq.ParquetWriter(str(saida_parquet), schema_ref, compression=compressao)
                            writer.write_table(tabela)
                        else:
                            df_tmp = df_tmp.reindex(columns=cols_parquet_ref)
                            tabela = pa.Table.from_pandas(df_tmp, preserve_index=False)
                            if tabela.schema != schema_ref:
                                tabela = tabela.cast(schema_ref, safe=False)
                            writer.write_table(tabela)
                except Exception as exc:
                    log(f"Aviso: falha ao adicionar {caminho.name} ao Parquet: {exc}")

            if gerar_parquet_por_ano and tem_pyarrow:
                import pyarrow as pa
                import pyarrow.parquet as pq
                try:
                    df_tmp = pd.read_csv(caminho)
                    if not df_tmp.empty:
                        tabela = pa.Table.from_pandas(df_tmp, preserve_index=False)
                        if coluna_particao in df_tmp.columns:
                            pq.write_to_dataset(tabela, root_path=str(saida_parquet_part),
                                                partition_cols=[coluna_particao], compression=compressao)
                        else:
                            pq.write_to_dataset(tabela, root_path=str(saida_parquet_part), compression=compressao)
                except Exception as exc:
                    log(f"Aviso: falha ao adicionar {caminho.name} ao Parquet por ano: {exc}")

            decorrido = time.time() - inicio
            if decorrido > 0:
                taxa = i / decorrido
                eta = (total - i) / taxa if taxa > 0 else 0.0
                set_progress(i, total_progresso, f"Consolidando {i}/{total} | tempo restante ~ {formatar_duracao_horas(eta)}")
            else:
                set_progress(i, total_progresso, f"Consolidando {i}/{total}")
    finally:
        if f_csv is not None:
            f_csv.close()
        if writer is not None:
            writer.close()

    if divergencias_header:
        log(f"Aviso: {divergencias_header} arquivo(s) tinham cabeçalho diferente e foram realinhados.")
    set_progress(total_progresso, total_progresso, "Consolidação concluída.")

    saidas = {}
    if gerar_csv:
        saidas["csv"] = saida_csv
        log(f"CSV consolidado: {saida_csv}")
    if gerar_parquet_unico and tem_pyarrow and writer is not None:
        saidas["parquet"] = saida_parquet
        log(f"Parquet único: {saida_parquet}")
    if gerar_parquet_por_ano and tem_pyarrow:
        saidas["parquet_por_ano"] = saida_parquet_part
        log(f"Parquet por ano: {saida_parquet_part}")

    return {"arquivos_lidos": total, "linhas": None, "saidas": saidas, "motor": "pyarrow/python", "divergencias_header": divergencias_header}


# =============================================================================
# Interface gráfica Tkinter
# =============================================================================

class AplicacaoGradeBuffer(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("GeoLand V1")
        self.geometry("1180x860")
        self.minsize(1060, 720)

        # Aba 1
        self.shp_poligono_var = tk.StringVar()
        self.pasta_saida_grade_var = tk.StringVar()
        self.distancias_grade_var = tk.StringVar(value="100, 250, 500")
        self.unidade_grade_var = tk.StringVar(value="m")
        self.stop_grade_event = threading.Event()

        # Aba 2
        self.shp_pontos_lista: List[str] = []
        self.shp_pontos_texto_var = tk.StringVar()
        self.pasta_saida_buffer_var = tk.StringVar()
        self.diametros_buffer_var = tk.StringVar(value="100, 250, 500")
        self.unidade_buffer_var = tk.StringVar(value="m")
        self.stop_buffer_event = threading.Event()

        # Aba 3
        self.shp_buffers_recorte_var = tk.StringVar()
        self.rasters_recorte_lista: List[str] = []
        self.rasters_recorte_texto_var = tk.StringVar()
        self.pasta_saida_recorte_var = tk.StringVar()
        self.campo_id_recorte_var = tk.StringVar(value="id_pt")
        self.campo_diam_recorte_var = tk.StringVar(value="diam_m")
        self.all_touched_recorte_var = tk.BooleanVar(value=True)
        self.workers_recorte_var = tk.StringVar(value=str(numero_padrao_de_nucleos()))
        self.stop_recorte_event = threading.Event()

        # Aba 6
        self.pasta_rasters_metricas_var = tk.StringVar()
        self.pasta_saida_metricas_var = tk.StringVar()
        self.nivel_patch_var = tk.BooleanVar(value=True)
        self.nivel_class_var = tk.BooleanVar(value=True)
        self.nivel_landscape_var = tk.BooleanVar(value=True)
        self.vizinhanca_metricas_var = tk.StringVar(value="8")
        self.recursivo_metricas_var = tk.BooleanVar(value=True)
        self.reprojetar_metricas_var = tk.BooleanVar(value=False)
        self.calcular_iep_var = tk.BooleanVar(value=False)
        self.workers_metricas_var = tk.StringVar(value=str(numero_padrao_de_nucleos()))
        # Motores de cálculo (o usuário pode marcar um, o outro, ou ambos)
        self.engine_pls_var = tk.BooleanVar(value=True)
        self.engine_r_var = tk.BooleanVar(value=False)
        # cellsize do landscapemetrics (R): options(landscapemetrics_cellsize = 30)
        self.cellsize_r_var = tk.StringVar(value="30")
        self.stop_metricas_event = threading.Event()

        # Aba 9
        self.shp_buffer_join_var = tk.StringVar()
        self.csv_metricas_join_var = tk.StringVar()
        self.tipo_entrada_join_var = tk.StringVar(value="CSV")
        self.campo_id_join_var = tk.StringVar(value="id_pt")
        self.pasta_saida_join_var = tk.StringVar()
        self.stop_join_event = threading.Event()

        # Aba 4 - Reprojeção de rasters
        self.rasters_reproj_lista: List[str] = []
        self.rasters_reproj_texto_var = tk.StringVar()
        self.pasta_saida_reproj_var = tk.StringVar()
        self.epsg_reproj_var = tk.StringVar()
        self.sufixo_epsg_reproj_var = tk.BooleanVar(value=True)
        self.resampling_reproj_var = tk.StringVar(value="nearest (categórico)")
        self.workers_reproj_var = tk.StringVar(value=str(numero_padrao_de_nucleos()))
        self.stop_reproj_event = threading.Event()

        # Aba 7 - Fundido em memória (recorte + métricas, sem salvar recorte)
        self.rasters_fundido_lista: List[str] = []
        self.rasters_fundido_texto_var = tk.StringVar()
        self.shp_buffers_fundido_lista: List[str] = []
        self.shp_buffers_fundido_texto_var = tk.StringVar()
        self.pasta_saida_fundido_var = tk.StringVar()
        self.campo_id_fundido_var = tk.StringVar(value="id_pt")
        self.campo_diam_fundido_var = tk.StringVar(value="diam_m")
        self.nivel_patch_fundido_var = tk.BooleanVar(value=True)
        self.nivel_class_fundido_var = tk.BooleanVar(value=True)
        self.nivel_landscape_fundido_var = tk.BooleanVar(value=True)
        self.vizinhanca_fundido_var = tk.StringVar(value="8")
        self.all_touched_fundido_var = tk.BooleanVar(value=True)
        self.workers_fundido_var = tk.StringVar(value=str(numero_padrao_de_nucleos()))
        self.usar_processos_fundido_var = tk.BooleanVar(value=True)
        self.salvar_csv_feicao_fundido_var = tk.BooleanVar(value=False)
        self.salvar_parquet_ano_fundido_var = tk.BooleanVar(value=True)
        self.parquet_unico_final_fundido_var = tk.BooleanVar(value=True)
        self.salvar_shp_fundido_var = tk.BooleanVar(value=False)
        self.shp_consolidado_fundido_var = tk.BooleanVar(value=False)
        self.pular_reproj_fundido_var = tk.BooleanVar(value=False)
        self.retomar_fundido_var = tk.BooleanVar(value=True)
        self.stop_fundido_event = threading.Event()

        # Aba 8 - Consolidar tabelas
        self.pasta_entrada_consol_var = tk.StringVar()
        self.pasta_saida_consol_var = tk.StringVar()
        self.nome_base_consol_var = tk.StringVar(value="tabela_consolidada")
        self.consol_csv_var = tk.BooleanVar(value=True)
        self.consol_parquet_var = tk.BooleanVar(value=True)
        self.consol_parquet_part_var = tk.BooleanVar(value=True)
        self.consol_coluna_part_var = tk.StringVar(value="year")
        self.consol_compressao_var = tk.StringVar(value="zstd")
        self.consol_entrada_parquet_var = tk.BooleanVar(value=False)
        self.workers_consol_var = tk.StringVar(value=str(numero_padrao_de_nucleos()))
        self.stop_consol_event = threading.Event()

        self._criar_interface()

    # -------------------------------------------------------------------------
    # Elementos gerais
    # -------------------------------------------------------------------------

    def _criar_interface(self):
        notebook = ttk.Notebook(self)
        notebook.pack(fill="both", expand=True, padx=10, pady=10)

        # Ordem das abas = sequência habitual de processamento. Cada aba tem rolagem vertical:
        # se a janela (ex.: maximizada numa tela baixa) for menor que o conteúdo, aparece a
        # barra de rolagem em vez de o Tk espremer listas e botões até sumirem.
        self.aba_grade = self._criar_aba_rolavel(notebook, "1. Grade de pontos")
        self.aba_buffer = self._criar_aba_rolavel(notebook, "2. Buffers")
        self.aba_recorte = self._criar_aba_rolavel(notebook, "3. Recortar rasters")
        self.aba_reproj = self._criar_aba_rolavel(notebook, "4. Reprojeção de rasters")
        self.aba_reclass = self._criar_aba_rolavel(notebook, "5. Reclassificar rasters")
        self.aba_metricas = self._criar_aba_rolavel(notebook, "6. Métricas e índice")
        self.aba_fundido = self._criar_aba_rolavel(notebook, "7. Fundido (memória)")
        self.aba_consol = self._criar_aba_rolavel(notebook, "8. Consolidar tabelas")
        self.aba_join = self._criar_aba_rolavel(notebook, "9. Buffer + Métricas (SHP)")
        self.aba_graficos = self._criar_aba_rolavel(notebook, "10. Gráficos")
        self.notebook = notebook

        self._montar_aba_grade()
        self._montar_aba_buffer()
        self._montar_aba_recorte_raster()
        self._montar_aba_metricas()
        self._montar_aba_join_metricas()
        self._montar_aba_reprojecao()
        self._montar_aba_fundido()
        self._montar_aba_consolidar()
        self._montar_aba_reclass()
        # Aba 10 (gráficos): fica em um módulo próprio, geoland_graficos.py.
        from geoland_graficos import AbaGraficos
        self.graficos = AbaGraficos(self, self.aba_graficos)

    def _criar_aba_rolavel(self, notebook, texto: str) -> ttk.Frame:
        """Adiciona ao notebook uma aba com rolagem vertical e devolve o quadro interno.

        O quadro interno ocupa toda a largura e, na altura, o maior valor entre o espaço
        da aba e o que o conteúdo pede: com espaço sobrando, as linhas com peso crescem;
        com falta de espaço, nada é espremido e a barra de rolagem aparece.
        """
        externo = ttk.Frame(notebook)
        externo.rowconfigure(0, weight=1)
        externo.columnconfigure(0, weight=1)
        canvas = tk.Canvas(externo, highlightthickness=0, borderwidth=0,
                           background=ttk.Style().lookup("TFrame", "background") or None)
        barra = ttk.Scrollbar(externo, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=barra.set)
        canvas.grid(row=0, column=0, sticky="nsew")
        interno = ttk.Frame(canvas)
        janela = canvas.create_window(0, 0, window=interno, anchor="nw")
        estado = {"tamanho": None}

        def ajustar(_evento=None):
            largura, altura = canvas.winfo_width(), canvas.winfo_height()
            if largura <= 1 or altura <= 1:
                return
            pedida = interno.winfo_reqheight()
            alvo = (largura, max(pedida, altura))
            if alvo != estado["tamanho"]:
                estado["tamanho"] = alvo
                canvas.itemconfigure(janela, width=alvo[0], height=alvo[1])
                canvas.configure(scrollregion=(0, 0, alvo[0], alvo[1]))
            if pedida > altura:
                barra.grid(row=0, column=1, sticky="ns")
            else:
                barra.grid_remove()
                canvas.yview_moveto(0)

        def rolar(evento):
            # Listas, tabelas e logs rolam sozinhos; a aba só rola se tiver barra visível.
            if not barra.winfo_ismapped() or isinstance(evento.widget, (tk.Text, tk.Listbox, ttk.Treeview, tk.Canvas)) \
                    and evento.widget is not canvas:
                return
            canvas.yview_scroll(int(-evento.delta / 120) or (-1 if evento.delta > 0 else 1), "units")

        canvas.bind("<Configure>", ajustar)
        # O conteúdo pode mudar de tamanho (ex.: quadros que aparecem/somem): reavalia.
        interno.bind("<Configure>", lambda _e: interno.after_idle(ajustar))
        externo.bind("<Enter>", lambda _e: canvas.bind_all("<MouseWheel>", rolar))
        externo.bind("<Leave>", lambda _e: canvas.unbind_all("<MouseWheel>"))
        notebook.add(externo, text=texto)
        notebook.bind("<<NotebookTabChanged>>", lambda _e: self.after_idle(ajustar), add="+")
        self._ajustes_abas = getattr(self, "_ajustes_abas", []) + [ajustar]
        return interno

    def _adicionar_linha_arquivo(self, frame, label, text_var, comando, row):
        ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w", padx=6, pady=6)
        ttk.Entry(frame, textvariable=text_var, width=90).grid(row=row, column=1, sticky="ew", padx=6, pady=6)
        ttk.Button(frame, text="Selecionar", command=comando).grid(row=row, column=2, sticky="ew", padx=6, pady=6)

    def _adicionar_linha_pasta(self, frame, label, text_var, comando, row):
        ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w", padx=6, pady=6)
        ttk.Entry(frame, textvariable=text_var, width=90).grid(row=row, column=1, sticky="ew", padx=6, pady=6)
        ttk.Button(frame, text="Selecionar", command=comando).grid(row=row, column=2, sticky="ew", padx=6, pady=6)

    def _selecionar_shp(self, titulo="Selecione um shapefile"):
        return filedialog.askopenfilename(
            title=titulo,
            filetypes=[("Shapefile", "*.shp"), ("Todos os arquivos", "*.*")]
        )

    def _selecionar_multiplos_shp(self, titulo="Selecione shapefiles"):
        return filedialog.askopenfilenames(
            title=titulo,
            filetypes=[("Shapefile", "*.shp"), ("Todos os arquivos", "*.*")]
        )

    def _selecionar_multiplos_rasters(self, titulo="Selecione arquivos raster"):
        return filedialog.askopenfilenames(
            title=titulo,
            filetypes=[
                ("GeoTIFF", "*.tif *.tiff"),
                ("Rasters", "*.tif *.tiff *.img *.vrt"),
                ("Todos os arquivos", "*.*"),
            ],
        )

    def _selecionar_pasta(self, titulo="Selecione uma pasta"):
        return filedialog.askdirectory(title=titulo)

    def _selecionar_csv(self, titulo="Selecione um arquivo CSV"):
        return filedialog.askopenfilename(
            title=titulo,
            filetypes=[("CSV", "*.csv"), ("Todos os arquivos", "*.*")],
        )

    def _ler_workers(self, var: tk.StringVar) -> int:
        """Lê o número de núcleos/threads informado, com limites seguros."""
        try:
            valor = int(str(var.get()).strip())
        except (ValueError, TypeError):
            valor = numero_padrao_de_nucleos()
        cpus = os.cpu_count() or 4
        return max(1, min(valor, cpus * 2))

    def _log(self, widget_texto, mensagem: str):
        def inserir():
            widget_texto.insert("end", mensagem + "\n")
            widget_texto.see("end")
            self.update_idletasks()
        widget_texto.after(0, inserir)

    def _limpar_log(self, widget_texto):
        widget_texto.delete("1.0", "end")

    def _iniciar_barra_carregamento(self, barra, label, mensagem: str):
        def atualizar():
            barra.configure(mode="indeterminate", maximum=100)
            barra["value"] = 0
            barra.start(12)
            label.config(text=mensagem)
            self.update_idletasks()
        barra.after(0, atualizar)

    def _parar_barra_carregamento(self, barra, label, mensagem: str):
        def atualizar():
            barra.stop()
            barra.configure(mode="determinate", maximum=100)
            barra["value"] = 100
            label.config(text=mensagem)
            self.update_idletasks()
        barra.after(0, atualizar)

    def _resetar_barras(self, barras_labels):
        for barra, label, texto in barras_labels:
            def atualizar(b=barra, l=label, t=texto):
                b.stop()
                b.configure(mode="determinate", maximum=100)
                b["value"] = 0
                l.config(text=t)
            barra.after(0, atualizar)

    def _atualizar_barra_processamento(self, barra, label, valor: float, maximo: float, mensagem: str = ""):
        def atualizar():
            maximo_seguro = max(float(maximo), 1.0)
            valor_seguro = min(max(float(valor), 0.0), maximo_seguro)
            barra.configure(mode="determinate", maximum=maximo_seguro)
            barra["value"] = valor_seguro
            percentual = (valor_seguro / maximo_seguro) * 100.0
            if mensagem:
                label.config(text=f"{percentual:5.1f}% - {mensagem}")
            else:
                label.config(text=f"{percentual:5.1f}%")
            self.update_idletasks()
        barra.after(0, atualizar)

    def _executar_thread(self, botao_iniciar, botao_parar, stop_event: threading.Event, funcao):
        def alvo():
            try:
                stop_event.clear()
                botao_iniciar.after(0, lambda: botao_iniciar.config(state="disabled"))
                botao_parar.after(0, lambda: botao_parar.config(state="normal"))
                funcao()
            except ProcessamentoInterrompido as exc:
                msg = str(exc)
                botao_iniciar.after(0, lambda m=msg: messagebox.showwarning("Processamento interrompido", m))
            except Exception as exc:
                msg = str(exc)
                botao_iniciar.after(0, lambda m=msg: messagebox.showerror("Erro no processamento", m))
            finally:
                botao_parar.after(0, lambda: botao_parar.config(state="disabled"))
                botao_iniciar.after(0, lambda: botao_iniciar.config(state="normal"))
                stop_event.clear()
        threading.Thread(target=alvo, daemon=True).start()

    def _solicitar_parada(self, stop_event: threading.Event, log_widget=None):
        stop_event.set()
        if log_widget is not None:
            self._log(log_widget, "Solicitação de parada enviada. O processamento será interrompido ao finalizar a etapa atual.")

    def _adicionar_botoes_iniciar_parar(self, frame, row, texto_iniciar, comando_iniciar, comando_parar):
        container = ttk.Frame(frame)
        container.grid(row=row, column=0, columnspan=3, sticky="ew", padx=6, pady=8)
        container.columnconfigure(0, weight=3)
        container.columnconfigure(1, weight=1)
        botao_iniciar = ttk.Button(container, text=texto_iniciar, command=comando_iniciar)
        botao_iniciar.grid(row=0, column=0, sticky="ew", padx=(0, 6))
        botao_parar = ttk.Button(container, text="Parar", command=comando_parar, state="disabled")
        botao_parar.grid(row=0, column=1, sticky="ew")
        return botao_iniciar, botao_parar

    # -------------------------------------------------------------------------
    # Aba 1: grade de pontos
    # -------------------------------------------------------------------------

    def _montar_aba_grade(self):
        frame = ttk.Frame(self.aba_grade, padding=10)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)

        titulo = ttk.Label(
            frame,
            text="Criar grade uniforme de pontos dentro de shapefile de polígonos",
            font=("Segoe UI", 12, "bold"),
        )
        titulo.grid(row=0, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 10))

        self._adicionar_linha_arquivo(frame, "Shapefile de polígonos (.shp):", self.shp_poligono_var, self._acao_selecionar_poligono, row=1)
        self._adicionar_linha_pasta(frame, "Pasta de saída:", self.pasta_saida_grade_var, self._acao_selecionar_pasta_grade, row=2)

        ttk.Label(frame, text="Distâncias entre pontos:").grid(row=3, column=0, sticky="w", padx=6, pady=6)
        ttk.Entry(frame, textvariable=self.distancias_grade_var, width=40).grid(row=3, column=1, sticky="w", padx=6, pady=6)
        ttk.Label(frame, text="Unidade:").grid(row=3, column=1, sticky="e", padx=(6, 80), pady=6)
        ttk.Combobox(frame, textvariable=self.unidade_grade_var, values=["m", "km"], width=6, state="readonly").grid(row=3, column=1, sticky="e", padx=6, pady=6)

        ajuda = "Exemplos: 100, 250, 500 | 0.5, 1, 2 se a unidade for km. Cada distância gera um shapefile .zip."
        ttk.Label(frame, text=ajuda, foreground="#555555").grid(row=4, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 8))

        self.botao_grade, self.botao_parar_grade = self._adicionar_botoes_iniciar_parar(
            frame,
            row=5,
            texto_iniciar="Gerar grades e salvar ZIPs",
            comando_iniciar=lambda: self._executar_thread(self.botao_grade, self.botao_parar_grade, self.stop_grade_event, self._processar_grade),
            comando_parar=lambda: self._solicitar_parada(self.stop_grade_event, self.log_grade),
        )

        self.label_load_grade = ttk.Label(frame, text="Carregamento: aguardando arquivo")
        self.label_load_grade.grid(row=6, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_load_grade = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_load_grade.grid(row=7, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 6))

        self.label_proc_grade = ttk.Label(frame, text="Processamento: aguardando início")
        self.label_proc_grade.grid(row=8, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_proc_grade = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_proc_grade.grid(row=9, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 8))

        ttk.Label(frame, text="Log do processamento:").grid(row=10, column=0, columnspan=3, sticky="w", padx=6)
        self.log_grade = tk.Text(frame, height=18, wrap="word")
        self.log_grade.grid(row=11, column=0, columnspan=3, sticky="nsew", padx=6, pady=6)
        barra = ttk.Scrollbar(frame, orient="vertical", command=self.log_grade.yview)
        barra.grid(row=11, column=3, sticky="ns")
        self.log_grade.configure(yscrollcommand=barra.set)
        frame.rowconfigure(11, weight=1)

    def _acao_selecionar_poligono(self):
        caminho = self._selecionar_shp("Selecione o shapefile de polígonos")
        if caminho:
            self.shp_poligono_var.set(caminho)
            if not self.pasta_saida_grade_var.get():
                self.pasta_saida_grade_var.set(str(Path(caminho).parent / "saida_grade_pontos"))

    def _acao_selecionar_pasta_grade(self):
        pasta = self._selecionar_pasta("Selecione a pasta de saída das grades")
        if pasta:
            self.pasta_saida_grade_var.set(pasta)

    def _processar_grade(self):
        self.log_grade.after(0, lambda: self._limpar_log(self.log_grade))
        self._resetar_barras([
            (self.progress_load_grade, self.label_load_grade, "Carregamento: aguardando arquivo"),
            (self.progress_proc_grade, self.label_proc_grade, "Processamento: aguardando início"),
        ])
        shp_poligono = self.shp_poligono_var.get().strip()
        pasta_saida = self.pasta_saida_grade_var.get().strip()
        unidade = self.unidade_grade_var.get().strip()
        if not shp_poligono:
            raise ValueError("Selecione um shapefile de polígonos.")
        if not pasta_saida:
            raise ValueError("Selecione uma pasta de saída.")
        distancias = parse_lista_numerica(self.distancias_grade_var.get())
        self._log(self.log_grade, "Iniciando geração de grade de pontos...")
        saidas = criar_grade_pontos(
            shp_poligono=shp_poligono,
            pasta_saida=pasta_saida,
            distancias_usuario=distancias,
            unidade=unidade,
            log_callback=lambda msg: self._log(self.log_grade, msg),
            load_callback=lambda acao, msg: self._iniciar_barra_carregamento(self.progress_load_grade, self.label_load_grade, msg) if acao == "start" else self._parar_barra_carregamento(self.progress_load_grade, self.label_load_grade, msg),
            progress_callback=lambda valor, maximo, msg: self._atualizar_barra_processamento(self.progress_proc_grade, self.label_proc_grade, valor, maximo, msg),
            stop_event=self.stop_grade_event,
        )
        self._log(self.log_grade, "=" * 70)
        self._log(self.log_grade, f"Processamento finalizado. Arquivos ZIP gerados: {len(saidas)}")
        self.log_grade.after(0, lambda: messagebox.showinfo("Concluído", f"Grades geradas com sucesso.\nArquivos ZIP: {len(saidas)}"))

    # -------------------------------------------------------------------------
    # Aba 2: buffers
    # -------------------------------------------------------------------------

    def _montar_aba_buffer(self):
        frame = ttk.Frame(self.aba_buffer, padding=10)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)

        ttk.Label(frame, text="Criar buffers a partir de um ou mais shapefiles de pontos", font=("Segoe UI", 12, "bold")).grid(row=0, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 10))
        self._adicionar_linha_arquivo(frame, "Shapefile(s) de pontos (.shp):", self.shp_pontos_texto_var, self._acao_selecionar_pontos, row=1)
        self._adicionar_linha_pasta(frame, "Pasta de saída:", self.pasta_saida_buffer_var, self._acao_selecionar_pasta_buffer, row=2)

        ttk.Label(frame, text="Diâmetros dos buffers:").grid(row=3, column=0, sticky="w", padx=6, pady=6)
        ttk.Entry(frame, textvariable=self.diametros_buffer_var, width=40).grid(row=3, column=1, sticky="w", padx=6, pady=6)
        ttk.Label(frame, text="Unidade:").grid(row=3, column=1, sticky="e", padx=(6, 80), pady=6)
        ttk.Combobox(frame, textvariable=self.unidade_buffer_var, values=["m", "km"], width=6, state="readonly").grid(row=3, column=1, sticky="e", padx=6, pady=6)

        ajuda = "Atenção: os valores informados são DIÂMETROS. O script usa raio = diâmetro / 2. Cada diâmetro gera um .zip."
        ttk.Label(frame, text=ajuda, foreground="#555555").grid(row=4, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 8))

        self.botao_buffer, self.botao_parar_buffer = self._adicionar_botoes_iniciar_parar(
            frame,
            row=5,
            texto_iniciar="Gerar buffers e salvar ZIPs",
            comando_iniciar=lambda: self._executar_thread(self.botao_buffer, self.botao_parar_buffer, self.stop_buffer_event, self._processar_buffer),
            comando_parar=lambda: self._solicitar_parada(self.stop_buffer_event, self.log_buffer),
        )

        self.label_load_buffer = ttk.Label(frame, text="Carregamento: aguardando arquivo")
        self.label_load_buffer.grid(row=6, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_load_buffer = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_load_buffer.grid(row=7, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 6))
        self.label_proc_buffer = ttk.Label(frame, text="Processamento: aguardando início")
        self.label_proc_buffer.grid(row=8, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_proc_buffer = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_proc_buffer.grid(row=9, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 8))

        ttk.Label(frame, text="Log do processamento:").grid(row=10, column=0, columnspan=3, sticky="w", padx=6)
        self.log_buffer = tk.Text(frame, height=18, wrap="word")
        self.log_buffer.grid(row=11, column=0, columnspan=3, sticky="nsew", padx=6, pady=6)
        barra = ttk.Scrollbar(frame, orient="vertical", command=self.log_buffer.yview)
        barra.grid(row=11, column=3, sticky="ns")
        self.log_buffer.configure(yscrollcommand=barra.set)
        frame.rowconfigure(11, weight=1)

    def _acao_selecionar_pontos(self):
        caminhos = self._selecionar_multiplos_shp("Selecione um ou mais shapefiles de pontos")
        if caminhos:
            self.shp_pontos_lista = list(caminhos)
            self.shp_pontos_texto_var.set(caminhos[0] if len(caminhos) == 1 else f"{len(caminhos)} shapefiles selecionados")
            if not self.pasta_saida_buffer_var.get():
                self.pasta_saida_buffer_var.set(str(Path(caminhos[0]).parent / "saida_buffers"))

    def _acao_selecionar_pasta_buffer(self):
        pasta = self._selecionar_pasta("Selecione a pasta de saída dos buffers")
        if pasta:
            self.pasta_saida_buffer_var.set(pasta)

    def _processar_buffer(self):
        self.log_buffer.after(0, lambda: self._limpar_log(self.log_buffer))
        self._resetar_barras([
            (self.progress_load_buffer, self.label_load_buffer, "Carregamento: aguardando arquivo"),
            (self.progress_proc_buffer, self.label_proc_buffer, "Processamento: aguardando início"),
        ])
        pasta_saida = self.pasta_saida_buffer_var.get().strip()
        unidade = self.unidade_buffer_var.get().strip()
        if not self.shp_pontos_lista:
            raise ValueError("Selecione um ou mais shapefiles de pontos.")
        if not pasta_saida:
            raise ValueError("Selecione uma pasta de saída.")
        diametros = parse_lista_numerica(self.diametros_buffer_var.get())
        self._log(self.log_buffer, "Iniciando geração de buffers...")
        saidas = gerar_buffers_pontos(
            shp_pontos_lista=self.shp_pontos_lista,
            pasta_saida=pasta_saida,
            diametros_usuario=diametros,
            unidade=unidade,
            log_callback=lambda msg: self._log(self.log_buffer, msg),
            load_callback=lambda acao, msg: self._iniciar_barra_carregamento(self.progress_load_buffer, self.label_load_buffer, msg) if acao == "start" else self._parar_barra_carregamento(self.progress_load_buffer, self.label_load_buffer, msg),
            progress_callback=lambda valor, maximo, msg: self._atualizar_barra_processamento(self.progress_proc_buffer, self.label_proc_buffer, valor, maximo, msg),
            stop_event=self.stop_buffer_event,
        )
        self._log(self.log_buffer, "=" * 70)
        self._log(self.log_buffer, f"Processamento finalizado. Arquivos ZIP gerados: {len(saidas)}")
        self.log_buffer.after(0, lambda: messagebox.showinfo("Concluído", f"Buffers gerados com sucesso.\nArquivos ZIP: {len(saidas)}"))

    # -------------------------------------------------------------------------
    # Aba 3: recorte de rasters por buffers
    # -------------------------------------------------------------------------

    def _montar_aba_recorte_raster(self):
        frame = ttk.Frame(self.aba_recorte, padding=10)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)

        ttk.Label(frame, text="Recortar rasters usando cada feição de um shapefile de buffers", font=("Segoe UI", 12, "bold")).grid(row=0, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 10))
        self._adicionar_linha_arquivo(frame, "Shapefile de buffers (.shp):", self.shp_buffers_recorte_var, self._acao_selecionar_buffers_recorte, row=1)
        self._adicionar_linha_arquivo(frame, "Arquivos raster:", self.rasters_recorte_texto_var, self._acao_selecionar_rasters_recorte, row=2)
        self._adicionar_linha_pasta(frame, "Pasta de saída dos rasters recortados:", self.pasta_saida_recorte_var, self._acao_selecionar_pasta_recorte, row=3)

        ttk.Label(frame, text="Campo ID:").grid(row=4, column=0, sticky="w", padx=6, pady=6)
        ttk.Entry(frame, textvariable=self.campo_id_recorte_var, width=22).grid(row=4, column=1, sticky="w", padx=6, pady=6)
        ttk.Label(frame, text="Campo diâmetro:").grid(row=4, column=1, sticky="e", padx=(6, 210), pady=6)
        ttk.Entry(frame, textvariable=self.campo_diam_recorte_var, width=22).grid(row=4, column=1, sticky="e", padx=6, pady=6)

        opcoes_rec = ttk.Frame(frame)
        opcoes_rec.grid(row=5, column=0, columnspan=3, sticky="ew", padx=6, pady=6)
        ttk.Checkbutton(opcoes_rec, text="Usar all_touched=True no recorte do raster", variable=self.all_touched_recorte_var).pack(side="left")
        ttk.Label(opcoes_rec, text="Núcleos paralelos:").pack(side="left", padx=(30, 6))
        ttk.Spinbox(opcoes_rec, from_=1, to=(os.cpu_count() or 4) * 2, width=5, textvariable=self.workers_recorte_var).pack(side="left")

        ajuda = (
            "Cada raster é recortado por cada feição do buffer e salvo em GeoTIFF UTM.\n"
            "Os recortes são gravados em uma subpasta com o ANO do raster (ex.: 1985) e nomeados como "
            "clip_id_<ID>_diam_<DIAMETRO>_<ANO>.tif (ex.: clip_id_100000_diam_1000_1985.tif). "
            "Use mais núcleos para acelerar o processamento."
        )
        ttk.Label(frame, text=ajuda, foreground="#555555", justify="left").grid(row=6, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 8))

        self.botao_recorte, self.botao_parar_recorte = self._adicionar_botoes_iniciar_parar(
            frame,
            row=7,
            texto_iniciar="Recortar rasters e salvar GeoTIFFs em UTM",
            comando_iniciar=lambda: self._executar_thread(self.botao_recorte, self.botao_parar_recorte, self.stop_recorte_event, self._processar_recorte_raster),
            comando_parar=lambda: self._solicitar_parada(self.stop_recorte_event, self.log_recorte),
        )

        self.label_load_recorte = ttk.Label(frame, text="Carregamento: aguardando arquivos")
        self.label_load_recorte.grid(row=8, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_load_recorte = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_load_recorte.grid(row=9, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 6))
        self.label_proc_recorte = ttk.Label(frame, text="Processamento: aguardando início")
        self.label_proc_recorte.grid(row=10, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_proc_recorte = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_proc_recorte.grid(row=11, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 8))

        ttk.Label(frame, text="Log do processamento:").grid(row=12, column=0, columnspan=3, sticky="w", padx=6)
        self.log_recorte = tk.Text(frame, height=16, wrap="word")
        self.log_recorte.grid(row=13, column=0, columnspan=3, sticky="nsew", padx=6, pady=6)
        barra = ttk.Scrollbar(frame, orient="vertical", command=self.log_recorte.yview)
        barra.grid(row=13, column=3, sticky="ns")
        self.log_recorte.configure(yscrollcommand=barra.set)
        frame.rowconfigure(13, weight=1)

    def _acao_selecionar_buffers_recorte(self):
        caminho = self._selecionar_shp("Selecione o shapefile de buffers")
        if caminho:
            self.shp_buffers_recorte_var.set(caminho)
            if not self.pasta_saida_recorte_var.get():
                self.pasta_saida_recorte_var.set(str(Path(caminho).parent / "saida_rasters_recortados"))

    def _acao_selecionar_rasters_recorte(self):
        caminhos = self._selecionar_multiplos_rasters("Selecione um ou mais arquivos raster")
        if caminhos:
            self.rasters_recorte_lista = list(caminhos)
            self.rasters_recorte_texto_var.set(caminhos[0] if len(caminhos) == 1 else f"{len(caminhos)} rasters selecionados")
            if not self.pasta_saida_recorte_var.get():
                self.pasta_saida_recorte_var.set(str(Path(caminhos[0]).parent / "saida_rasters_recortados"))

    def _acao_selecionar_pasta_recorte(self):
        pasta = self._selecionar_pasta("Selecione a pasta de saída dos rasters recortados")
        if pasta:
            self.pasta_saida_recorte_var.set(pasta)

    def _processar_recorte_raster(self):
        self.log_recorte.after(0, lambda: self._limpar_log(self.log_recorte))
        self._resetar_barras([
            (self.progress_load_recorte, self.label_load_recorte, "Carregamento: aguardando arquivos"),
            (self.progress_proc_recorte, self.label_proc_recorte, "Processamento: aguardando início"),
        ])
        shp_buffers = self.shp_buffers_recorte_var.get().strip()
        pasta_saida = self.pasta_saida_recorte_var.get().strip()
        campo_id = self.campo_id_recorte_var.get().strip() or "id_pt"
        campo_diametro = self.campo_diam_recorte_var.get().strip() or "diam_m"
        all_touched = bool(self.all_touched_recorte_var.get())
        if not shp_buffers:
            raise ValueError("Selecione o shapefile de buffers.")
        if not self.rasters_recorte_lista:
            raise ValueError("Selecione um ou mais rasters.")
        if not pasta_saida:
            raise ValueError("Selecione a pasta de saída.")
        self._log(self.log_recorte, "Iniciando recorte dos rasters por buffers...")
        saidas = recortar_rasters_por_buffers(
            shp_buffers=shp_buffers,
            rasters_lista=self.rasters_recorte_lista,
            pasta_saida=pasta_saida,
            campo_id=campo_id,
            campo_diametro=campo_diametro,
            all_touched=all_touched,
            max_workers=self._ler_workers(self.workers_recorte_var),
            log_callback=lambda msg: self._log(self.log_recorte, msg),
            load_callback=lambda acao, msg: self._iniciar_barra_carregamento(self.progress_load_recorte, self.label_load_recorte, msg) if acao == "start" else self._parar_barra_carregamento(self.progress_load_recorte, self.label_load_recorte, msg),
            progress_callback=lambda valor, maximo, msg: self._atualizar_barra_processamento(self.progress_proc_recorte, self.label_proc_recorte, valor, maximo, msg),
            stop_event=self.stop_recorte_event,
        )
        self._log(self.log_recorte, "=" * 70)
        self._log(self.log_recorte, f"Processamento finalizado. GeoTIFFs gerados: {len(saidas)}")
        self.log_recorte.after(0, lambda: messagebox.showinfo("Concluído", f"Recortes raster gerados com sucesso.\nGeoTIFFs: {len(saidas)}"))

    # -------------------------------------------------------------------------
    # Aba 6: métricas PyLandStats
    # -------------------------------------------------------------------------

    def _montar_aba_metricas(self):
        frame = ttk.Frame(self.aba_metricas, padding=10)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)

        ttk.Label(frame, text="Calcular métricas da estrutura da paisagem com PyLandStats", font=("Segoe UI", 12, "bold")).grid(row=0, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 10))
        self._adicionar_linha_pasta(frame, "Pasta com rasters recortados/reprojetados:", self.pasta_rasters_metricas_var, self._acao_selecionar_pasta_rasters_metricas, row=1)
        self._adicionar_linha_pasta(frame, "Pasta de saída dos CSVs:", self.pasta_saida_metricas_var, self._acao_selecionar_pasta_saida_metricas, row=2)

        opcoes = ttk.LabelFrame(frame, text="Níveis escalares e vizinhança", padding=8)
        opcoes.grid(row=3, column=0, columnspan=3, sticky="ew", padx=6, pady=6)
        ttk.Checkbutton(opcoes, text="Patch", variable=self.nivel_patch_var, command=self._atualizar_estado_metricas).grid(row=0, column=0, sticky="w", padx=6)
        ttk.Checkbutton(opcoes, text="Classe", variable=self.nivel_class_var, command=self._atualizar_estado_metricas).grid(row=0, column=1, sticky="w", padx=6)
        ttk.Checkbutton(opcoes, text="Paisagem", variable=self.nivel_landscape_var, command=self._atualizar_estado_metricas).grid(row=0, column=2, sticky="w", padx=6)
        ttk.Label(opcoes, text="Vizinhança:").grid(row=0, column=3, sticky="e", padx=(30, 6))
        ttk.Combobox(opcoes, textvariable=self.vizinhanca_metricas_var, values=["4", "8"], width=6, state="readonly").grid(row=0, column=4, sticky="w", padx=6)
        ttk.Checkbutton(opcoes, text="Ler todas as subpastas (recursivo)", variable=self.recursivo_metricas_var).grid(row=0, column=5, sticky="w", padx=(30, 6))
        ttk.Label(opcoes, text="Núcleos paralelos:").grid(row=0, column=6, sticky="e", padx=(30, 6))
        ttk.Spinbox(opcoes, from_=1, to=(os.cpu_count() or 4) * 2, width=5, textvariable=self.workers_metricas_var).grid(row=0, column=7, sticky="w", padx=6)

        # Motores de cálculo (PyLandStats e/ou landscapemetrics R) + cellsize do R
        ttk.Label(opcoes, text="Motores:").grid(row=1, column=0, sticky="e", padx=6, pady=(8, 0))
        ttk.Checkbutton(opcoes, text="PyLandStats", variable=self.engine_pls_var,
                        command=self._atualizar_estado_metricas).grid(row=1, column=1, sticky="w", padx=6, pady=(8, 0))
        ttk.Checkbutton(opcoes, text="landscapemetrics (R)", variable=self.engine_r_var,
                        command=self._atualizar_estado_metricas).grid(row=1, column=2, columnspan=2, sticky="w", padx=6, pady=(8, 0))
        ttk.Label(opcoes, text="cellsize (R):").grid(row=1, column=4, sticky="e", padx=(30, 6), pady=(8, 0))
        ttk.Entry(opcoes, textvariable=self.cellsize_r_var, width=8).grid(row=1, column=5, sticky="w", padx=6, pady=(8, 0))
        ttk.Label(opcoes, text="options(landscapemetrics_cellsize)", foreground="#777777").grid(row=1, column=6, columnspan=2, sticky="w", padx=6, pady=(8, 0))

        # Projeção: botão de verificação e opção de reprojetar para metros antes do cálculo
        self.botao_verificar_proj = ttk.Button(
            opcoes, text="Verificar projeção (metros)",
            command=lambda: self._executar_thread(
                self.botao_verificar_proj, self.botao_parar_metricas, self.stop_metricas_event,
                lambda: self._com_botao_bloqueado(self.botao_metricas, self._verificar_projecao_metricas)))
        self.botao_verificar_proj.grid(row=2, column=0, columnspan=3, sticky="w", padx=6, pady=(8, 0))
        ttk.Checkbutton(
            opcoes, text="Reprojetar para metros (UTM) antes do cálculo os rasters que não estiverem em metros",
            variable=self.reprojetar_metricas_var,
        ).grid(row=2, column=3, columnspan=5, sticky="w", padx=(30, 6), pady=(8, 0))

        # Índice de Estabilidade da Paisagem: caixa separada das métricas, calculado ao final.
        caixa_iep = ttk.LabelFrame(opcoes, text="Índice de Estabilidade da Paisagem (IEP) - por buffer, classe e raio", padding=4)
        caixa_iep.grid(row=3, column=0, columnspan=8, sticky="ew", padx=6, pady=(8, 0))
        ttk.Checkbutton(
            caixa_iep, variable=self.calcular_iep_var,
            text="Calcular o IEP ao terminar as métricas de todos os anos (usa CA, PLAND, NP, LPI e ENN_MN "
                 "de classe, incluídas automaticamente)",
        ).pack(anchor="w")

        metricas_nb = ttk.Notebook(frame)
        metricas_nb.grid(row=4, column=0, columnspan=3, sticky="nsew", padx=6, pady=6)

        # Página PyLandStats (mesmo painel por nível já existente)
        pag_pls = ttk.Frame(metricas_nb, padding=6)
        for c in range(3):
            pag_pls.columnconfigure(c, weight=1)
        pag_pls.rowconfigure(0, weight=1)
        metricas_nb.add(pag_pls, text="Métricas PyLandStats")
        self.listbox_patch = self._criar_listbox_metricas(pag_pls, "Patch", PATCH_METRICS, 0)
        self.listbox_class = self._criar_listbox_metricas(pag_pls, "Classe", CLASS_METRICS, 1)
        self.listbox_landscape = self._criar_listbox_metricas(pag_pls, "Paisagem", LANDSCAPE_METRICS, 2)

        # Página landscapemetrics (R) — mesmo painel por nível, com funções lsm_*
        pag_r = ttk.Frame(metricas_nb, padding=6)
        for c in range(3):
            pag_r.columnconfigure(c, weight=1)
        pag_r.rowconfigure(0, weight=1)
        metricas_nb.add(pag_r, text="Métricas landscapemetrics (R)")
        self.listbox_r_patch = self._criar_listbox_metricas_r(pag_r, "Patch (lsm_p_*)", R_PATCH_METRICS, 0)
        self.listbox_r_class = self._criar_listbox_metricas_r(pag_r, "Classe (lsm_c_*)", R_CLASS_METRICS, 1)
        self.listbox_r_landscape = self._criar_listbox_metricas_r(pag_r, "Paisagem (lsm_l_*)", R_LANDSCAPE_METRICS, 2)

        ajuda = (
            "Escolha um ou os dois motores (PyLandStats e/ou landscapemetrics do R). A pasta é lida "
            "recursivamente (todas as subpastas). Cada raster gera um CSV por motor, salvo em subpasta do ANO "
            "(ex.: 1985); um CSV geral por motor é salvo na raiz. Colunas: source_file, source_stem, source_path, "
            "year, level, neighborhood_rule e engine. O motor R usa 'directions' (4/8) e options(landscapemetrics_cellsize)."
        )
        ttk.Label(frame, text=ajuda, foreground="#555555", wraplength=1080, justify="left").grid(row=5, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 8))

        self.botao_metricas, self.botao_parar_metricas = self._adicionar_botoes_iniciar_parar(
            frame,
            row=6,
            texto_iniciar="Calcular métricas e salvar CSVs",
            comando_iniciar=lambda: self._executar_thread(
                self.botao_metricas, self.botao_parar_metricas, self.stop_metricas_event,
                lambda: self._com_botao_bloqueado(self.botao_verificar_proj, self._processar_metricas)),
            comando_parar=lambda: self._solicitar_parada(self.stop_metricas_event, self.log_metricas),
        )

        self.label_load_metricas = ttk.Label(frame, text="Carregamento: aguardando pasta de rasters")
        self.label_load_metricas.grid(row=7, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_load_metricas = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_load_metricas.grid(row=8, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 6))

        self.label_proc_metricas = ttk.Label(frame, text="Processamento: aguardando início")
        self.label_proc_metricas.grid(row=9, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_proc_metricas = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_proc_metricas.grid(row=10, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 8))

        ttk.Label(frame, text="Log do processamento:").grid(row=11, column=0, columnspan=3, sticky="w", padx=6)
        self.log_metricas = tk.Text(frame, height=10, wrap="word")
        self.log_metricas.grid(row=12, column=0, columnspan=3, sticky="nsew", padx=6, pady=6)
        barra = ttk.Scrollbar(frame, orient="vertical", command=self.log_metricas.yview)
        barra.grid(row=12, column=3, sticky="ns")
        self.log_metricas.configure(yscrollcommand=barra.set)
        frame.rowconfigure(4, weight=3)  # as listas de métricas cedem mais espaço que o log
        frame.rowconfigure(12, weight=1)

        self._selecionar_todas_metricas(self.listbox_patch)
        self._selecionar_todas_metricas(self.listbox_class)
        self._selecionar_todas_metricas(self.listbox_landscape)
        self._selecionar_todas_metricas(self.listbox_r_patch)
        self._selecionar_todas_metricas(self.listbox_r_class)
        self._selecionar_todas_metricas(self.listbox_r_landscape)
        self._atualizar_estado_metricas()

    def _criar_listbox_metricas(self, parent, titulo: str, metricas: Sequence[str], coluna: int):
        container = ttk.Frame(parent)
        container.grid(row=0, column=coluna, sticky="nsew", padx=5, pady=5)
        container.columnconfigure(0, weight=1)
        ttk.Label(container, text=titulo, font=("Segoe UI", 10, "bold")).grid(row=0, column=0, columnspan=2, sticky="w")
        listbox = tk.Listbox(container, selectmode="extended", height=12, exportselection=False)
        listbox.grid(row=1, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(container, orient="vertical", command=listbox.yview)
        scroll.grid(row=1, column=1, sticky="ns")
        listbox.configure(yscrollcommand=scroll.set)
        for metrica in metricas:
            listbox.insert("end", metrica)
        botoes = ttk.Frame(container)
        botoes.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(4, 0))
        ttk.Button(botoes, text="Todas", command=lambda lb=listbox: self._selecionar_todas_metricas(lb)).pack(side="left", fill="x", expand=True, padx=(0, 2))
        ttk.Button(botoes, text="Limpar", command=lambda lb=listbox: self._limpar_metricas(lb)).pack(side="left", fill="x", expand=True, padx=(2, 0))
        container.rowconfigure(1, weight=1)
        return listbox

    def _criar_listbox_metricas_r(self, parent, titulo: str, pares: Sequence, coluna: int):
        """Listbox de métricas do R. Exibe 'CÓDIGO — lsm_x_...' e guarda a função
        lsm_* correspondente em listbox._rfuncs (índice -> função)."""
        container = ttk.Frame(parent)
        container.grid(row=0, column=coluna, sticky="nsew", padx=5, pady=5)
        container.columnconfigure(0, weight=1)
        ttk.Label(container, text=titulo, font=("Segoe UI", 10, "bold")).grid(row=0, column=0, columnspan=2, sticky="w")
        listbox = tk.Listbox(container, selectmode="extended", height=12, exportselection=False)
        listbox.grid(row=1, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(container, orient="vertical", command=listbox.yview)
        scroll.grid(row=1, column=1, sticky="ns")
        listbox.configure(yscrollcommand=scroll.set)
        funcs = []
        for codigo, func in pares:
            listbox.insert("end", "%s  —  %s" % (codigo, func))
            funcs.append(func)
        listbox._rfuncs = funcs
        botoes = ttk.Frame(container)
        botoes.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(4, 0))
        ttk.Button(botoes, text="Todas", command=lambda lb=listbox: self._selecionar_todas_metricas(lb)).pack(side="left", fill="x", expand=True, padx=(0, 2))
        ttk.Button(botoes, text="Limpar", command=lambda lb=listbox: self._limpar_metricas(lb)).pack(side="left", fill="x", expand=True, padx=(2, 0))
        container.rowconfigure(1, weight=1)
        return listbox

    def _obter_metricas_r(self, listbox) -> List[str]:
        """Funções lsm_* selecionadas numa listbox do R (vazio se desabilitada)."""
        if str(listbox.cget("state")) == "disabled":
            return []
        funcs = getattr(listbox, "_rfuncs", [])
        return [funcs[i] for i in listbox.curselection() if 0 <= i < len(funcs)]

    def _selecionar_todas_metricas(self, listbox):
        listbox.select_set(0, "end")

    def _limpar_metricas(self, listbox):
        listbox.selection_clear(0, "end")

    def _atualizar_estado_metricas(self):
        pls_on = bool(self.engine_pls_var.get())
        r_on = bool(self.engine_r_var.get())
        niveis = {
            "patch": bool(self.nivel_patch_var.get()),
            "class": bool(self.nivel_class_var.get()),
            "landscape": bool(self.nivel_landscape_var.get()),
        }
        pls_map = {"patch": getattr(self, "listbox_patch", None),
                   "class": getattr(self, "listbox_class", None),
                   "landscape": getattr(self, "listbox_landscape", None)}
        r_map = {"patch": getattr(self, "listbox_r_patch", None),
                 "class": getattr(self, "listbox_r_class", None),
                 "landscape": getattr(self, "listbox_r_landscape", None)}
        for nivel, lb in pls_map.items():
            if lb is not None:
                lb.configure(state="normal" if (pls_on and niveis[nivel]) else "disabled")
        for nivel, lb in r_map.items():
            if lb is not None:
                lb.configure(state="normal" if (r_on and niveis[nivel]) else "disabled")

    def _obter_metricas_listbox(self, listbox) -> List[str]:
        if str(listbox.cget("state")) == "disabled":
            return []
        return [listbox.get(i) for i in listbox.curselection()]

    def _acao_selecionar_pasta_rasters_metricas(self):
        pasta = self._selecionar_pasta("Selecione a pasta com os rasters recortados/reprojetados")
        if pasta:
            self.pasta_rasters_metricas_var.set(pasta)
            if not self.pasta_saida_metricas_var.get():
                self.pasta_saida_metricas_var.set(str(Path(pasta) / "csv_metricas_paisagem"))

    def _acao_selecionar_pasta_saida_metricas(self):
        pasta = self._selecionar_pasta("Selecione a pasta onde os CSVs serão salvos")
        if pasta:
            self.pasta_saida_metricas_var.set(pasta)

    def _com_botao_bloqueado(self, botao, funcao):
        """Executa `funcao` mantendo `botao` desabilitado (um processamento por vez na aba)."""
        botao.after(0, lambda: botao.config(state="disabled"))
        try:
            funcao()
        finally:
            botao.after(0, lambda: botao.config(state="normal"))

    def _verificar_projecao_metricas(self):
        """Botão "Verificar projeção (metros)": informa no log cada raster que não está em metros."""
        self.log_metricas.after(0, lambda: self._limpar_log(self.log_metricas))
        self._resetar_barras([
            (self.progress_load_metricas, self.label_load_metricas, "Carregamento: aguardando pasta de rasters"),
            (self.progress_proc_metricas, self.label_proc_metricas, "Processamento: aguardando início"),
        ])
        pasta_rasters = self.pasta_rasters_metricas_var.get().strip()
        if not pasta_rasters:
            raise ValueError("Selecione a pasta com os rasters recortados/reprojetados.")
        log_cb = lambda msg: self._log(self.log_metricas, msg)
        rasters = listar_rasters_pasta(pasta_rasters, recursivo=bool(self.recursivo_metricas_var.get()))
        self._parar_barra_carregamento(self.progress_load_metricas, self.label_load_metricas,
                                       f"Rasters encontrados: {len(rasters)}")
        if not rasters:
            raise ValueError("Nenhum raster .tif/.tiff foi encontrado na pasta selecionada.")
        log_cb(f"Verificando a projeção de {len(rasters)} raster(s) em: {pasta_rasters}")
        resultado = verificar_projecao_rasters(
            rasters,
            progress_callback=lambda valor, maximo, msg: self._atualizar_barra_processamento(
                self.progress_proc_metricas, self.label_proc_metricas, valor, maximo, msg),
            stop_event=self.stop_metricas_event,
        )
        resumo = registrar_verificacao_projecao(resultado, pasta_rasters, log_cb)
        texto = resumo.replace(" | ", "\n")
        if resultado["nao_metros"]:
            dica = 'Marque "Reprojetar para metros (UTM) antes do cálculo" para reprojetá-los automaticamente.'
            log_cb(dica)
            self.log_metricas.after(0, lambda: messagebox.showwarning("Projeção dos rasters", texto + "\n\n" + dica))
        else:
            self.log_metricas.after(0, lambda: messagebox.showinfo("Projeção dos rasters", texto))

    def _processar_metricas(self):
        self.log_metricas.after(0, lambda: self._limpar_log(self.log_metricas))
        self._resetar_barras([
            (self.progress_load_metricas, self.label_load_metricas, "Carregamento: aguardando pasta de rasters"),
            (self.progress_proc_metricas, self.label_proc_metricas, "Processamento: aguardando início"),
        ])

        pasta_rasters = self.pasta_rasters_metricas_var.get().strip()
        pasta_saida = self.pasta_saida_metricas_var.get().strip()
        if not pasta_rasters:
            raise ValueError("Selecione a pasta com os rasters recortados/reprojetados.")
        if not pasta_saida:
            raise ValueError("Selecione a pasta de saída dos CSVs.")

        niveis = {
            "patch": bool(self.nivel_patch_var.get()),
            "class": bool(self.nivel_class_var.get()),
            "landscape": bool(self.nivel_landscape_var.get()),
        }
        if not any(niveis.values()):
            raise ValueError("Selecione ao menos um nível escalar (patch, classe ou paisagem).")

        usar_pls = bool(self.engine_pls_var.get())
        usar_r = bool(self.engine_r_var.get())
        if not (usar_pls or usar_r):
            raise ValueError("Selecione ao menos um motor: PyLandStats e/ou landscapemetrics (R).")

        viz = self.vizinhanca_metricas_var.get()
        recursivo = bool(self.recursivo_metricas_var.get())
        workers = self._ler_workers(self.workers_metricas_var)

        # Callbacks de log/carregamento/progresso reaproveitados pelos dois motores.
        log_cb = lambda msg: self._log(self.log_metricas, msg)
        load_cb = lambda acao, msg: (
            self._iniciar_barra_carregamento(self.progress_load_metricas, self.label_load_metricas, msg)
            if acao == "start" else
            self._parar_barra_carregamento(self.progress_load_metricas, self.label_load_metricas, msg))
        prog_cb = lambda valor, maximo, msg: self._atualizar_barra_processamento(
            self.progress_proc_metricas, self.label_proc_metricas, valor, maximo, msg)

        self._log(self.log_metricas, "Iniciando cálculo das métricas da paisagem...")
        self._log(self.log_metricas, f"Níveis: {', '.join([k for k, v in niveis.items() if v])}")
        self._log(self.log_metricas, f"Vizinhança: {viz} | recursivo: {recursivo} | núcleos: {workers}")
        motores_txt = " + ".join([m for m, on in [("PyLandStats", usar_pls), ("landscapemetrics (R)", usar_r)] if on])
        self._log(self.log_metricas, f"Motores selecionados: {motores_txt}")

        # Índice de estabilidade (IEP): precisa do nível de classe e de cinco métricas de classe.
        import geoland_indice
        calcular_iep = bool(self.calcular_iep_var.get())
        if calcular_iep:
            niveis["class"] = True
            self._log(self.log_metricas, "IEP marcado: nível Classe e as métricas CA, PLAND, NP, LPI e ENN_MN "
                                         "entram no cálculo; o índice sai ao final, com todos os anos.")

        # Opcional: verifica a projeção e reprojeta para metros os rasters que não estiverem em metros.
        rasters_lista = None
        if bool(self.reprojetar_metricas_var.get()):
            self._log(self.log_metricas, "=" * 70)
            rasters_lista = preparar_rasters_em_metros(
                pasta_rasters=pasta_rasters,
                pasta_saida=pasta_saida,
                recursivo=recursivo,
                max_workers=workers,
                log_callback=log_cb,
                load_callback=load_cb,
                progress_callback=prog_cb,
                stop_event=self.stop_metricas_event,
            )

        resumo = []

        # ---------------- Motor PyLandStats ----------------
        if usar_pls:
            metricas_por_nivel = {
                "patch": self._obter_metricas_listbox(self.listbox_patch),
                "class": self._obter_metricas_listbox(self.listbox_class),
                "landscape": self._obter_metricas_listbox(self.listbox_landscape),
            }
            if calcular_iep:
                metricas_por_nivel["class"] = list(dict.fromkeys(
                    list(metricas_por_nivel["class"]) + list(geoland_indice.METRICAS_PYLANDSTATS)))
            self._log(self.log_metricas, "=" * 70)
            self._log(self.log_metricas, "Motor: PyLandStats")
            res_pls = calcular_metricas_paisagem_pylandstats(
                pasta_rasters=pasta_rasters,
                pasta_saida=pasta_saida,
                niveis=niveis,
                metricas_por_nivel=metricas_por_nivel,
                vizinhanca=viz,
                recursivo=recursivo,
                max_workers=workers,
                engine_label="pylandstats",
                log_callback=log_cb,
                load_callback=load_cb,
                progress_callback=prog_cb,
                stop_event=self.stop_metricas_event,
                rasters_lista=rasters_lista,
            )
            resumo.append(("PyLandStats", res_pls))
            self._log(self.log_metricas,
                      f"[PyLandStats] rasters={res_pls['rasters_processados']} | "
                      f"CSVs individuais={len(res_pls['csvs_individuais'])} | CSV geral={res_pls['csv_geral']}")

        # ------------- Motor landscapemetrics (R) -------------
        if usar_r:
            try:
                cellsize = float(str(self.cellsize_r_var.get()).replace(",", "."))
            except Exception:
                cellsize = 30.0
            metricas_func_por_nivel = {
                "patch": self._obter_metricas_r(self.listbox_r_patch),
                "class": self._obter_metricas_r(self.listbox_r_class),
                "landscape": self._obter_metricas_r(self.listbox_r_landscape),
            }
            if calcular_iep:
                metricas_func_por_nivel["class"] = list(dict.fromkeys(
                    list(metricas_func_por_nivel["class"]) + list(geoland_indice.FUNCOES_LSM)))
            self._log(self.log_metricas, "=" * 70)
            self._log(self.log_metricas, f"Motor: landscapemetrics (R)  |  cellsize={cellsize}")
            res_r = calcular_metricas_paisagem_r(
                pasta_rasters=pasta_rasters,
                pasta_saida=pasta_saida,
                niveis=niveis,
                metricas_func_por_nivel=metricas_func_por_nivel,
                vizinhanca=viz,
                cellsize=cellsize,
                recursivo=recursivo,
                max_workers=workers,
                log_callback=log_cb,
                load_callback=load_cb,
                progress_callback=prog_cb,
                stop_event=self.stop_metricas_event,
                rasters_lista=rasters_lista,
            )
            resumo.append(("landscapemetrics (R)", res_r))
            self._log(self.log_metricas,
                      f"[R] rasters={res_r['rasters_processados']} | "
                      f"CSVs individuais={len(res_r['csvs_individuais'])} | CSV geral={res_r['csv_geral']}")

        # Índice de Estabilidade da Paisagem: só agora, com as métricas de toda a série temporal.
        if calcular_iep:
            self._log(self.log_metricas, "=" * 70)
            self._log(self.log_metricas, "Índice de Estabilidade da Paisagem (IEP)")
            for nome, r in resumo:
                if not r.get("csv_geral"):
                    continue
                motor = "pylandstats" if nome == "PyLandStats" else "landscapemetrics"
                destino = str(Path(pasta_saida) / f"indice_estabilidade_IEP_{motor}_{datetime.now():%Y%m%d_%H%M%S}.csv")
                try:
                    info = geoland_indice.salvar_iep_do_csv_geral(str(r["csv_geral"]), destino)
                    self._log(self.log_metricas,
                              f"[IEP - {nome}] {info['linhas']} combinações de buffer x classe x raio | até "
                              f"{info['anos']} anos | IEP médio = {info['iep_medio']:.3f} | {destino}")
                    if info["anos"] < 2:
                        self._log(self.log_metricas, "ATENÇÃO: há um único ano; sem série temporal o CV é 0 e o IEP vale 1.")
                except Exception as exc:
                    self._log(self.log_metricas, f"[IEP - {nome}] não calculado: {exc}")

        self._log(self.log_metricas, "=" * 70)
        linhas = "\n".join(
            f"- {nome}: {r['rasters_processados']} rasters, {len(r['csvs_individuais'])} CSVs, geral: {r['csv_geral']}"
            for nome, r in resumo)
        self._log(self.log_metricas, "Concluído.\n" + linhas)
        self.log_metricas.after(0, lambda: messagebox.showinfo("Concluído", "Métricas calculadas com sucesso.\n\n" + linhas))

    # -------------------------------------------------------------------------
    # Aba 9: buffer + métricas -> shapefiles (.shp)
    # -------------------------------------------------------------------------

    def _montar_aba_join_metricas(self):
        frame = ttk.Frame(self.aba_join, padding=10)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)

        ttk.Label(
            frame,
            text="Gerar shapefiles (.shp) do buffer enriquecidos com as métricas da paisagem (CSV ou Parquet)",
            font=("Segoe UI", 12, "bold"),
        ).grid(row=0, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 10))

        self._adicionar_linha_arquivo(frame, "Shapefile de buffer (.shp):", self.shp_buffer_join_var, self._acao_selecionar_buffer_join, row=1)
        self._adicionar_linha_arquivo(frame, "Métricas (CSV ou Parquet):", self.csv_metricas_join_var, self._acao_selecionar_csv_join, row=2)
        self._adicionar_linha_pasta(frame, "Pasta de saída:", self.pasta_saida_join_var, self._acao_selecionar_pasta_join, row=3)

        linha = ttk.Frame(frame)
        linha.grid(row=4, column=0, columnspan=3, sticky="w", padx=6, pady=6)
        ttk.Label(linha, text="Campo ID:").pack(side="left")
        ttk.Entry(linha, textvariable=self.campo_id_join_var, width=18).pack(side="left", padx=(4, 20))
        ttk.Label(linha, text="Tipo de entrada:").pack(side="left")
        ttk.Combobox(linha, textvariable=self.tipo_entrada_join_var, values=["CSV", "Parquet"], width=10, state="readonly").pack(side="left", padx=(4, 0))

        ajuda = (
            "Cada feição do buffer é pareada pelo 'Campo ID' (id_pt) com o ID das métricas. Para CSV usa a coluna "
            "'id_pt'/'id' se existir, senão o ID extraído de 'source_file' (id_<N>); para Parquet usa 'id_pt'/'id'. "
            "Como há várias classes e anos por ID, a geometria é replicada por linha de métrica (formato longo). "
            "Em Parquet, pode-se apontar um arquivo .parquet ou uma pasta com vários (lidos juntos).\n"
            "Saídas (.shp): buffer_metricas_<ano>.shp por ano, buffer_metricas_TODOS_ANOS.shp consolidado, "
            "e buffer_metricas_mapa_campos.csv mapeando nomes curtos (<=10) -> nomes originais. O CRS do buffer é preservado."
        )
        ttk.Label(frame, text=ajuda, foreground="#555555", justify="left", wraplength=980).grid(row=5, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 8))

        self.botao_join, self.botao_parar_join = self._adicionar_botoes_iniciar_parar(
            frame,
            row=6,
            texto_iniciar="Processar e gerar shapefiles (.shp)",
            comando_iniciar=lambda: self._executar_thread(self.botao_join, self.botao_parar_join, self.stop_join_event, self._processar_join_metricas),
            comando_parar=lambda: self._solicitar_parada(self.stop_join_event, self.log_join),
        )

        self.label_load_join = ttk.Label(frame, text="Carregamento: aguardando arquivos")
        self.label_load_join.grid(row=7, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_load_join = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_load_join.grid(row=8, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 6))
        self.label_proc_join = ttk.Label(frame, text="Processamento: aguardando início")
        self.label_proc_join.grid(row=9, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_proc_join = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_proc_join.grid(row=10, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 8))

        ttk.Label(frame, text="Log do processamento:").grid(row=11, column=0, columnspan=3, sticky="w", padx=6)
        self.log_join = tk.Text(frame, height=16, wrap="word")
        self.log_join.grid(row=12, column=0, columnspan=3, sticky="nsew", padx=6, pady=6)
        barra = ttk.Scrollbar(frame, orient="vertical", command=self.log_join.yview)
        barra.grid(row=12, column=3, sticky="ns")
        self.log_join.configure(yscrollcommand=barra.set)
        frame.rowconfigure(12, weight=1)

    def _acao_selecionar_buffer_join(self):
        caminho = self._selecionar_shp("Selecione o shapefile de buffer")
        if caminho:
            self.shp_buffer_join_var.set(caminho)
            if not self.pasta_saida_join_var.get():
                self.pasta_saida_join_var.set(str(Path(caminho).parent / "saida_buffer_metricas"))

    def _acao_selecionar_csv_join(self):
        if self.tipo_entrada_join_var.get().lower().startswith("parquet"):
            from tkinter import filedialog
            caminho = filedialog.askopenfilename(
                title="Selecione o Parquet de métricas",
                filetypes=[("Parquet", "*.parquet"), ("Todos os arquivos", "*.*")],
            )
        else:
            caminho = self._selecionar_csv("Selecione o CSV de métricas")
        if caminho:
            self.csv_metricas_join_var.set(caminho)
            if not self.pasta_saida_join_var.get():
                self.pasta_saida_join_var.set(str(Path(caminho).parent / "saida_buffer_metricas"))

    def _acao_selecionar_pasta_join(self):
        pasta = self._selecionar_pasta("Selecione a pasta de saída dos shapefiles")
        if pasta:
            self.pasta_saida_join_var.set(pasta)

    def _processar_join_metricas(self):
        self.log_join.after(0, lambda: self._limpar_log(self.log_join))
        self._resetar_barras([
            (self.progress_load_join, self.label_load_join, "Carregamento: aguardando arquivos"),
            (self.progress_proc_join, self.label_proc_join, "Processamento: aguardando início"),
        ])

        shp_buffer = self.shp_buffer_join_var.get().strip()
        csv_metricas = self.csv_metricas_join_var.get().strip()
        pasta_saida = self.pasta_saida_join_var.get().strip()
        campo_id = self.campo_id_join_var.get().strip() or "id_pt"
        tipo_entrada = "parquet" if self.tipo_entrada_join_var.get().lower().startswith("parquet") else "csv"

        if not shp_buffer:
            raise ValueError("Selecione o shapefile de buffer.")
        if not csv_metricas:
            raise ValueError("Selecione o arquivo de métricas (CSV ou Parquet).")
        if not pasta_saida:
            raise ValueError("Selecione a pasta de saída.")

        self._log(self.log_join, "Iniciando junção do buffer com as métricas da paisagem...")
        self._log(self.log_join, f"Tipo de entrada: {tipo_entrada.upper()} | Campo ID de pareamento: '{campo_id}'")

        resultado = gerar_shapefiles_buffer_metricas(
            shp_buffer=shp_buffer,
            csv_metricas=csv_metricas,
            pasta_saida=pasta_saida,
            campo_id=campo_id,
            tipo_entrada=tipo_entrada,
            log_callback=lambda msg: self._log(self.log_join, msg),
            load_callback=lambda acao, msg: self._iniciar_barra_carregamento(self.progress_load_join, self.label_load_join, msg) if acao == "start" else self._parar_barra_carregamento(self.progress_load_join, self.label_load_join, msg),
            progress_callback=lambda valor, maximo, msg: self._atualizar_barra_processamento(self.progress_proc_join, self.label_proc_join, valor, maximo, msg),
            stop_event=self.stop_join_event,
        )

        shapefiles = resultado["shapefiles"]
        anos = resultado["anos"]
        self._log(self.log_join, "=" * 70)
        self._log(self.log_join, f"Shapefiles gerados: {len(shapefiles)} (anos: {', '.join(anos)} + consolidado)")
        self._log(self.log_join, f"Feições no consolidado: {resultado['feicoes_total']:,}".replace(",", "."))
        self._log(self.log_join, f"Campos renomeados (>10 caracteres): {len(resultado['campos_renomeados'])}")
        self._log(self.log_join, f"IDs do CSV sem buffer: {len(resultado['csv_sem_buffer'])} | IDs do buffer sem CSV: {len(resultado['buffer_sem_csv'])}")
        self._log(self.log_join, f"Mapa de campos: {resultado['mapa_campos']}")
        self.log_join.after(0, lambda: messagebox.showinfo(
            "Concluído",
            "Shapefiles gerados com sucesso.\n"
            f"Anuais: {len(anos)} | Consolidado: 1\n"
            f"Feições no consolidado: {resultado['feicoes_total']}\n"
            f"Mapa de campos: {Path(resultado['mapa_campos']).name}",
        ))

    # -------------------------------------------------------------------------
    # Aba 4: reprojeção de rasters em lote
    # -------------------------------------------------------------------------

    def _montar_aba_reprojecao(self):
        frame = ttk.Frame(self.aba_reproj, padding=10)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)

        ttk.Label(
            frame,
            text="Reprojetar rasters em lote (sugestão automática de EPSG, editável)",
            font=("Segoe UI", 12, "bold"),
        ).grid(row=0, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 10))

        self._adicionar_linha_pasta(frame, "Pasta dos rasters de entrada:", self.rasters_reproj_texto_var, self._acao_selecionar_rasters_reproj, row=1)
        self._adicionar_linha_pasta(frame, "Pasta de saída:", self.pasta_saida_reproj_var, self._acao_selecionar_pasta_reproj, row=2)

        opcoes = ttk.LabelFrame(frame, text="EPSG de destino e opções", padding=8)
        opcoes.grid(row=3, column=0, columnspan=3, sticky="ew", padx=6, pady=6)
        ttk.Label(opcoes, text="EPSG (editável):").grid(row=0, column=0, sticky="e", padx=6)
        ttk.Entry(opcoes, textvariable=self.epsg_reproj_var, width=12).grid(row=0, column=1, sticky="w", padx=6)
        ttk.Button(opcoes, text="Sugerir EPSG", command=self._acao_sugerir_epsg_reproj).grid(row=0, column=2, sticky="w", padx=6)
        ttk.Label(opcoes, text="Reamostragem:").grid(row=0, column=3, sticky="e", padx=(30, 6))
        ttk.Combobox(opcoes, textvariable=self.resampling_reproj_var, values=list(RESAMPLING_REPROJECAO.keys()), width=18, state="readonly").grid(row=0, column=4, sticky="w", padx=6)
        ttk.Label(opcoes, text="Núcleos:").grid(row=0, column=5, sticky="e", padx=(30, 6))
        ttk.Spinbox(opcoes, from_=1, to=(os.cpu_count() or 4) * 2, width=5, textvariable=self.workers_reproj_var).grid(row=0, column=6, sticky="w", padx=6)
        ttk.Checkbutton(opcoes, text="Acrescentar sufixo do EPSG ao nome (ex.: _31984)", variable=self.sufixo_epsg_reproj_var).grid(row=1, column=0, columnspan=5, sticky="w", padx=6, pady=(6, 0))

        ajuda = (
            "Selecione a PASTA dos rasters: todos os .tif/.tiff da pasta e das subpastas são carregados "
            "automaticamente. Ao selecionar, o EPSG é sugerido automaticamente "
            "(SIRGAS 2000 / UTM pela zona do centróide; ex.: 31978–31985 no Sul), mas o campo é editável. "
            "Cada raster reprojetado é gravado na pasta de saída, na mesma estrutura de subpastas da entrada, "
            "preservando o nome (com sufixo opcional do EPSG). Erros em um arquivo não interrompem o lote."
        )
        ttk.Label(frame, text=ajuda, foreground="#555555", justify="left", wraplength=1040).grid(row=4, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 8))

        self.botao_reproj, self.botao_parar_reproj = self._adicionar_botoes_iniciar_parar(
            frame,
            row=5,
            texto_iniciar="Reprojetar rasters em lote",
            comando_iniciar=lambda: self._executar_thread(self.botao_reproj, self.botao_parar_reproj, self.stop_reproj_event, self._processar_reprojecao),
            comando_parar=lambda: self._solicitar_parada(self.stop_reproj_event, self.log_reproj),
        )

        self.label_load_reproj = ttk.Label(frame, text="Carregamento: aguardando arquivos")
        self.label_load_reproj.grid(row=6, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_load_reproj = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_load_reproj.grid(row=7, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 6))
        self.label_proc_reproj = ttk.Label(frame, text="Processamento: aguardando início")
        self.label_proc_reproj.grid(row=8, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_proc_reproj = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_proc_reproj.grid(row=9, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 8))

        ttk.Label(frame, text="Log do processamento:").grid(row=10, column=0, columnspan=3, sticky="w", padx=6)
        self.log_reproj = tk.Text(frame, height=16, wrap="word")
        self.log_reproj.grid(row=11, column=0, columnspan=3, sticky="nsew", padx=6, pady=6)
        barra = ttk.Scrollbar(frame, orient="vertical", command=self.log_reproj.yview)
        barra.grid(row=11, column=3, sticky="ns")
        self.log_reproj.configure(yscrollcommand=barra.set)
        frame.rowconfigure(11, weight=1)

    def _acao_selecionar_rasters_reproj(self):
        pasta = self._selecionar_pasta("Selecione a pasta com os rasters para reprojetar (as subpastas são lidas automaticamente)")
        if not pasta:
            return
        self.rasters_reproj_texto_var.set(pasta)
        if not self.pasta_saida_reproj_var.get():
            # Pasta irmã da entrada: fora dela, para a saída não ser relida como entrada.
            p = Path(pasta)
            self.pasta_saida_reproj_var.set(str(p.parent / f"{p.name}_reprojetado"))
        self._carregar_rasters_reproj(avisar=True)
        # Sugestão automática de EPSG a partir do primeiro raster.
        if self.rasters_reproj_lista and not self.epsg_reproj_var.get().strip():
            self._acao_sugerir_epsg_reproj()

    def _carregar_rasters_reproj(self, avisar: bool = False) -> List[str]:
        """Lê a pasta de entrada da aba 4 e todas as subpastas e carrega a lista de rasters."""
        pasta = self.rasters_reproj_texto_var.get().strip()
        self.rasters_reproj_lista = []
        if not pasta:
            return []
        try:
            rasters = listar_rasters_pasta(pasta, recursivo=True)
        except ValueError:
            rasters = []
        # Não relê como entrada o que já estiver dentro da pasta de saída.
        saida = self.pasta_saida_reproj_var.get().strip()
        if saida:
            saida_abs = os.path.normcase(os.path.abspath(saida)) + os.sep
            rasters = [r for r in rasters if not os.path.normcase(os.path.abspath(r)).startswith(saida_abs)]
        self.rasters_reproj_lista = [str(r) for r in rasters]
        n_pastas = len({str(r.parent) for r in rasters})
        self._log(self.log_reproj, f"Pasta de entrada: {pasta}")
        self._log(self.log_reproj, f"Rasters carregados (pasta e subpastas): {len(rasters)} em {n_pastas} pasta(s)")
        self._parar_barra_carregamento(self.progress_load_reproj, self.label_load_reproj,
                                       f"Rasters carregados: {len(rasters)} em {n_pastas} pasta(s)")
        if avisar and not rasters:
            messagebox.showwarning("Rasters de entrada", "Nenhum raster .tif/.tiff foi encontrado na pasta nem nas subpastas.")
        return self.rasters_reproj_lista

    def _acao_selecionar_pasta_reproj(self):
        pasta = self._selecionar_pasta("Selecione a pasta de saída dos rasters reprojetados")
        if pasta:
            self.pasta_saida_reproj_var.set(pasta)

    def _acao_sugerir_epsg_reproj(self):
        if not self.rasters_reproj_lista:
            self._carregar_rasters_reproj()  # a pasta pode ter sido digitada no campo
        if not self.rasters_reproj_lista:
            messagebox.showinfo("Sugerir EPSG", "Selecione primeiro a pasta com os rasters.")
            return
        epsg = sugerir_epsg_de_raster(self.rasters_reproj_lista[0])
        if epsg:
            self.epsg_reproj_var.set(str(epsg))
            self._log(self.log_reproj, f"EPSG sugerido a partir de {Path(self.rasters_reproj_lista[0]).name}: {epsg}")
        else:
            messagebox.showwarning("Sugerir EPSG", "Não foi possível sugerir o EPSG (raster sem CRS?). Informe manualmente.")

    def _processar_reprojecao(self):
        self.log_reproj.after(0, lambda: self._limpar_log(self.log_reproj))
        self._resetar_barras([
            (self.progress_load_reproj, self.label_load_reproj, "Carregamento: aguardando arquivos"),
            (self.progress_proc_reproj, self.label_proc_reproj, "Processamento: aguardando início"),
        ])
        pasta_entrada = self.rasters_reproj_texto_var.get().strip()
        if not pasta_entrada:
            raise ValueError("Selecione a pasta com os rasters de entrada.")
        pasta_saida = self.pasta_saida_reproj_var.get().strip()
        if not pasta_saida:
            raise ValueError("Selecione a pasta de saída.")
        if os.path.normcase(os.path.abspath(pasta_saida)) == os.path.normcase(os.path.abspath(pasta_entrada)):
            raise ValueError("A pasta de saída deve ser diferente da pasta de entrada.")
        # Relê a pasta e as subpastas na hora de processar (o conteúdo pode ter mudado).
        if not self._carregar_rasters_reproj():
            raise ValueError("Nenhum raster .tif/.tiff foi encontrado na pasta de entrada nem nas subpastas.")
        epsg_txt = self.epsg_reproj_var.get().strip()
        if not epsg_txt:
            raise ValueError("Informe (ou sugira) o EPSG de destino.")
        try:
            epsg_destino = int(epsg_txt)
        except ValueError:
            raise ValueError(f"EPSG inválido: {epsg_txt}")

        resampling = RESAMPLING_REPROJECAO.get(self.resampling_reproj_var.get(), Resampling.nearest)
        self._log(self.log_reproj, "Iniciando reprojeção em lote...")
        resultado = reprojetar_rasters_lote(
            rasters_lista=self.rasters_reproj_lista,
            pasta_saida=pasta_saida,
            epsg_destino=epsg_destino,
            sufixo_epsg=bool(self.sufixo_epsg_reproj_var.get()),
            resampling=resampling,
            max_workers=self._ler_workers(self.workers_reproj_var),
            log_callback=lambda msg: self._log(self.log_reproj, msg),
            load_callback=lambda acao, msg: self._iniciar_barra_carregamento(self.progress_load_reproj, self.label_load_reproj, msg) if acao == "start" else self._parar_barra_carregamento(self.progress_load_reproj, self.label_load_reproj, msg),
            progress_callback=lambda valor, maximo, msg: self._atualizar_barra_processamento(self.progress_proc_reproj, self.label_proc_reproj, valor, maximo, msg),
            stop_event=self.stop_reproj_event,
            pasta_base=pasta_entrada,
        )
        self._log(self.log_reproj, "=" * 70)
        self._log(self.log_reproj, f"Reprojetados: {len(resultado['saidas'])}/{resultado['total']} | erros: {len(resultado['erros'])}")
        self.log_reproj.after(0, lambda: messagebox.showinfo(
            "Concluído",
            f"Rasters reprojetados: {len(resultado['saidas'])}/{resultado['total']}\nErros: {len(resultado['erros'])}",
        ))

    # -------------------------------------------------------------------------
    # Aba 5: reclassificar rasters por regras (portada do GeoRasterTools)
    # -------------------------------------------------------------------------

    def _montar_aba_reclass(self):
        self.pasta_entrada_reclass_var = tk.StringVar()
        self.pasta_saida_reclass_var = tk.StringVar()
        self.recursivo_reclass_var = tk.BooleanVar(value=True)
        self.prefixo_reclass_var = tk.StringVar(value="Reclass_")
        self.dtype_reclass_var = tk.StringVar(value=rotulo_dtype_reclass("INT1U"))
        self.outros_reclass_var = tk.StringVar(value="NA (sem valor)")
        self.intervalo_reclass_var = tk.StringVar(value="Ambos os lados (NA)")
        self.incluir_minimo_reclass_var = tk.BooleanVar(value=True)
        self.sobrescrever_reclass_var = tk.BooleanVar(value=True)
        self.workers_reclass_var = tk.StringVar(value=str(numero_padrao_de_nucleos()))
        self.regra_ativa_reclass_var = tk.StringVar(value=NOME_REGRA_PADRAO)
        self.nova_regra_reclass_var = tk.StringVar()
        self.stop_reclass_event = threading.Event()
        # Estado: arquivos escaneados ({"Arquivo", "Caminho", "Regra"}) e conjuntos de regras por nome.
        self.reclass_arquivos: List[Dict[str, str]] = []
        self.reclass_regras: Dict[str, list] = {NOME_REGRA_PADRAO: regras_reclass_padrao()}
        self.reclass_regra_ativa = NOME_REGRA_PADRAO

        frame = ttk.Frame(self.aba_reclass, padding=10)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)

        ttk.Label(
            frame, text="Reclassificar rasters por regras (Val_Min, Val_Max, Novo_Val)",
            font=("Segoe UI", 12, "bold"),
        ).grid(row=0, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 6))
        self._adicionar_linha_pasta(frame, "Pasta dos rasters de entrada:", self.pasta_entrada_reclass_var, self._acao_selecionar_pasta_entrada_reclass, row=1)
        self._adicionar_linha_pasta(frame, "Pasta de saída:", self.pasta_saida_reclass_var, self._acao_selecionar_pasta_saida_reclass, row=2)

        corpo = ttk.Frame(frame)
        corpo.grid(row=3, column=0, columnspan=3, sticky="nsew", padx=6, pady=4)
        corpo.columnconfigure(1, weight=3)
        corpo.columnconfigure(2, weight=2)
        corpo.rowconfigure(0, weight=1)
        corpo.rowconfigure(1, weight=1)

        # Parâmetros de saída (coluna da esquerda)
        params = ttk.LabelFrame(corpo, text="Parâmetros de saída", padding=8)
        params.grid(row=0, column=0, rowspan=2, sticky="nsw", padx=(0, 6))
        ttk.Label(params, text="Tipo de dado de saída:").pack(anchor="w")
        ttk.Combobox(params, textvariable=self.dtype_reclass_var, width=26, state="readonly",
                     values=[rotulo_dtype_reclass(c) for c in RECLASS_DTYPES]).pack(fill="x", pady=2)
        ttk.Label(params, text="Pixels não mapeados:").pack(anchor="w")
        ttk.Combobox(params, textvariable=self.outros_reclass_var, state="readonly",
                     values=list(RECLASS_NAO_MAPEADOS)).pack(fill="x", pady=2)
        ttk.Label(params, text="Intervalo fechado em:").pack(anchor="w")
        ttk.Combobox(params, textvariable=self.intervalo_reclass_var, state="readonly",
                     values=list(RECLASS_INTERVALOS)).pack(fill="x", pady=2)
        ttk.Checkbutton(params, text="Incluir valor mínimo (include.lowest)",
                        variable=self.incluir_minimo_reclass_var).pack(anchor="w", pady=2)
        ttk.Checkbutton(params, text="Sobrescrever arquivos existentes",
                        variable=self.sobrescrever_reclass_var).pack(anchor="w", pady=2)
        ttk.Checkbutton(params, text="Ler subpastas (replica a estrutura na saída)",
                        variable=self.recursivo_reclass_var).pack(anchor="w", pady=2)
        linha_prefixo = ttk.Frame(params)
        linha_prefixo.pack(fill="x", pady=2)
        ttk.Label(linha_prefixo, text="Prefixo do nome de saída:").pack(side="left")
        ttk.Entry(linha_prefixo, textvariable=self.prefixo_reclass_var, width=10).pack(side="left", padx=4)
        linha_nucleos = ttk.Frame(params)
        linha_nucleos.pack(fill="x", pady=2)
        ttk.Label(linha_nucleos, text="Núcleos:").pack(side="left")
        ttk.Spinbox(linha_nucleos, from_=1, to=(os.cpu_count() or 4) * 2, width=5,
                    textvariable=self.workers_reclass_var).pack(side="left", padx=4)
        ttk.Button(params, text="Escanear arquivos", command=self._acao_escanear_reclass).pack(fill="x", pady=(8, 0))

        # Arquivos escaneados e a regra de cada um
        arquivos_box = ttk.LabelFrame(corpo, text="Arquivos escaneados (selecione linhas para atribuir regras)", padding=4)
        arquivos_box.grid(row=0, column=1, sticky="nsew", padx=(0, 6))
        arquivos_box.columnconfigure(0, weight=1)
        arquivos_box.rowconfigure(0, weight=1)
        self.tree_arquivos_reclass = ttk.Treeview(arquivos_box, columns=("Arquivo", "Regra"), show="headings",
                                                  height=6, selectmode="extended")
        self.tree_arquivos_reclass.heading("Arquivo", text="Arquivo")
        self.tree_arquivos_reclass.heading("Regra", text="Regra")
        self.tree_arquivos_reclass.column("Arquivo", width=320)
        self.tree_arquivos_reclass.column("Regra", width=90, anchor="center")
        self.tree_arquivos_reclass.grid(row=0, column=0, sticky="nsew")
        barra_arquivos = ttk.Scrollbar(arquivos_box, orient="vertical", command=self.tree_arquivos_reclass.yview)
        barra_arquivos.grid(row=0, column=1, sticky="ns")
        self.tree_arquivos_reclass.configure(yscrollcommand=barra_arquivos.set)

        # Conjuntos de regras
        regras_box = ttk.LabelFrame(corpo, text="Conjuntos de regras", padding=6)
        regras_box.grid(row=0, column=2, sticky="nsew")
        ttk.Label(regras_box, text="Regra ativa no editor:").pack(anchor="w")
        self.combo_regra_reclass = ttk.Combobox(regras_box, textvariable=self.regra_ativa_reclass_var,
                                                values=[NOME_REGRA_PADRAO], state="readonly")
        self.combo_regra_reclass.pack(fill="x", pady=2)
        self.combo_regra_reclass.bind("<<ComboboxSelected>>", self._acao_selecionar_regra_reclass)
        linha_nova = ttk.Frame(regras_box)
        linha_nova.pack(fill="x", pady=2)
        ttk.Entry(linha_nova, textvariable=self.nova_regra_reclass_var).pack(side="left", fill="x", expand=True)
        ttk.Button(linha_nova, text="Criar", width=7, command=self._acao_criar_regra_reclass).pack(side="left", padx=(4, 0))
        ttk.Button(regras_box, text="Atribuir regra ativa aos selecionados",
                   command=self._acao_atribuir_selecionados_reclass).pack(fill="x", pady=2)
        ttk.Button(regras_box, text="Atribuir regra ativa a todos os arquivos",
                   command=self._acao_atribuir_todos_reclass).pack(fill="x", pady=2)
        linha_txt = ttk.Frame(regras_box)
        linha_txt.pack(fill="x", pady=2)
        ttk.Button(linha_txt, text="Carregar .txt", command=self._acao_carregar_txt_reclass).pack(
            side="left", fill="x", expand=True, padx=(0, 2))
        ttk.Button(linha_txt, text="Salvar .txt", command=self._acao_salvar_txt_reclass).pack(
            side="left", fill="x", expand=True, padx=(2, 0))

        # Editor da regra ativa
        editor_box = ttk.LabelFrame(corpo, text="Editor de regras", padding=4)
        editor_box.grid(row=1, column=1, columnspan=2, sticky="nsew", pady=(6, 0))
        editor_box.columnconfigure(0, weight=1)
        editor_box.rowconfigure(1, weight=1)
        cabecalho = ttk.Frame(editor_box)
        cabecalho.grid(row=0, column=0, columnspan=2, sticky="ew", pady=(0, 2))
        self.label_editor_reclass = ttk.Label(cabecalho, text=f"EDITOR - {NOME_REGRA_PADRAO}", font=("Segoe UI", 9, "bold"))
        self.label_editor_reclass.pack(side="left")
        ttk.Button(cabecalho, text="+ Linha", command=self._acao_adicionar_linha_reclass).pack(side="left", padx=4)
        ttk.Button(cabecalho, text="- Linha selecionada", command=self._acao_remover_linha_reclass).pack(side="left", padx=4)
        ttk.Button(cabecalho, text="Salvar regra", command=self._acao_salvar_regra_reclass).pack(side="left", padx=4)
        ttk.Label(cabecalho, text="Duplo clique em uma célula para editar.", foreground="#555555").pack(side="left", padx=8)
        self.tree_regras_reclass = ttk.Treeview(editor_box, columns=COLUNAS_REGRA_RECLASS, show="headings", height=6)
        for coluna in COLUNAS_REGRA_RECLASS:
            self.tree_regras_reclass.heading(coluna, text=coluna)
            self.tree_regras_reclass.column(coluna, width=100, anchor="center")
        self.tree_regras_reclass.grid(row=1, column=0, sticky="nsew")
        barra_regras = ttk.Scrollbar(editor_box, orient="vertical", command=self.tree_regras_reclass.yview)
        barra_regras.grid(row=1, column=1, sticky="ns")
        self.tree_regras_reclass.configure(yscrollcommand=barra_regras.set)
        self.tree_regras_reclass.bind("<Double-1>", self._acao_editar_celula_reclass)

        self.botao_reclass, self.botao_parar_reclass = self._adicionar_botoes_iniciar_parar(
            frame,
            row=4,
            texto_iniciar="Executar reclassificação",
            comando_iniciar=lambda: self._executar_thread(self.botao_reclass, self.botao_parar_reclass, self.stop_reclass_event, self._processar_reclass),
            comando_parar=lambda: self._solicitar_parada(self.stop_reclass_event, self.log_reclass),
        )
        self.label_proc_reclass = ttk.Label(frame, text="Processamento: aguardando início")
        self.label_proc_reclass.grid(row=5, column=0, columnspan=3, sticky="w", padx=6, pady=(2, 0))
        self.progress_proc_reclass = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_proc_reclass.grid(row=6, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 6))

        ttk.Label(frame, text="Log do processamento:").grid(row=7, column=0, columnspan=3, sticky="w", padx=6)
        self.log_reclass = tk.Text(frame, height=7, wrap="word")
        self.log_reclass.grid(row=8, column=0, columnspan=3, sticky="nsew", padx=6, pady=6)
        barra = ttk.Scrollbar(frame, orient="vertical", command=self.log_reclass.yview)
        barra.grid(row=8, column=3, sticky="ns")
        self.log_reclass.configure(yscrollcommand=barra.set)
        frame.rowconfigure(3, weight=3)
        frame.rowconfigure(8, weight=1)

        self._atualizar_editor_reclass()

    def _atualizar_editor_reclass(self):
        self.label_editor_reclass.config(text=f"EDITOR - {self.reclass_regra_ativa}")
        tree = self.tree_regras_reclass
        tree.delete(*tree.get_children())
        for regra in self.reclass_regras[self.reclass_regra_ativa]:
            tree.insert("", "end", values=[f"{v:g}" for v in regra])

    def _atualizar_arquivos_reclass(self):
        tree = self.tree_arquivos_reclass
        tree.delete(*tree.get_children())
        for arquivo in self.reclass_arquivos:
            tree.insert("", "end", values=(arquivo["Arquivo"], arquivo["Regra"]))

    def _linhas_editor_reclass(self) -> List[tuple]:
        tree = self.tree_regras_reclass
        return [tuple(float(v) for v in tree.item(item, "values")) for item in tree.get_children()]

    def _ativar_regra_reclass(self, nome: str):
        self.reclass_regra_ativa = nome
        self.regra_ativa_reclass_var.set(nome)
        self.combo_regra_reclass.config(values=list(self.reclass_regras))
        self._atualizar_editor_reclass()

    def _acao_selecionar_pasta_entrada_reclass(self):
        pasta = self._selecionar_pasta("Selecione a pasta com os rasters para reclassificar")
        if not pasta:
            return
        self.pasta_entrada_reclass_var.set(pasta)
        if not self.pasta_saida_reclass_var.get():
            # Pasta irmã da entrada: fora dela, para a saída não ser relida como entrada.
            p = Path(pasta)
            self.pasta_saida_reclass_var.set(str(p.parent / f"{p.name}_reclassificado"))
        self._acao_escanear_reclass()

    def _acao_selecionar_pasta_saida_reclass(self):
        pasta = self._selecionar_pasta("Selecione a pasta de saída dos rasters reclassificados")
        if pasta:
            self.pasta_saida_reclass_var.set(pasta)

    def _acao_escanear_reclass(self):
        pasta = self.pasta_entrada_reclass_var.get().strip()
        if not pasta:
            messagebox.showwarning("Aviso", "Selecione a pasta de entrada.")
            return
        try:
            rasters = listar_rasters_pasta(pasta, recursivo=bool(self.recursivo_reclass_var.get()))
        except ValueError:
            rasters = []
        # Não relê como entrada o que já estiver dentro da pasta de saída.
        saida = self.pasta_saida_reclass_var.get().strip()
        if saida:
            saida_abs = os.path.normcase(os.path.abspath(saida)) + os.sep
            rasters = [r for r in rasters if not os.path.normcase(os.path.abspath(r)).startswith(saida_abs)]
        # Ao reescanear, cada arquivo mantém a regra que já tinha. Arquivos novos recebem a regra
        # que está no editor (antes recebiam sempre a "Padrao": quem carregava o .txt ANTES de
        # escolher a pasta tinha os arquivos reclassificados pela Padrao, que gera o valor 1).
        self._guardar_editor_reclass()
        regra_nova = self.reclass_regra_ativa
        anteriores = {a["Caminho"]: a["Regra"] for a in self.reclass_arquivos if a["Regra"] in self.reclass_regras}
        self.reclass_arquivos = [
            {"Arquivo": os.path.relpath(r, pasta), "Caminho": str(r), "Regra": anteriores.get(str(r), regra_nova)}
            for r in rasters
        ]
        self._atualizar_arquivos_reclass()
        if not rasters:
            messagebox.showwarning("Aviso", "Nenhum .tif/.tiff encontrado na pasta de entrada.")
            return
        n_pastas = len({str(r.parent) for r in rasters})
        self._log(self.log_reclass, f"=== ESCANEADO: {len(rasters)} arquivo(s) em {n_pastas} pasta(s) ===")
        self._log(self.log_reclass, f"Arquivos novos recebem a regra do editor ('{regra_nova}'). "
                                    "Selecione arquivos e atribua outras regras se necessário.")

    def _guardar_editor_reclass(self):
        """Grava o conteúdo do editor na regra ativa (antes de trocar de regra ou executar)."""
        try:
            self.reclass_regras[self.reclass_regra_ativa] = self._linhas_editor_reclass()
        except (ValueError, KeyError):
            pass

    def _acao_selecionar_regra_reclass(self, _evento=None):
        nome = self.regra_ativa_reclass_var.get()
        if nome in self.reclass_regras:
            self._guardar_editor_reclass()  # não perde edições ao trocar de regra
            self.reclass_regra_ativa = nome
            self._atualizar_editor_reclass()

    def _acao_criar_regra_reclass(self):
        nome = self.nova_regra_reclass_var.get().strip()
        if not nome:
            messagebox.showwarning("Aviso", "Digite um nome para a nova regra.")
            return
        if nome not in self.reclass_regras:
            self.reclass_regras[nome] = regras_reclass_padrao()
        self._ativar_regra_reclass(nome)
        self.nova_regra_reclass_var.set("")
        messagebox.showinfo("Regra criada", f"Regra criada: {nome}")

    def _acao_atribuir_selecionados_reclass(self):
        tree = self.tree_arquivos_reclass
        selecao = tree.selection()
        if not selecao:
            messagebox.showwarning("Aviso", "Selecione arquivos na tabela antes de atribuir.")
            return
        for item in selecao:
            self.reclass_arquivos[tree.index(item)]["Regra"] = self.reclass_regra_ativa
        self._atualizar_arquivos_reclass()
        messagebox.showinfo("Atribuído", f"Regra '{self.reclass_regra_ativa}' atribuída a {len(selecao)} arquivo(s).")

    def _acao_atribuir_todos_reclass(self):
        if not self.reclass_arquivos:
            return
        for arquivo in self.reclass_arquivos:
            arquivo["Regra"] = self.reclass_regra_ativa
        self._atualizar_arquivos_reclass()
        messagebox.showinfo("Atribuído", f"Regra '{self.reclass_regra_ativa}' atribuída a todos os "
                                         f"{len(self.reclass_arquivos)} arquivo(s).")

    def _acao_editar_celula_reclass(self, evento):
        """Duplo clique em uma célula do editor de regras: abre um campo de edição sobre ela."""
        tree = self.tree_regras_reclass
        if tree.identify("region", evento.x, evento.y) != "cell":
            return
        linha, coluna = tree.identify_row(evento.y), tree.identify_column(evento.x)
        caixa = tree.bbox(linha, coluna) if linha else None
        if not caixa:
            return
        nome_coluna = tree["columns"][int(coluna.replace("#", "")) - 1]
        editor = ttk.Entry(tree)
        editor.insert(0, tree.set(linha, nome_coluna))
        editor.select_range(0, "end")
        editor.place(x=caixa[0], y=caixa[1], width=caixa[2], height=caixa[3])
        editor.focus_set()

        def confirmar(_evento=None):
            if not editor.winfo_exists():
                return
            texto = editor.get().strip().replace(",", ".")
            editor.destroy()
            try:
                valor = float(texto)
            except ValueError:
                messagebox.showwarning("Aviso", "Valor inválido. Use apenas números.")
                return
            tree.set(linha, nome_coluna, f"{valor:g}")

        editor.bind("<Return>", confirmar)
        editor.bind("<FocusOut>", confirmar)
        editor.bind("<Escape>", lambda _evento: editor.destroy())

    def _acao_adicionar_linha_reclass(self):
        self.tree_regras_reclass.insert("", "end", values=(0, 0, 0))

    def _acao_remover_linha_reclass(self):
        tree = self.tree_regras_reclass
        selecao = tree.selection()
        if not selecao:
            messagebox.showwarning("Aviso", "Selecione uma linha para remover.")
            return
        if len(tree.get_children()) - len(selecao) < 1:
            messagebox.showwarning("Aviso", "A regra deve ter ao menos uma linha.")
            return
        tree.delete(*selecao)

    def _acao_salvar_regra_reclass(self):
        linhas = self._linhas_editor_reclass()
        self.reclass_regras[self.reclass_regra_ativa] = linhas
        messagebox.showinfo("Salvo", f"Regra '{self.reclass_regra_ativa}' salva com {len(linhas)} linha(s).")

    def _acao_carregar_txt_reclass(self):
        caminho = filedialog.askopenfilename(title="Selecionar arquivo de regras",
                                             filetypes=[("Texto", "*.txt"), ("Todos os arquivos", "*.*")])
        if not caminho:
            return
        try:
            regras = ler_regras_reclass_txt(caminho)
        except OSError as exc:
            messagebox.showerror("Erro ao ler", str(exc))
            return
        if not regras:
            messagebox.showwarning("Aviso", "Nenhuma regra válida encontrada no arquivo.")
            return
        nome = Path(caminho).stem
        self._guardar_editor_reclass()
        self.reclass_regras[nome] = regras
        self._ativar_regra_reclass(nome)
        # Correção: antes, a regra carregada ficava só no editor e os arquivos continuavam
        # na regra "Padrao" (agregação MapBiomas), que era a aplicada na execução.
        # Agora a regra carregada passa a valer para os arquivos que estavam na "Padrao".
        n = 0
        for arquivo in self.reclass_arquivos:
            if arquivo["Regra"] == NOME_REGRA_PADRAO:
                arquivo["Regra"] = nome
                n += 1
        self._atualizar_arquivos_reclass()
        self._log(self.log_reclass, f"Regra '{nome}' carregada de {Path(caminho).name} ({len(regras)} linhas) "
                                    f"e atribuída a {n} arquivo(s) que estavam na regra '{NOME_REGRA_PADRAO}'.")
        messagebox.showinfo("Carregado", f"{len(regras)} regra(s) carregada(s) de '{Path(caminho).name}' como '{nome}'.\n"
                                         f"A regra foi atribuída a {n} arquivo(s) que estavam na regra '{NOME_REGRA_PADRAO}'.")

    def _acao_salvar_txt_reclass(self):
        caminho = filedialog.asksaveasfilename(title="Salvar regras em arquivo", defaultextension=".txt",
                                               filetypes=[("Texto", "*.txt")])
        if not caminho:
            return
        try:
            salvar_regras_reclass_txt(self._linhas_editor_reclass(), caminho, nome_regra=self.reclass_regra_ativa)
            messagebox.showinfo("Salvo", f"Regras salvas em: {Path(caminho).name}")
        except Exception as exc:
            messagebox.showerror("Erro ao salvar", str(exc))

    def _processar_reclass(self):
        self.log_reclass.after(0, lambda: self._limpar_log(self.log_reclass))
        self._resetar_barras([
            (self.progress_proc_reclass, self.label_proc_reclass, "Processamento: aguardando início"),
        ])
        if not self.reclass_arquivos:
            raise ValueError("Escaneie os arquivos primeiro (selecione a pasta de entrada).")
        pasta_entrada = self.pasta_entrada_reclass_var.get().strip()
        pasta_saida = self.pasta_saida_reclass_var.get().strip()
        if not pasta_saida:
            raise ValueError("Selecione a pasta de saída.")
        prefixo = self.prefixo_reclass_var.get().strip()
        if re.search(r'[\\/:*?"<>|]', prefixo):
            raise ValueError('O prefixo não pode conter \\ / : * ? " < > |')
        if not prefixo and os.path.normcase(os.path.abspath(pasta_saida)) == os.path.normcase(os.path.abspath(pasta_entrada)):
            raise ValueError("Com o prefixo vazio, a pasta de saída deve ser diferente da pasta de entrada.")

        # O que está no editor vale para a regra ativa, mesmo sem clicar em "Salvar regra".
        self.reclass_regras[self.reclass_regra_ativa] = self._linhas_editor_reclass()
        dtype = self.dtype_reclass_var.get().split("-")[0].strip().upper()
        outros = RECLASS_NAO_MAPEADOS.get(self.outros_reclass_var.get(), "NA")
        arquivos = [(a["Caminho"], a["Regra"]) for a in self.reclass_arquivos]

        self._log(self.log_reclass, "=== RECLASSIFICAÇÃO INICIADA ===")
        # Mostra exatamente qual regra vai para quantos arquivos (evita aplicar a regra errada sem perceber).
        contagem: Dict[str, int] = {}
        for _caminho, nome_regra in arquivos:
            contagem[nome_regra] = contagem.get(nome_regra, 0) + 1
        for nome_regra, n in contagem.items():
            regras = self.reclass_regras.get(nome_regra) or []
            resumo = "; ".join(f"{a:g}-{b:g}→{c:g}" if a != b else f"{a:g}→{c:g}" for a, b, c in regras[:12])
            self._log(self.log_reclass, f"Regra '{nome_regra}' em {n} arquivo(s): {resumo}{' ...' if len(regras) > 12 else ''}")
        if self.reclass_regra_ativa not in contagem:
            self._log(self.log_reclass, f"AVISO: a regra do editor ('{self.reclass_regra_ativa}') não está atribuída a "
                                        "nenhum arquivo. Use 'Atribuir regra ativa...' se ela é a que deve ser aplicada.")
        self._log(self.log_reclass, f"Arquivos: {len(arquivos)} | Tipo de dado: {self.dtype_reclass_var.get()} | "
                                    f"Não mapeados: {outros} | Intervalo: {self.intervalo_reclass_var.get()}")
        resultado = reclassificar_rasters_lote(
            arquivos=arquivos,
            regras_por_nome=dict(self.reclass_regras),
            pasta_saida=pasta_saida,
            dtype=dtype,
            outros=outros,
            direita=RECLASS_INTERVALOS.get(self.intervalo_reclass_var.get()),
            incluir_minimo=bool(self.incluir_minimo_reclass_var.get()),
            sobrescrever=bool(self.sobrescrever_reclass_var.get()),
            prefixo=prefixo,
            pasta_base=pasta_entrada,
            max_workers=self._ler_workers(self.workers_reclass_var),
            log_callback=lambda msg: self._log(self.log_reclass, msg),
            progress_callback=lambda valor, maximo, msg: self._atualizar_barra_processamento(self.progress_proc_reclass, self.label_proc_reclass, valor, maximo, msg),
            stop_event=self.stop_reclass_event,
        )
        n_ok, n_erro = len(resultado["saidas"]), len(resultado["erros"])
        self._log(self.log_reclass, f"=== CONCLUÍDO | OK: {n_ok} | ERRO: {n_erro} ===")
        presos = resultado.get("auxiliares_presos") or []
        aviso = (f"\n\nATENÇÃO: em {len(presos)} arquivo(s) a tabela de atributos antiga (.vat.dbf/.aux.xml) "
                 "não pôde ser apagada (raster aberto no ArcGIS?). Veja o log." if presos else "")
        self.log_reclass.after(0, lambda: messagebox.showinfo(
            "Concluído", f"Reclassificação - OK: {n_ok} | Erros: {n_erro}{aviso}"))

    # -------------------------------------------------------------------------
    # Aba 7: fundido em memória (recorte + métricas, sem salvar recorte)
    # -------------------------------------------------------------------------

    def _montar_aba_fundido(self):
        frame = ttk.Frame(self.aba_fundido, padding=10)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)

        ttk.Label(
            frame,
            text="Recorte + métricas EM MEMÓRIA (sem salvar rasters recortados)",
            font=("Segoe UI", 12, "bold"),
        ).grid(row=0, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 10))

        self._adicionar_linha_arquivo(frame, "Rasters a processar:", self.rasters_fundido_texto_var, self._acao_selecionar_rasters_fundido, row=1)
        self._adicionar_linha_arquivo(frame, "Shapefiles de buffers:", self.shp_buffers_fundido_texto_var, self._acao_selecionar_buffers_fundido, row=2)
        self._adicionar_linha_pasta(frame, "Pasta de saída dos CSVs:", self.pasta_saida_fundido_var, self._acao_selecionar_pasta_fundido, row=3)

        campos = ttk.Frame(frame)
        campos.grid(row=4, column=0, columnspan=3, sticky="ew", padx=6, pady=6)
        ttk.Label(campos, text="Campo ID:").pack(side="left")
        ttk.Entry(campos, textvariable=self.campo_id_fundido_var, width=18).pack(side="left", padx=(4, 20))
        ttk.Label(campos, text="Campo diâmetro:").pack(side="left")
        ttk.Entry(campos, textvariable=self.campo_diam_fundido_var, width=18).pack(side="left", padx=(4, 0))

        opcoes = ttk.LabelFrame(frame, text="Níveis, vizinhança e opções", padding=8)
        opcoes.grid(row=5, column=0, columnspan=3, sticky="ew", padx=6, pady=6)
        ttk.Checkbutton(opcoes, text="Patch", variable=self.nivel_patch_fundido_var, command=self._atualizar_estado_metricas_fundido).grid(row=0, column=0, sticky="w", padx=6)
        ttk.Checkbutton(opcoes, text="Classe", variable=self.nivel_class_fundido_var, command=self._atualizar_estado_metricas_fundido).grid(row=0, column=1, sticky="w", padx=6)
        ttk.Checkbutton(opcoes, text="Paisagem", variable=self.nivel_landscape_fundido_var, command=self._atualizar_estado_metricas_fundido).grid(row=0, column=2, sticky="w", padx=6)
        ttk.Label(opcoes, text="Vizinhança:").grid(row=0, column=3, sticky="e", padx=(20, 6))
        ttk.Combobox(opcoes, textvariable=self.vizinhanca_fundido_var, values=["4", "8"], width=6, state="readonly").grid(row=0, column=4, sticky="w", padx=6)
        ttk.Checkbutton(opcoes, text="all_touched", variable=self.all_touched_fundido_var).grid(row=0, column=5, sticky="w", padx=(20, 6))
        ttk.Label(opcoes, text="Núcleos:").grid(row=0, column=6, sticky="e", padx=(20, 6))
        ttk.Spinbox(opcoes, from_=1, to=(os.cpu_count() or 4) * 2, width=5, textvariable=self.workers_fundido_var).grid(row=0, column=7, sticky="w", padx=6)
        ttk.Checkbutton(opcoes, text="Usar processos (ProcessPoolExecutor)", variable=self.usar_processos_fundido_var).grid(row=1, column=0, columnspan=3, sticky="w", padx=6, pady=(6, 0))
        ttk.Checkbutton(opcoes, text="Salvar CSV por feição (id/diam/ano)", variable=self.salvar_csv_feicao_fundido_var).grid(row=1, column=3, columnspan=2, sticky="w", padx=6, pady=(6, 0))
        ttk.Checkbutton(opcoes, text="Pular reprojeção se origem já em metros", variable=self.pular_reproj_fundido_var).grid(row=1, column=5, columnspan=3, sticky="w", padx=6, pady=(6, 0))
        ttk.Checkbutton(opcoes, text="Salvar Parquet por ano (1 por pasta de ano)", variable=self.salvar_parquet_ano_fundido_var).grid(row=2, column=0, columnspan=3, sticky="w", padx=6, pady=(6, 0))
        ttk.Checkbutton(opcoes, text="Gerar Parquet único ao final", variable=self.parquet_unico_final_fundido_var).grid(row=2, column=3, columnspan=2, sticky="w", padx=6, pady=(6, 0))
        ttk.Checkbutton(opcoes, text="Retomar de onde parou (pular rasters já processados)", variable=self.retomar_fundido_var).grid(row=2, column=5, columnspan=3, sticky="w", padx=6, pady=(6, 0))
        ttk.Checkbutton(opcoes, text="Salvar shapefile com métricas (como a aba 9)", variable=self.salvar_shp_fundido_var).grid(row=3, column=0, columnspan=3, sticky="w", padx=6, pady=(6, 0))
        ttk.Checkbutton(opcoes, text="Incluir shapefile consolidado (todos os anos)", variable=self.shp_consolidado_fundido_var).grid(row=3, column=3, columnspan=3, sticky="w", padx=6, pady=(6, 0))

        metricas_frame = ttk.LabelFrame(frame, text="Seleção de métricas por nível", padding=8)
        metricas_frame.grid(row=6, column=0, columnspan=3, sticky="nsew", padx=6, pady=6)
        metricas_frame.columnconfigure(0, weight=1)
        metricas_frame.columnconfigure(1, weight=1)
        metricas_frame.columnconfigure(2, weight=1)
        self.listbox_patch_fundido = self._criar_listbox_metricas(metricas_frame, "Patch", PATCH_METRICS, 0)
        self.listbox_class_fundido = self._criar_listbox_metricas(metricas_frame, "Classe", CLASS_METRICS, 1)
        self.listbox_landscape_fundido = self._criar_listbox_metricas(metricas_frame, "Paisagem", LANDSCAPE_METRICS, 2)

        ajuda = (
            "Mesma lógica de recorte e métricas das abas 3 e 6, sem gravar o recorte em disco. Uma tarefa por raster, "
            "em paralelo. Marque 'Usar processos' para escapar do GIL; desmarque para usar threads. "
            "'Salvar Parquet por ano' grava 1 Parquet em cada pasta de ano (em fluxo, baixa memória) — recomendado; "
            "'Gerar Parquet único ao final' junta esses Parquets num só. 'Salvar CSV por feição' é opcional. "
            "'Retomar de onde parou' pula o que já existe (Parquet do ano ou CSVs), com barra/contagem/ETA globais; "
            "se a configuração mudar, a retomada é desativada. Saídas em pastas por ano; o Parquet único fica na raiz. "
            "'Salvar shapefile com métricas' gera, em 'shapefiles_metricas/', um .shp por ano no mesmo formato da aba 9 "
            "(geometria replicada por linha de métrica, campos com até 10 caracteres). "
            "A aba '8. Consolidar' também junta vários Parquets/CSVs num único."
        )
        ttk.Label(frame, text=ajuda, foreground="#555555", justify="left", wraplength=1040).grid(row=7, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 8))

        self.botao_fundido, self.botao_parar_fundido = self._adicionar_botoes_iniciar_parar(
            frame,
            row=8,
            texto_iniciar="Processar em memória e salvar CSVs",
            comando_iniciar=lambda: self._executar_thread(self.botao_fundido, self.botao_parar_fundido, self.stop_fundido_event, self._processar_fundido),
            comando_parar=lambda: self._solicitar_parada(self.stop_fundido_event, self.log_fundido),
        )

        self.label_load_fundido = ttk.Label(frame, text="Carregamento: aguardando arquivos")
        self.label_load_fundido.grid(row=9, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_load_fundido = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_load_fundido.grid(row=10, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 6))
        self.label_proc_fundido = ttk.Label(frame, text="Processamento: aguardando início")
        self.label_proc_fundido.grid(row=11, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_proc_fundido = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_proc_fundido.grid(row=12, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 8))

        ttk.Label(frame, text="Log do processamento:").grid(row=13, column=0, columnspan=3, sticky="w", padx=6)
        self.log_fundido = tk.Text(frame, height=9, wrap="word")
        self.log_fundido.grid(row=14, column=0, columnspan=3, sticky="nsew", padx=6, pady=6)
        barra = ttk.Scrollbar(frame, orient="vertical", command=self.log_fundido.yview)
        barra.grid(row=14, column=3, sticky="ns")
        self.log_fundido.configure(yscrollcommand=barra.set)
        frame.rowconfigure(6, weight=2)
        frame.rowconfigure(14, weight=1)

        self._selecionar_todas_metricas(self.listbox_patch_fundido)
        self._selecionar_todas_metricas(self.listbox_class_fundido)
        self._selecionar_todas_metricas(self.listbox_landscape_fundido)
        self._atualizar_estado_metricas_fundido()

    def _atualizar_estado_metricas_fundido(self):
        controles = [
            (self.listbox_patch_fundido, self.nivel_patch_fundido_var.get()),
            (self.listbox_class_fundido, self.nivel_class_fundido_var.get()),
            (self.listbox_landscape_fundido, self.nivel_landscape_fundido_var.get()),
        ]
        for lb, ativo in controles:
            lb.configure(state="normal" if ativo else "disabled")

    def _acao_selecionar_rasters_fundido(self):
        caminhos = self._selecionar_multiplos_rasters("Selecione um ou mais rasters")
        if caminhos:
            self.rasters_fundido_lista = list(caminhos)
            self.rasters_fundido_texto_var.set(caminhos[0] if len(caminhos) == 1 else f"{len(caminhos)} rasters selecionados")
            if not self.pasta_saida_fundido_var.get():
                self.pasta_saida_fundido_var.set(str(Path(caminhos[0]).parent / "saida_metricas_fundido"))

    def _acao_selecionar_buffers_fundido(self):
        caminhos = self._selecionar_multiplos_shp("Selecione um ou mais shapefiles de buffers")
        if caminhos:
            self.shp_buffers_fundido_lista = list(caminhos)
            self.shp_buffers_fundido_texto_var.set(caminhos[0] if len(caminhos) == 1 else f"{len(caminhos)} shapefiles selecionados")
            if not self.pasta_saida_fundido_var.get():
                self.pasta_saida_fundido_var.set(str(Path(caminhos[0]).parent / "saida_metricas_fundido"))

    def _acao_selecionar_pasta_fundido(self):
        pasta = self._selecionar_pasta("Selecione a pasta de saída dos CSVs")
        if pasta:
            self.pasta_saida_fundido_var.set(pasta)

    def _processar_fundido(self):
        self.log_fundido.after(0, lambda: self._limpar_log(self.log_fundido))
        self._resetar_barras([
            (self.progress_load_fundido, self.label_load_fundido, "Carregamento: aguardando arquivos"),
            (self.progress_proc_fundido, self.label_proc_fundido, "Processamento: aguardando início"),
        ])
        if not self.rasters_fundido_lista:
            raise ValueError("Selecione um ou mais rasters.")
        if not self.shp_buffers_fundido_lista:
            raise ValueError("Selecione um ou mais shapefiles de buffers.")
        pasta_saida = self.pasta_saida_fundido_var.get().strip()
        if not pasta_saida:
            raise ValueError("Selecione a pasta de saída dos CSVs.")

        niveis = {
            "patch": bool(self.nivel_patch_fundido_var.get()),
            "class": bool(self.nivel_class_fundido_var.get()),
            "landscape": bool(self.nivel_landscape_fundido_var.get()),
        }
        metricas_por_nivel = {
            "patch": self._obter_metricas_listbox(self.listbox_patch_fundido),
            "class": self._obter_metricas_listbox(self.listbox_class_fundido),
            "landscape": self._obter_metricas_listbox(self.listbox_landscape_fundido),
        }

        self._log(self.log_fundido, "Iniciando processamento fundido em memória (sem salvar recortes)...")
        resultado = processar_buffers_rasters_em_memoria(
            shp_buffers_lista=self.shp_buffers_fundido_lista,
            rasters_lista=self.rasters_fundido_lista,
            pasta_saida=pasta_saida,
            niveis=niveis,
            metricas_por_nivel=metricas_por_nivel,
            vizinhanca=self.vizinhanca_fundido_var.get(),
            campo_id=self.campo_id_fundido_var.get().strip() or "id_pt",
            campo_diametro=self.campo_diam_fundido_var.get().strip() or "diam_m",
            all_touched=bool(self.all_touched_fundido_var.get()),
            max_workers=self._ler_workers(self.workers_fundido_var),
            usar_processos=bool(self.usar_processos_fundido_var.get()),
            salvar_csv_por_feicao=bool(self.salvar_csv_feicao_fundido_var.get()),
            salvar_parquet_por_ano=bool(self.salvar_parquet_ano_fundido_var.get()),
            gerar_parquet_unico_final=bool(self.parquet_unico_final_fundido_var.get()),
            pular_reproj_metrico=bool(self.pular_reproj_fundido_var.get()),
            retomar=bool(self.retomar_fundido_var.get()),
            log_callback=lambda msg: self._log(self.log_fundido, msg),
            load_callback=lambda acao, msg: self._iniciar_barra_carregamento(self.progress_load_fundido, self.label_load_fundido, msg) if acao == "start" else self._parar_barra_carregamento(self.progress_load_fundido, self.label_load_fundido, msg),
            progress_callback=lambda valor, maximo, msg: self._atualizar_barra_processamento(self.progress_proc_fundido, self.label_proc_fundido, valor, maximo, msg),
            stop_event=self.stop_fundido_event,
        )
        n_csvs = resultado["n_csvs_por_feicao"]
        parquet_unico = resultado.get("parquet_unico")
        self._log(self.log_fundido, "=" * 70)
        self._log(self.log_fundido, f"Rasters processados agora: {resultado['rasters_processados']} | "
                                    f"pulados (já prontos): {resultado.get('rasters_pulados', 0)} | "
                                    f"total: {resultado.get('rasters_total', resultado['rasters_processados'])}")
        self._log(self.log_fundido, f"Feições já existentes (retomadas): {resultado.get('feicoes_ja_existentes', 0)} | "
                                    f"concluídas no total: {resultado.get('feicoes_concluidas', 0)}/{resultado.get('total_global', 0)}")
        self._log(self.log_fundido, f"Anos: {', '.join(resultado['anos'])}")
        self._log(self.log_fundido, f"CSVs por feição gravados agora: {n_csvs}")
        self._log(self.log_fundido, f"Parquets por ano gravados agora: {resultado.get('n_parquet_ano', 0)}")
        self._log(self.log_fundido, f"Parquet único final: {parquet_unico or '—'}")

        # ---- Pós-etapa opcional: shapefiles de métricas (como a aba 9) ----
        shp_info = None
        if bool(self.salvar_shp_fundido_var.get()):
            self._log(self.log_fundido, "-" * 70)
            self._log(self.log_fundido, "Gerando shapefiles de métricas (formato da aba 9)...")
            shp_info = gerar_shapefiles_metricas_fundido(
                shp_buffers_lista=list(self.shp_buffers_fundido_lista),
                pasta_saida=self.pasta_saida_fundido_var.get().strip(),
                campo_id=self.campo_id_fundido_var.get().strip() or "id_pt",
                incluir_consolidado=bool(self.shp_consolidado_fundido_var.get()),
                log_callback=lambda msg: self._log(self.log_fundido, msg),
                load_callback=lambda acao, msg: self._iniciar_barra_carregamento(self.progress_load_fundido, self.label_load_fundido, msg) if acao == "start" else self._parar_barra_carregamento(self.progress_load_fundido, self.label_load_fundido, msg),
                progress_callback=lambda valor, maximo, msg: self._atualizar_barra_processamento(self.progress_proc_fundido, self.label_proc_fundido, valor, maximo, msg),
                stop_event=self.stop_fundido_event,
            )
            self._log(self.log_fundido, f"Shapefiles de métricas: {len(shp_info['shapefiles'])} por ano "
                                        f"em {shp_info['pasta']}")

        linha_shp = (f"\nShapefiles de métricas: {len(shp_info['shapefiles'])} (pasta shapefiles_metricas)"
                     if shp_info else "")
        self.log_fundido.after(0, lambda: messagebox.showinfo(
            "Concluído",
            f"Processamento fundido concluído.\n"
            f"Rasters processados agora: {resultado['rasters_processados']} | pulados: {resultado.get('rasters_pulados', 0)}\n"
            f"Parquets por ano: {resultado.get('n_parquet_ano', 0)} | CSVs por feição: {n_csvs}\n"
            f"Parquet único final: {Path(parquet_unico).name if parquet_unico else '—'}" + linha_shp,
        ))

    # -------------------------------------------------------------------------
    # Aba 8: consolidação rápida das tabelas (CSV / Parquet)
    # -------------------------------------------------------------------------

    def _montar_aba_consolidar(self):
        frame = ttk.Frame(self.aba_consol, padding=10)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)

        ttk.Label(
            frame,
            text="Consolidar todas as tabelas (CSV das subpastas) em uma tabela única",
            font=("Segoe UI", 12, "bold"),
        ).grid(row=0, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 10))

        self._adicionar_linha_pasta(frame, "Pasta de entrada (com subpastas de CSV):", self.pasta_entrada_consol_var, self._acao_selecionar_pasta_entrada_consol, row=1)
        self._adicionar_linha_pasta(frame, "Pasta de saída da tabela final:", self.pasta_saida_consol_var, self._acao_selecionar_pasta_saida_consol, row=2)

        nome = ttk.Frame(frame)
        nome.grid(row=3, column=0, columnspan=3, sticky="ew", padx=6, pady=6)
        ttk.Label(nome, text="Nome base do arquivo:").pack(side="left")
        ttk.Entry(nome, textvariable=self.nome_base_consol_var, width=32).pack(side="left", padx=(4, 0))

        formatos = ttk.LabelFrame(frame, text="Formatos de saída e desempenho", padding=8)
        formatos.grid(row=4, column=0, columnspan=3, sticky="ew", padx=6, pady=6)
        ttk.Checkbutton(formatos, text="CSV único", variable=self.consol_csv_var).grid(row=0, column=0, sticky="w", padx=6, pady=2)
        ttk.Checkbutton(formatos, text="Parquet único (todos os dados)", variable=self.consol_parquet_var).grid(row=1, column=0, sticky="w", padx=6, pady=2)
        ttk.Checkbutton(formatos, text="Parquet por ano (separado/particionado)", variable=self.consol_parquet_part_var).grid(row=2, column=0, sticky="w", padx=6, pady=2)
        ttk.Checkbutton(formatos, text="Entrada em Parquet (juntar Parquets, em vez de ler CSVs)", variable=self.consol_entrada_parquet_var).grid(row=3, column=0, columnspan=3, sticky="w", padx=6, pady=(6, 2))
        ttk.Label(formatos, text="Coluna de partição:").grid(row=2, column=1, sticky="e", padx=(20, 6))
        ttk.Entry(formatos, textvariable=self.consol_coluna_part_var, width=10).grid(row=2, column=2, sticky="w", padx=6)
        ttk.Label(formatos, text="Núcleos (DuckDB):").grid(row=0, column=1, sticky="e", padx=(20, 6))
        ttk.Spinbox(formatos, from_=1, to=(os.cpu_count() or 4) * 2, width=5, textvariable=self.workers_consol_var).grid(row=0, column=2, sticky="w", padx=6)
        ttk.Label(formatos, text="Compressão Parquet:").grid(row=1, column=1, sticky="e", padx=(20, 6))
        ttk.Combobox(formatos, textvariable=self.consol_compressao_var, values=["zstd", "snappy", "gzip"], width=8, state="readonly").grid(row=1, column=2, sticky="w", padx=6)

        ajuda = (
            "Lê todos os CSVs da pasta e de TODAS as subpastas. Se o DuckDB estiver instalado (pip install duckdb), "
            "ele lê todos os arquivos de uma vez, em paralelo (use o campo Núcleos), e grava diretamente o Parquet único, "
            "o Parquet por ano e/ou o CSV — é o caminho mais rápido e eficiente, ideal para o seu volume. Sem DuckDB, "
            "usa um caminho alternativo (CSV por concatenação de bytes e Parquet via pyarrow), mais lento. Parquet costuma "
            "ficar 5–15× menor que o CSV; o particionado por ano permite à app web carregar só o(s) ano(s) filtrado(s). "
            "A leitura mostra o número de arquivos e a estimativa de tempo (lê em lotes para a barra avançar)."
        )
        ttk.Label(frame, text=ajuda, foreground="#555555", justify="left", wraplength=1040).grid(row=5, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 8))

        self.botao_consol, self.botao_parar_consol = self._adicionar_botoes_iniciar_parar(
            frame,
            row=6,
            texto_iniciar="Consolidar tabelas",
            comando_iniciar=lambda: self._executar_thread(self.botao_consol, self.botao_parar_consol, self.stop_consol_event, self._processar_consolidar),
            comando_parar=lambda: self._solicitar_parada(self.stop_consol_event, self.log_consol),
        )

        self.label_load_consol = ttk.Label(frame, text="Carregamento: aguardando")
        self.label_load_consol.grid(row=7, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_load_consol = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_load_consol.grid(row=8, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 6))
        self.label_proc_consol = ttk.Label(frame, text="Processamento: aguardando início")
        self.label_proc_consol.grid(row=9, column=0, columnspan=3, sticky="w", padx=6, pady=(4, 0))
        self.progress_proc_consol = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.progress_proc_consol.grid(row=10, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 8))

        ttk.Label(frame, text="Log do processamento:").grid(row=11, column=0, columnspan=3, sticky="w", padx=6)
        self.log_consol = tk.Text(frame, height=14, wrap="word")
        self.log_consol.grid(row=12, column=0, columnspan=3, sticky="nsew", padx=6, pady=6)
        barra = ttk.Scrollbar(frame, orient="vertical", command=self.log_consol.yview)
        barra.grid(row=12, column=3, sticky="ns")
        self.log_consol.configure(yscrollcommand=barra.set)
        frame.rowconfigure(12, weight=1)

    def _acao_selecionar_pasta_entrada_consol(self):
        pasta = self._selecionar_pasta("Selecione a pasta com os subdiretórios de CSV")
        if pasta:
            self.pasta_entrada_consol_var.set(pasta)
            if not self.pasta_saida_consol_var.get():
                self.pasta_saida_consol_var.set(str(Path(pasta) / "tabela_consolidada"))

    def _acao_selecionar_pasta_saida_consol(self):
        pasta = self._selecionar_pasta("Selecione a pasta de saída da tabela consolidada")
        if pasta:
            self.pasta_saida_consol_var.set(pasta)

    def _processar_consolidar(self):
        self.log_consol.after(0, lambda: self._limpar_log(self.log_consol))
        self._resetar_barras([
            (self.progress_load_consol, self.label_load_consol, "Carregamento: aguardando"),
            (self.progress_proc_consol, self.label_proc_consol, "Processamento: aguardando início"),
        ])
        pasta_entrada = self.pasta_entrada_consol_var.get().strip()
        pasta_saida = self.pasta_saida_consol_var.get().strip()
        if not pasta_entrada:
            raise ValueError("Selecione a pasta de entrada (com os subdiretórios de CSV).")
        if not pasta_saida:
            raise ValueError("Selecione a pasta de saída da tabela final.")

        self._log(self.log_consol, "Iniciando consolidação das tabelas...")
        resultado = consolidar_tabelas(
            pasta_entrada=pasta_entrada,
            pasta_saida=pasta_saida,
            nome_base=self.nome_base_consol_var.get().strip() or "tabela_consolidada",
            gerar_csv=bool(self.consol_csv_var.get()),
            gerar_parquet_unico=bool(self.consol_parquet_var.get()),
            gerar_parquet_por_ano=bool(self.consol_parquet_part_var.get()),
            coluna_particao=self.consol_coluna_part_var.get().strip() or "year",
            compressao=self.consol_compressao_var.get().strip() or "zstd",
            entrada_parquet=bool(self.consol_entrada_parquet_var.get()),
            max_workers=self._ler_workers(self.workers_consol_var),
            log_callback=lambda msg: self._log(self.log_consol, msg),
            load_callback=lambda acao, msg: self._iniciar_barra_carregamento(self.progress_load_consol, self.label_load_consol, msg) if acao == "start" else self._parar_barra_carregamento(self.progress_load_consol, self.label_load_consol, msg),
            progress_callback=lambda valor, maximo, msg: self._atualizar_barra_processamento(self.progress_proc_consol, self.label_proc_consol, valor, maximo, msg),
            stop_event=self.stop_consol_event,
        )
        saidas = resultado["saidas"]
        self._log(self.log_consol, "=" * 70)
        self._log(self.log_consol, f"Motor: {resultado.get('motor', '—')}")
        if resultado.get("arquivos_lidos") is not None:
            self._log(self.log_consol, f"Arquivos lidos: {resultado['arquivos_lidos']}")
        if resultado.get("linhas"):
            self._log(self.log_consol, f"Linhas consolidadas: {resultado['linhas']}")
        for tipo, caminho in saidas.items():
            self._log(self.log_consol, f"Saída {tipo}: {caminho}")
        resumo = "\n".join(f"{t}: {Path(c).name}" for t, c in saidas.items()) or "nenhuma saída gerada"
        self.log_consol.after(0, lambda: messagebox.showinfo(
            "Concluído",
            f"Consolidação concluída (motor: {resultado.get('motor','—')}).\n{resumo}",
        ))


# =============================================================================
# Inicialização
# =============================================================================

if __name__ == "__main__":
    import multiprocessing as _mp
    _mp.freeze_support()  # necessário para ProcessPoolExecutor (spawn) no Windows/exe
    app = AplicacaoGradeBuffer()
    app.mainloop()