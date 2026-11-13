#!/usr/bin/env python3
"""
geoland_indice.py
=================

Indice de Estabilidade da Paisagem (IEP), portado de
script_calculo_estabilidade.R mantendo a mesma base de calculo:

1. Para cada paisagem (buffer) x classe x escala (raio), o coeficiente de
   variacao (CV = desvio padrao amostral / media) de cada metrica ao longo
   dos anos: CA, PLAND, NP, LPI e ENN_MN. Valores ausentes sao ignorados;
   com 0 ou 1 valor, ou media zero, o CV e 0.
2. Cada CV e normalizado entre 0 e 1 (min-max) considerando todas as
   paisagens, classes e escalas da tabela; sem variacao, o normalizado e 0.
3. Variabilidade por componente (media dos CVs normalizados):
   composicao (CA, PLAND), configuracao (NP, LPI) e conectividade (ENN_MN).
4. Variabilidade global V_u = media ponderada dos componentes (1/3 cada) e
   IEP_u = 1 - V_u: quanto mais perto de 1, mais estavel a classe no periodo.

O indice precisa da serie temporal completa; por isso e calculado depois
das metricas de todos os anos.
"""

from __future__ import annotations

import re
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

# Metricas do indice e os nomes equivalentes nas tabelas do GeoLand
# (landscapemetrics em minusculas; PyLandStats por extenso).
METRICAS_IEP = ("CA", "PLAND", "NP", "LPI", "ENN_MN")
APELIDOS = {
    "CA": ("ca", "total_area"),
    "PLAND": ("pland", "proportion_of_landscape"),
    "NP": ("np", "number_of_patches"),
    "LPI": ("lpi", "largest_patch_index"),
    "ENN_MN": ("enn_mn", "euclidean_nearest_neighbor_mn"),
}
COMPONENTES = {"comp": ("CA", "PLAND"), "conf": ("NP", "LPI"), "conn": ("ENN_MN",)}
PESOS_COMPONENTES = {"comp": 1 / 3, "conf": 1 / 3, "conn": 1 / 3}
# Metricas a pedir a cada motor para o indice poder ser calculado.
METRICAS_PYLANDSTATS = tuple(APELIDOS[m][1] for m in METRICAS_IEP)
FUNCOES_LSM = tuple("lsm_c_" + APELIDOS[m][0] for m in METRICAS_IEP)
# Colunas de resultado que podem virar grafico, com o nome por extenso.
NOMES_INDICE = {
    "IEP_u": "Índice de Estabilidade da Paisagem (IEP)",
    "V_u": "Variabilidade global (V)",
    "V_comp": "Variabilidade da composição (CA, PLAND)",
    "V_conf": "Variabilidade da configuração (NP, LPI)",
    "V_conn": "Variabilidade da conectividade (ENN_MN)",
}

PADRAO_RECORTE = re.compile(r"(?i)^(?:.+_)?clip_id_(?P<buffer>.+)_diam_(?P<raio>[^_]+)_(?P<ano>[^_]+)$")


def cv(valores) -> float:
    """CV robusto do script: ignora ausentes; 0 se houver ate 1 valor ou se a media for zero."""
    x = pd.to_numeric(pd.Series(valores), errors="coerce").dropna().to_numpy(dtype=float)
    if len(x) <= 1:
        return 0.0
    media = x.mean()
    if media == 0:
        return 0.0
    return float(x.std(ddof=1) / media)


def colunas_das_metricas(df: pd.DataFrame) -> Dict[str, str]:
    """Acha, na tabela, a coluna de cada metrica do indice ({'CA': 'ca', ...}); so as encontradas."""
    minusculas = {str(c).strip().lower(): c for c in df.columns}
    achadas = {}
    for metrica, apelidos in APELIDOS.items():
        for nome in (metrica.lower(),) + apelidos:
            if nome in minusculas:
                achadas[metrica] = minusculas[nome]
                break
    return achadas


def identificar_paisagem_e_escala(df: pd.DataFrame, coluna_nome: str = "source_stem") -> pd.DataFrame:
    """
    Cria as colunas "Buffer" (paisagem) e "Raio" (escala) a partir do nome do
    raster. Nos recortes do GeoLand (clip_id_<ID>_diam_<RAIO>_<ANO>) usa esses
    campos; em outros nomes segue o script original: paisagem = texto ate o
    primeiro "_" e escala = primeiro numero entre dois "_".
    """
    df = df.copy()
    nomes = df[coluna_nome].astype(str).str.replace(r"\.tiff?$", "", regex=True, case=False)
    nomes = nomes.str.replace(r"^.*[\\/]", "", regex=True)
    partes = nomes.str.extract(PADRAO_RECORTE)
    if partes["buffer"].notna().mean() >= 0.9:
        df["Buffer"], escala = partes["buffer"], partes["raio"]
    else:
        df["Buffer"] = nomes.str.extract(r"^([^_]+)")[0]
        escala = nomes.str.extract(r"(?<=_)([0-9]+)(?=_)")[0]
    df["Raio"] = pd.to_numeric(escala, errors="coerce")
    return df


def calcular_iep(df: pd.DataFrame, paisagem: str = "Buffer", classe: str = "class_val",
                 escala: Optional[str] = "Raio", ano: Optional[str] = "year",
                 metricas: Optional[Dict[str, str]] = None,
                 pesos_componentes: Optional[Dict[str, float]] = None) -> pd.DataFrame:
    """
    Calcula o IEP por paisagem x classe x escala a partir de uma tabela com
    uma linha por paisagem, classe, escala e ano.

    Devolve: as colunas de identificacao, n_anos, CV_<metrica>,
    NORM_CV_<metrica>, V_comp, V_conf, V_conn, V_u e IEP_u.
    """
    metricas = metricas or colunas_das_metricas(df)
    faltam = [m for m in METRICAS_IEP if m not in metricas]
    if faltam:
        raise ValueError("Para o índice de estabilidade faltam na tabela as métricas de classe: " + ", ".join(faltam))
    for coluna in (paisagem, classe):
        if coluna not in df.columns:
            raise ValueError(f"Coluna necessária para o índice não encontrada: {coluna}")
    ids = [paisagem, classe] + ([escala] if escala and escala in df.columns else [])
    dados = df.dropna(subset=[classe])

    grupos = dados.groupby(ids, sort=True, dropna=False)
    resultado = grupos.size().rename("n_linhas").reset_index()
    if ano and ano in dados.columns:
        resultado["n_anos"] = grupos[ano].nunique().to_numpy()
    else:
        resultado["n_anos"] = resultado["n_linhas"]
    resultado = resultado.drop(columns="n_linhas")
    for metrica in METRICAS_IEP:  # 1) CV temporal
        resultado["CV_" + metrica] = grupos[metricas[metrica]].agg(cv).to_numpy()
    for metrica in METRICAS_IEP:  # 2) normalizacao min-max entre todos os grupos
        valores = resultado["CV_" + metrica]
        minimo, maximo = valores.min(), valores.max()
        if not np.isfinite(minimo) or not np.isfinite(maximo) or maximo - minimo == 0:
            resultado["NORM_CV_" + metrica] = 0.0
        else:
            resultado["NORM_CV_" + metrica] = (valores - minimo) / (maximo - minimo)
    for componente, mets in COMPONENTES.items():  # 3) media por componente (ausentes ignorados; tudo ausente = 0)
        resultado["V_" + componente] = resultado[["NORM_CV_" + m for m in mets]].mean(axis=1, skipna=True).fillna(0.0)
    pesos = dict(pesos_componentes or PESOS_COMPONENTES)  # 4) variabilidade global e indice
    total = sum(pesos.values())
    resultado["V_u"] = sum(pesos[c] / total * resultado["V_" + c] for c in COMPONENTES)
    resultado["IEP_u"] = 1 - resultado["V_u"]
    return resultado


def calcular_iep_de_tabela_de_metricas(df: pd.DataFrame) -> pd.DataFrame:
    """IEP a partir de uma tabela geral de metricas do GeoLand (nivel de classe)."""
    minusculas = {str(c).lower(): c for c in df.columns}
    if "level" in minusculas:
        df = df[df[minusculas["level"]].astype(str).str.lower().isin(("class", "classe"))]
    if "Buffer" not in df.columns or "Raio" not in df.columns:
        origem = next((minusculas[n] for n in ("source_stem", "source_file", "lid") if n in minusculas), None)
        if origem is None:
            raise ValueError("A tabela não tem a coluna com o nome do raster (source_stem) para identificar buffer e raio.")
        df = identificar_paisagem_e_escala(df, origem)
    classe = next((minusculas[n] for n in ("class_val", "class", "classe", "classe_num") if n in minusculas), None)
    if classe is None:
        raise ValueError("A tabela não tem a coluna da classe (class_val).")
    ano = next((minusculas[n] for n in ("year", "ano") if n in minusculas), None)
    return calcular_iep(df, "Buffer", classe, "Raio", ano)


def salvar_iep_do_csv_geral(csv_geral: str, destino: str) -> Dict[str, object]:
    """Lê o CSV geral de métricas de um motor, calcula o IEP e grava o CSV do índice."""
    tabela = pd.read_csv(csv_geral, low_memory=False)
    indice = calcular_iep_de_tabela_de_metricas(tabela)
    indice.to_csv(destino, index=False, encoding="utf-8-sig")
    return {"arquivo": destino, "linhas": len(indice), "anos": int(indice["n_anos"].max()) if len(indice) else 0,
            "iep_medio": float(indice["IEP_u"].mean()) if len(indice) else float("nan")}
