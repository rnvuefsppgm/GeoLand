#!/usr/bin/env python3
"""
geoland_graficos.py
===================

Aba "10. Gráficos" do GeoLand: lê uma tabela qualquer (CSV, TXT ou Parquet),
estrutura os campos automaticamente e gera gráficos de linhas e pontos no
padrão visual do ggplot2, um a um (pré-visualização) ou em lote.

O arquivo tem duas partes independentes:

1. Motor (sem interface): `ler_tabela`, `estruturar_tabela`, `ler_legenda`,
   `Plano` e `gerar_graficos`. Pode ser usado em scripts.
2. Interface: `AbaGraficos`, montada pelo aplicativo principal.

Modelo de dados: tabela "larga", com uma coluna para o eixo X (ex.: year),
uma coluna opcional que separa as linhas (ex.: class_val), colunas que
identificam cada gráfico (ex.: Buffer e Raio) e uma coluna numérica por
métrica. Nas tabelas do GeoLand, Buffer e Raio são criados a partir do
nome do raster (clip_id_<ID>_diam_<DIAM>_<ANO>).
"""

from __future__ import annotations

import base64
import io
import math
import re
import textwrap
import threading
import unicodedata
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from nomes_metricas_lsm import NOMES_LSM
except ImportError:  # arquivo gerado; sem ele o rótulo é o nome da coluna
    NOMES_LSM = {}
from geoland_indice import NOMES_INDICE, calcular_iep_de_tabela_de_metricas

# =============================================================================
# Constantes
# =============================================================================

PADRAO_RECORTE = re.compile(r"(?i)^(?:.+_)?clip_id_(?P<buffer>.+)_diam_(?P<diam>[^_]+)_(?P<ano>[^_]+)$")
# Leitura flexível do nome do arquivo (aceita prefixos e sufixos, ex.: Reclass_clip_id_APA_1_diam_500_2025,
# clip_id_100_diam_1_ano_2024_C10_1_reclass_2024_EPSG32724): ID do buffer, diâmetro/raio e ano.
PADRAO_ID_DIAM = re.compile(r"(?i)(?:^|_)id_(?P<buffer>.+?)_(?:diam|diametro|raio|radius)_?(?P<diam>\d+(?:[.,]\d+)?)(?P<resto>.*)$")
PADRAO_ANO = re.compile(r"(?i)(?:^|_)(?:ano_?|year_?)?(?P<ano>(?:19|20)\d{2})(?=_|$)")


def extrair_campos_do_nome(nomes: pd.Series) -> pd.DataFrame:
    """Colunas buffer, diam e ano a partir do nome do arquivo (sem pasta e sem .tif).
    O ano é o 1º ano (19xx/20xx) depois do diâmetro; nomes fora do padrão ficam vazios."""
    base = nomes.astype(str).str.replace(r"^.*[\\/]", "", regex=True).str.replace(r"\.tiff?$", "", regex=True, case=False)
    partes = base.str.extract(PADRAO_ID_DIAM)
    partes["diam"] = partes["diam"].str.replace(",", ".", regex=False)
    partes["ano"] = partes["resto"].fillna("").str.extract(PADRAO_ANO)["ano"]
    return partes[["buffer", "diam", "ano"]]
NOMES_X = ("year", "ano", "anos", "data", "date", "tempo", "time", "periodo")
NOMES_SERIE = ("class_val", "class", "classe", "class_name", "classe_nome", "categoria", "grupo", "group")
ARQUIVOS_LEGENDA = ("cores.txt", "cores.csv", "legenda.txt", "legenda.csv")
SEM_VALOR = "(sem valor)"
NA_INTEIRO_R = -2147483648  # NA_integer_ do R
NENHUMA = "(nenhuma)"
AGREGACOES = {"Média": "mean", "Mediana": "median", "Soma": "sum", "Mínimo": "min", "Máximo": "max"}
MARCADORES = ("o", "^", "s", "D", "v", "P", "X", "<", ">", "p", "h", "*")
FONTES = ["Arial", "Helvetica", "DejaVu Sans"]
SUFIXOS_LSM = {"_mn": "mean", "_sd": "standard deviation", "_cv": "coefficient of variation"}
ROTULOS_PADRAO = {"year": "Ano", "ano": "Ano", "class_val": "Classe", "class": "Classe", "classe": "Classe"}

# Temas equivalentes aos do ggplot2: fundo do painel, cor da grade e contorno.
TEMA_PADRAO = "Cinza do ggplot2 (theme_gray)"
TEMAS = {
    TEMA_PADRAO: {"painel": "#EBEBEB", "grade": "#FFFFFF", "contorno": ""},
    "Branco com grade (theme_bw)": {"painel": "#FFFFFF", "grade": "#EBEBEB", "contorno": "caixa"},
    "Clássico sem grade (theme_classic)": {"painel": "#FFFFFF", "grade": "", "contorno": "eixos"},
}
LEGENDAS = ("À direita", "Abaixo", "Sem legenda")

# Palavras que indicam que uma coluna numérica é uma dimensão (não uma métrica).
_DIMENSAO_EXATA = {"id", "fid", "cod", "codigo", "code", "ano", "year", "level", "nivel", "rule", "engine"}
_DIMENSAO_INICIO = ("class", "diam", "zona", "zone", "buffer", "grup", "group", "tipo", "type", "categ", "regra")
_IDENTIFICADOR = {"id", "fid", "objectid", "gid"}


# =============================================================================
# Utilidades
# =============================================================================

def nome_seguro(texto) -> str:
    """Texto sem acentos nem símbolos, para nomes de arquivo."""
    texto = unicodedata.normalize("NFKD", str(texto)).encode("ascii", "ignore").decode("ascii")
    return re.sub(r"_+", "_", re.sub(r"[^\w\-.]+", "_", texto)).strip("_") or "sem_nome"


def texto_valor(valor) -> str:
    """Valor como texto de legenda: 3.0 -> '3'; vazio -> '(sem valor)'."""
    if valor is None or (isinstance(valor, (float, np.floating)) and math.isnan(valor)):
        return SEM_VALOR
    if isinstance(valor, (float, np.floating)) and float(valor).is_integer():
        return str(int(valor))
    if isinstance(valor, np.integer):
        return str(int(valor))
    return str(valor)


def chave_natural(texto) -> list:
    """Chave de ordenação em que 'APA_2' vem antes de 'APA_10'."""
    return [int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", str(texto))]


def ordenar_valores(valores) -> list:
    valores = [v for v in valores if texto_valor(v) != SEM_VALOR]
    if all(isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool) for v in valores):
        return sorted(valores)
    return sorted(valores, key=lambda v: chave_natural(texto_valor(v)))


def nome_bonito(coluna: str) -> str:
    """'number_of_patches' -> 'Number of patches'."""
    chave = str(coluna).strip().lower()
    if chave in ROTULOS_PADRAO:
        return ROTULOS_PADRAO[chave]
    texto = str(coluna).replace("_", " ").strip()
    return texto[:1].upper() + texto[1:]


def nome_metrica(coluna: str) -> str:
    """Nome por extenso da métrica: siglas do landscapemetrics viram
    'Number of patches (NP)'; outras colunas, 'Number of patches'."""
    if str(coluna) in NOMES_INDICE:
        return NOMES_INDICE[str(coluna)]
    chave = str(coluna).strip().lower()
    if chave in NOMES_LSM:
        nome = NOMES_LSM[chave]
        for sufixo, estatistica in SUFIXOS_LSM.items():
            if chave.endswith(sufixo):
                nome = f"{nome}, {estatistica}"
                break
        return f"{nome[:1].upper()}{nome[1:]} ({str(coluna).upper()})"
    return nome_bonito(coluna)


def _inteiros(serie: pd.Series) -> bool:
    valores = serie.dropna()
    if valores.empty or not pd.api.types.is_numeric_dtype(valores) or pd.api.types.is_bool_dtype(valores):
        return False
    return bool(np.all(np.mod(valores.to_numpy(dtype=float), 1) == 0))


def _parece_dimensao(coluna: str) -> bool:
    partes = [p for p in re.split(r"[^a-z0-9]+", str(coluna).lower()) if p]
    return any(p in _DIMENSAO_EXATA or p.startswith(_DIMENSAO_INICIO) for p in partes)


# ----------------------------- cores -----------------------------------------

def _hcl_para_hex(h: float, c: float, l: float) -> str:
    """Cor HCL (CIELUV polar, branco D65) em hexadecimal sRGB, como grDevices::hcl do R."""
    xn, yn, zn = 95.047, 100.0, 108.883
    u, v = c * math.cos(math.radians(h)), c * math.sin(math.radians(h))
    y = yn * ((l + 16) / 116) ** 3 if l > 8 else yn * l / 903.3
    un, vn = 4 * xn / (xn + 15 * yn + 3 * zn), 9 * yn / (xn + 15 * yn + 3 * zn)
    ul, vl = u / (13 * l) + un, v / (13 * l) + vn
    x, z = 9 * y * ul / (4 * vl), (12 - 3 * ul - 20 * vl) * y / (4 * vl)
    x, y, z = x / 100, y / 100, z / 100
    canais = (3.240479 * x - 1.537150 * y - 0.498535 * z,
              -0.969256 * x + 1.875992 * y + 0.041556 * z,
              0.055648 * x - 0.204043 * y + 1.057311 * z)
    saida = []
    for canal in canais:
        canal = 12.92 * canal if canal <= 0.0031308 else 1.055 * canal ** (1 / 2.4) - 0.055
        saida.append(int(round(255 * min(1.0, max(0.0, canal)))))
    return "#%02X%02X%02X" % tuple(saida)


def paleta_hue(n: int) -> List[str]:
    """Paleta discreta padrão do ggplot2 (scales::hue_pal): n matizes igualmente espaçados."""
    return [_hcl_para_hex((15 + 360.0 * i / n) % 360, 100, 65) for i in range(max(0, n))]


def _rgb(cor: str) -> Tuple[float, float, float]:
    cor = cor.lstrip("#")
    return tuple(int(cor[i:i + 2], 16) / 255 for i in (0, 2, 4))


def _clara(cor: str) -> bool:
    """Cores muito claras somem no fundo do gráfico; elas ganham um contorno."""
    r, g, b = _rgb(cor)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.78


def _escurecer(cor: str, fator: float = 0.62) -> str:
    return "#%02X%02X%02X" % tuple(int(round(255 * canal * fator)) for canal in _rgb(cor))


# =============================================================================
# Leitura e estruturação da tabela
# =============================================================================

def ler_tabela(caminho: str) -> pd.DataFrame:
    """Lê CSV/TXT (detecta separador, vírgula decimal e codificação) ou Parquet."""
    if Path(caminho).suffix.lower() in (".parquet", ".pq"):
        return pd.read_parquet(caminho)
    amostra, codificacao = "", "utf-8-sig"
    for codificacao in ("utf-8-sig", "latin-1"):
        try:
            with open(caminho, "r", encoding=codificacao, newline="") as fh:
                amostra = fh.read(65536)
            break
        except UnicodeDecodeError:
            continue
    primeira = amostra.splitlines()[0] if amostra else ""
    separador = max((",", ";", "\t", "|"), key=primeira.count)
    decimal = "," if separador != "," and re.search(r"\d,\d", amostra) and not re.search(r"\d\.\d", amostra) else "."
    return pd.read_csv(caminho, sep=separador, decimal=decimal, encoding=codificacao, low_memory=False)


def ler_legenda(caminho: str) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
    """
    Lê um arquivo de cores e nomes das séries: uma linha por valor, no formato
    `valor, "#rrggbb", nome` (a cor ou o nome podem faltar; o nome pode conter
    vírgulas). Devolve {valor em texto: (cor, nome)}.
    """
    texto = ""
    for codificacao in ("utf-8-sig", "latin-1"):
        try:
            texto = Path(caminho).read_text(encoding=codificacao)
            break
        except UnicodeDecodeError:
            continue
    legenda: Dict[str, Tuple[Optional[str], Optional[str]]] = {}
    for linha in texto.splitlines():
        partes = [p.strip().strip('"').strip() for p in re.split(r"[,;\t]", linha.strip(), maxsplit=2)]
        if len(partes) < 2 or not partes[0]:
            continue
        cor = nome = None
        if re.fullmatch(r"#[0-9A-Fa-f]{6}", partes[1]):
            cor, nome = partes[1].upper(), (partes[2] if len(partes) > 2 and partes[2] else None)
        else:
            nome = re.split(r"[,;\t]", linha.strip(), maxsplit=1)[1].strip().strip('"').strip() or None
        try:
            chave = texto_valor(float(partes[0]))
        except ValueError:
            chave = partes[0]
        legenda[chave] = (cor, nome)
    return legenda


def salvar_legenda(caminho: str, legenda: Dict[str, Tuple[Optional[str], Optional[str]]]) -> None:
    """Grava o dicionário de cores e nomes no formato lido por ler_legenda: valor, "#rrggbb", nome."""
    linhas = []
    for chave, (cor, nome) in legenda.items():
        partes = [str(chave)] + ([f'"{cor}"'] if cor else []) + ([str(nome)] if nome else [])
        if len(partes) > 1:
            linhas.append(", ".join(partes))
    Path(caminho).write_text("\n".join(linhas) + "\n", encoding="utf-8")


# Paletas para preencher o dicionário de cores de uma vez (None = matizes do ggplot2).
PALETAS = {"ggplot2 (matizes)": None, "Viridis": "viridis", "Turbo (arco-íris)": "turbo", "Espectral": "Spectral",
           "Terreno": "gist_earth", "Qualitativa (20 cores)": "tab20"}


def cores_da_paleta(nome: str, n: int) -> List[str]:
    """n cores tiradas de uma paleta: ao longo do contínuo de cores ou, na qualitativa, em sequência."""
    mapa = PALETAS.get(nome)
    if not mapa:
        return paleta_hue(n)
    import matplotlib

    cmap = matplotlib.colormaps[mapa]
    qualitativa = getattr(cmap, "N", 256) <= 20
    cores = []
    for i in range(n):
        r, g, b = cmap(i % cmap.N if qualitativa else i / max(1, n - 1))[:3]
        cores.append("#%02X%02X%02X" % (int(round(255 * r)), int(round(255 * g)), int(round(255 * b))))
    return cores


@dataclass
class Estrutura:
    """Campos de uma tabela, como a aba os entende."""
    dimensoes: List[str] = field(default_factory=list)   # colunas de categoria (filtros, séries, painéis)
    metricas: List[str] = field(default_factory=list)    # colunas numéricas que podem ir para o eixo Y
    x: Optional[str] = None
    serie: Optional[str] = None
    painel: List[str] = field(default_factory=list)      # colunas que definem "um gráfico por"
    derivadas: List[str] = field(default_factory=list)   # colunas criadas a partir do nome do arquivo
    origem_derivadas: Optional[str] = None
    filtros_iniciais: Dict[str, list] = field(default_factory=dict)


def estruturar_tabela(df: pd.DataFrame) -> Tuple[pd.DataFrame, Estrutura]:
    """
    Descobre o papel de cada coluna e sugere a estrutura do gráfico.

    - Se alguma coluna de texto traz nomes no padrão dos recortes do GeoLand
      (clip_id_<ID>_diam_<DIAM>_<ANO>), cria as colunas "Buffer" e "Raio"
      (e "Ano", se a tabela não tiver ano).
    - Dimensões: colunas de texto com até 2000 valores e colunas de números
      inteiros com até 100 valores cujo nome indica categoria (ano, classe...).
    - Métricas: as demais colunas numéricas que variam.
    - Sugestões: eixo X = ano; linhas = classe; um gráfico por Buffer e
      Raio (ou pelas demais colunas de texto, em tabelas genéricas).
    """
    df = df.copy()
    est = Estrutura()
    minusculas = {str(c).lower(): c for c in df.columns}

    # 0) O R grava o "sem valor" de inteiros como -2147483648 (ex.: class_val nas linhas
    #    de nível paisagem do landscapemetrics): vira ausente, para não parecer uma classe.
    for coluna in df.columns:
        if pd.api.types.is_integer_dtype(df[coluna]) and (df[coluna] == NA_INTEIRO_R).any():
            df[coluna] = df[coluna].where(df[coluna] != NA_INTEIRO_R)

    # 1) Campos embutidos no nome do arquivo de origem.
    for coluna in df.columns:
        if not (pd.api.types.is_object_dtype(df[coluna]) or pd.api.types.is_string_dtype(df[coluna])):
            continue
        nomes = df[coluna].dropna()
        if nomes.empty or extrair_campos_do_nome(nomes.head(2000))["diam"].notna().mean() < 0.9:
            continue
        partes = extrair_campos_do_nome(df[coluna].where(df[coluna].notna(), ""))
        novas = {"Buffer": partes["buffer"], "Raio": partes["diam"]}
        if partes["ano"].notna().any() and not any(n in minusculas for n in NOMES_X):
            novas["Ano"] = partes["ano"]
        for nome, valores in novas.items():
            if nome in df.columns:
                continue
            numerico = pd.to_numeric(valores, errors="coerce")
            df[nome] = numerico if numerico.notna().sum() == valores.notna().sum() else valores
            est.derivadas.append(nome)
        est.origem_derivadas = str(coluna)
        break

    # 2) Dimensões e métricas.
    n_linhas = max(1, len(df))
    for coluna in df.columns:
        serie = df[coluna]
        n = serie.nunique(dropna=True)
        numerica = pd.api.types.is_numeric_dtype(serie) and not pd.api.types.is_bool_dtype(serie)
        if not numerica:
            if 1 <= n <= 2000 and (n <= 50 or n < 0.5 * n_linhas):
                est.dimensoes.append(coluna)
        elif coluna in est.derivadas or (_inteiros(serie) and n <= 100 and _parece_dimensao(coluna)):
            est.dimensoes.append(coluna)
        elif n >= 2 and not (set(re.split(r"[^a-z0-9]+", str(coluna).lower())) & _IDENTIFICADOR):
            est.metricas.append(coluna)
    est.dimensoes.sort(key=lambda c: c not in est.derivadas)  # campos derivados primeiro

    # 3) Sugestões.
    minusculas = {str(c).lower(): c for c in df.columns}
    est.x = next((minusculas[n] for n in NOMES_X if n in minusculas), None)
    if est.x is None:
        est.x = next((c for c in est.dimensoes if pd.api.types.is_numeric_dtype(df[c])), None) or \
            (est.metricas[0] if est.metricas else None)
    est.serie = next((minusculas[n] for n in NOMES_SERIE if n in minusculas and minusculas[n] in est.dimensoes), None)
    variam = [c for c in est.dimensoes if df[c].nunique(dropna=True) >= 2 and c != est.x]
    textos = [c for c in variam if not pd.api.types.is_numeric_dtype(df[c]) and c not in est.derivadas]
    if est.serie is None:  # tabela genérica: a última coluna de texto com poucos valores vira as linhas
        est.serie = next((c for c in reversed(textos) if df[c].nunique(dropna=True) <= 20), None)
    if "Buffer" in est.derivadas:
        est.painel = [c for c in ("Buffer", "Raio") if c in est.derivadas and df[c].nunique(dropna=True) >= 2]
    else:
        est.painel = [c for c in textos if c != est.serie and df[c].nunique(dropna=True) <= 200 and
                      str(c).lower() not in ("level", "nivel", "engine")]
    nivel = next((minusculas[n] for n in ("level", "nivel") if n in minusculas), None)
    if nivel is not None and est.serie is not None and df[nivel].nunique(dropna=True) >= 2:
        classe = [v for v in df[nivel].dropna().unique() if str(v).lower() in ("class", "classe")]
        if classe:  # as linhas são classes: começa mostrando só o nível de classe
            est.filtros_iniciais[nivel] = classe
    return df, est


# =============================================================================
# Desenho
# =============================================================================

@dataclass
class OpcoesGrafico:
    tema: str = TEMA_PADRAO
    largura_cm: float = 18.0
    altura_cm: float = 10.0
    fonte_pt: float = 11.0
    legenda: str = LEGENDAS[0]
    marcadores_distintos: bool = True
    ligar_lacunas: bool = False   # False: a linha é interrompida onde a série não tem dado
    virgula_decimal: bool = False
    usar_cores_do_arquivo: bool = True
    titulo: str = "{painel}"
    subtitulo: str = ""
    rotulo_x: str = ""          # vazio = nome da coluna
    rotulo_y: str = "{metrica}"
    titulo_legenda: str = ""    # vazio = nome da coluna


_MPL = {}


def _mpl() -> dict:
    """Importa o matplotlib só quando o primeiro gráfico é pedido (a janela abre mais rápido)."""
    if not _MPL:
        import matplotlib
        import matplotlib.backends.backend_pdf  # noqa: F401  (formatos vetoriais, usados por savefig)
        import matplotlib.backends.backend_svg  # noqa: F401
        from matplotlib import patheffects, ticker
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure

        matplotlib.rcParams.update({
            "font.family": "sans-serif", "font.sans-serif": FONTES,
            "pdf.fonttype": 42, "svg.fonttype": "none", "axes.unicode_minus": True,
        })
        _MPL.update(Figure=Figure, Canvas=FigureCanvasAgg, ticker=ticker, pe=patheffects)
    return _MPL


def desenhar_grafico(series: Sequence[dict], opcoes: OpcoesGrafico, categorias_x: Optional[List[str]] = None,
                     x_inteiro: bool = False, titulo: str = "", subtitulo: str = "", rotulo_x: str = "",
                     rotulo_y: str = "", titulo_legenda: str = "", marcas_x: Optional[Sequence[float]] = None):
    """
    Desenha um gráfico de linhas e pontos no padrão do ggplot2 e devolve a Figure.

    `series` é uma lista de dicionários {"rotulo", "cor", "marcador", "x", "y"};
    com `categorias_x`, o eixo X é categórico e "x" traz posições (0, 1, 2...).
    """
    m = _mpl()
    ticker, pe = m["ticker"], m["pe"]
    tema = TEMAS.get(opcoes.tema, TEMAS[TEMA_PADRAO])
    base = float(opcoes.fonte_pt)
    largura_pt, altura_pt = opcoes.largura_cm / 2.54 * 72, opcoes.altura_cm / 2.54 * 72

    fig = m["Figure"](figsize=(opcoes.largura_cm / 2.54, opcoes.altura_cm / 2.54), layout="constrained",
                      facecolor="white")
    m["Canvas"](fig)
    fig.get_layout_engine().set(w_pad=0.08, h_pad=0.08)
    ax = fig.add_subplot()
    ax.set_facecolor(tema["painel"])
    ax.set_axisbelow(True)
    if tema["grade"]:
        ax.grid(True, which="major", color=tema["grade"], linewidth=1.07)
        ax.grid(True, which="minor", color=tema["grade"], linewidth=0.53)
    else:
        ax.grid(False)
    for lado, borda in ax.spines.items():
        borda.set_visible(tema["contorno"] == "caixa" or (tema["contorno"] == "eixos" and lado in ("left", "bottom")))
        borda.set_color("#333333")
        borda.set_linewidth(0.8)
    ax.tick_params(axis="both", which="major", length=2.75, width=0.8, color="#333333", labelcolor="#4D4D4D",
                   labelsize=base * 0.8)
    ax.tick_params(axis="both", which="minor", length=0)

    for serie in series:
        cor = serie["cor"]
        extras = {}
        if _clara(cor):  # contorno cinza para a linha clara não sumir no fundo
            extras["path_effects"] = [pe.Stroke(linewidth=2.1, foreground="#8C8C8C"), pe.Normal()]
        ax.plot(serie["x"], serie["y"], color=cor, linewidth=1.1, marker=serie["marcador"], markersize=4.3,
                markeredgecolor=_escurecer(cor), markeredgewidth=0.5, label=serie["rotulo"], zorder=3,
                solid_joinstyle="round", **extras)

    def formatar(valor, _posicao=None):
        texto = f"{valor:.10g}"
        return texto.replace(".", ",") if opcoes.virgula_decimal else texto

    if categorias_x is not None:
        ax.set_xticks(range(len(categorias_x)), categorias_x)
        if len(categorias_x) > 8:
            ax.tick_params(axis="x", labelrotation=45)
            for rotulo in ax.get_xticklabels():
                rotulo.set_horizontalalignment("right")
    else:
        if marcas_x is not None:  # poucos valores de X (ex.: raios 500, 1000, 1500): uma marca em cada um
            ax.set_xticks(list(marcas_x))
        else:
            if x_inteiro:
                ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=9, integer=True, min_n_ticks=2))
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(formatar))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6, min_n_ticks=3))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(formatar))
    ax.margins(x=0.04, y=0.06)

    # Legenda abaixo: quantas colunas cabem na largura e quanto da altura ela ocupa.
    tamanho = base * 0.8
    com_legenda = opcoes.legenda != "Sem legenda" and any(s["rotulo"] for s in series) and \
        not (len(series) == 1 and series[0]["rotulo"] in ("", SEM_VALOR))
    colunas_abaixo, altura_legenda = 1, 0.0
    if com_legenda and opcoes.legenda == "Abaixo":
        maior = max(len(str(s["rotulo"])) for s in series)
        colunas_abaixo = max(1, min(len(series), int((largura_pt - 20) / (26 + maior * tamanho * 0.47))))
        altura_legenda = math.ceil(len(series) / colunas_abaixo) * tamanho * 1.6 + (base * 1.7 if titulo_legenda else 0)

    # Rótulos longos são quebrados para caber na altura/largura do gráfico.
    ax.set_xlabel(rotulo_x, fontsize=base, color="black", labelpad=5)
    altura_eixos = max(60.0, altura_pt - 75 - altura_legenda)
    ax.set_ylabel(textwrap.fill(rotulo_y, max(14, int(altura_eixos / (base * 0.5)))), fontsize=base,
                  color="black", labelpad=6)
    largura_texto = max(30, int(largura_pt * 0.62 / (base * 0.55)))
    deslocamento = 6.0
    if subtitulo:
        ax.annotate(textwrap.fill(subtitulo, largura_texto), xy=(0, 1), xycoords="axes fraction",
                    xytext=(0, deslocamento), textcoords="offset points", ha="left", va="bottom",
                    fontsize=base, color="#4D4D4D")
        deslocamento += base * 1.35 * (1 + subtitulo.count("\n"))
    if titulo:
        ax.annotate(textwrap.fill(titulo, int(largura_texto / 1.2)), xy=(0, 1), xycoords="axes fraction",
                    xytext=(0, deslocamento), textcoords="offset points", ha="left", va="bottom",
                    fontsize=base * 1.2, color="black")

    if com_legenda:
        comuns = dict(title=titulo_legenda or None, frameon=False, fontsize=tamanho, title_fontsize=base,
                      handlelength=1.9, handletextpad=0.6, labelspacing=0.45, alignment="left")
        if opcoes.legenda == "Abaixo":
            fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center", ncols=colunas_abaixo,
                       columnspacing=1.2, **comuns)
        else:
            linhas = max(1, int((altura_pt - 60) / (tamanho * 1.6)))
            ax.legend(loc="center left", bbox_to_anchor=(1.015, 0.5), borderaxespad=0,
                      ncols=math.ceil(len(series) / linhas), columnspacing=1.0, **comuns)
    return fig


# =============================================================================
# Plano: da tabela filtrada aos gráficos
# =============================================================================

class _Campos(dict):
    def __missing__(self, chave):
        return "{" + chave + "}"


def _preencher(modelo: str, campos: dict) -> str:
    try:
        return str(modelo or "").format_map(_Campos(campos)).strip()
    except (ValueError, IndexError):
        return str(modelo or "").strip()


class Plano:
    """
    Reúne uma tabela (já filtrada) e as escolhas do usuário e sabe produzir
    cada gráfico: um por combinação das colunas de `painel` e por métrica.
    As cores, os marcadores e os nomes das séries são definidos uma vez para
    a tabela toda, para que a mesma classe tenha a mesma cor em todos os gráficos.
    """

    def __init__(self, df: pd.DataFrame, x: str, serie: Optional[str] = None, painel: Sequence[str] = (),
                 opcoes: Optional[OpcoesGrafico] = None, legenda: Optional[dict] = None, agregacao: str = "Média",
                 chaves_series: Optional[Sequence[str]] = None):
        if x not in df.columns:
            raise ValueError("Escolha a coluna do eixo X.")
        self.df, self.x, self.serie = df, x, (serie or None)
        self.painel = [c for c in painel if c not in (x, serie)]
        self.opcoes = opcoes or OpcoesGrafico()
        self.funcao = AGREGACOES.get(agregacao, "mean")
        legenda = legenda or {}

        # `chaves_series` fixa a lista completa de valores da série (a da tabela sem filtros),
        # para a cor de cada classe não mudar quando um filtro esconde outras classes.
        chaves = []
        if self.serie:
            chaves = list(chaves_series) if chaves_series is not None else \
                [texto_valor(v) for v in ordenar_valores(df[self.serie].dropna().unique())]
            if df[self.serie].isna().any() and SEM_VALOR not in chaves:
                chaves.append(SEM_VALOR)
        paleta = paleta_hue(len(chaves))
        self.chaves = chaves
        self.cores = {k: paleta[i] for i, k in enumerate(chaves)}
        self.marcadores = {k: (MARCADORES[i % len(MARCADORES)] if self.opcoes.marcadores_distintos else "o")
                           for i, k in enumerate(chaves)}
        self.rotulos = {k: k for k in chaves}
        for k in chaves:
            cor, nome = legenda.get(k, (None, None))
            if cor and self.opcoes.usar_cores_do_arquivo:
                self.cores[k] = cor
            if nome:
                self.rotulos[k] = nome

        self.x_numerico = pd.api.types.is_numeric_dtype(df[x]) and not pd.api.types.is_bool_dtype(df[x])
        self.x_inteiro = self.x_numerico and _inteiros(df[x])
        if self.painel:
            indices = df.groupby(self.painel, dropna=False, sort=False).indices
            self._grupos = {(k if isinstance(k, tuple) else (k,)): v for k, v in indices.items()}
            self.combos = sorted(self._grupos, key=lambda c: [chave_natural(texto_valor(v)) for v in c])
        else:
            self._grupos, self.combos = {(): None}, [()]

    # -- textos --------------------------------------------------------------
    def texto_combo(self, combo: tuple) -> str:
        return " | ".join(f"{nome_bonito(c)}: {texto_valor(v)}" for c, v in zip(self.painel, combo))

    def nome_arquivo(self, metrica: str, combo: tuple) -> Path:
        """Caminho relativo (sem extensão): <métrica>/<métrica>__<coluna>-<valor>_..."""
        base = nome_seguro(metrica)
        if not self.painel:
            return Path(base)
        sufixo = "_".join(f"{nome_seguro(c)}-{nome_seguro(texto_valor(v))}" for c, v in zip(self.painel, combo))
        return Path(base) / f"{base}__{sufixo}"

    # -- dados ---------------------------------------------------------------
    def bloco(self, combo: tuple) -> pd.DataFrame:
        posicoes = self._grupos[combo]
        return self.df if posicoes is None else self.df.iloc[posicoes]

    def dados(self, bloco: pd.DataFrame, metrica: str) -> pd.DataFrame:
        """Um ponto por (série, X): linhas repetidas são agregadas (média por padrão)."""
        chaves = ([self.serie] if self.serie else []) + [self.x]
        sub = bloco[chaves + [metrica]].dropna(subset=[metrica, self.x])
        if sub.empty:
            return sub
        return sub.groupby(chaves, dropna=False, sort=False)[metrica].agg(self.funcao).reset_index()

    # -- figura --------------------------------------------------------------
    def figura(self, dados: pd.DataFrame, metrica: str, combo: tuple = ()):
        categorias = None
        if self.x_numerico:
            dados = dados.sort_values(self.x, kind="stable")
            posicao_x = dados[self.x].to_numpy(dtype=float)
        else:
            categorias = sorted({texto_valor(v) for v in dados[self.x]}, key=chave_natural)
            indice = {c: i for i, c in enumerate(categorias)}
            posicao_x = np.array([indice[texto_valor(v)] for v in dados[self.x]], dtype=float)
        valores = dados[metrica].to_numpy(dtype=float)
        todos_x = np.unique(posicao_x)

        def com_lacunas(xs, ys):
            """Interrompe a linha onde a série pula valores de X que existem no gráfico
            (ex.: classe ausente em alguns anos), em vez de ligar pontos distantes."""
            if self.opcoes.ligar_lacunas or len(xs) < 2:
                return xs, ys
            saltos = np.where(np.diff(np.searchsorted(todos_x, xs)) > 1)[0] + 1
            return np.insert(xs, saltos, np.nan), np.insert(ys, saltos, np.nan)

        series = []
        if self.serie:
            chaves = dados[self.serie].map(texto_valor).to_numpy()
            for k in self.chaves:
                marca = chaves == k
                if marca.any():
                    ordem = np.argsort(posicao_x[marca], kind="stable")
                    xs, ys = com_lacunas(posicao_x[marca][ordem], valores[marca][ordem])
                    series.append({"rotulo": self.rotulos[k], "cor": self.cores[k], "marcador": self.marcadores[k],
                                   "x": xs, "y": ys})
        else:
            ordem = np.argsort(posicao_x, kind="stable")
            series.append({"rotulo": "", "cor": "#333333", "marcador": "o", "x": posicao_x[ordem], "y": valores[ordem]})

        campos = {"metrica": nome_metrica(metrica), "sigla": str(metrica), "painel": self.texto_combo(combo)}
        campos.update({str(c): texto_valor(v) for c, v in zip(self.painel, combo)})
        op = self.opcoes
        return desenhar_grafico(
            series, op, categorias_x=categorias, x_inteiro=self.x_inteiro,
            marcas_x=todos_x if categorias is None and len(todos_x) <= 12 else None,
            titulo=_preencher(op.titulo, campos), subtitulo=_preencher(op.subtitulo, campos),
            rotulo_x=_preencher(op.rotulo_x, campos) or nome_bonito(self.x),
            rotulo_y=_preencher(op.rotulo_y, campos) or nome_metrica(metrica),
            titulo_legenda=(_preencher(op.titulo_legenda, campos) or nome_bonito(self.serie)) if self.serie else "",
        )

    def salvar(self, metrica: str, combo: tuple, pasta_saida: str, formatos: Sequence[str] = ("png",),
               dpi: int = 300) -> List[Path]:
        """Gera e grava um gráfico; devolve os arquivos criados (lista vazia se não houver dados)."""
        dados = self.dados(self.bloco(combo), metrica)
        if dados.empty:
            return []
        fig = self.figura(dados, metrica, combo)
        base = Path(pasta_saida) / self.nome_arquivo(metrica, combo)
        base.parent.mkdir(parents=True, exist_ok=True)
        arquivos = []
        for formato in formatos:
            arquivo = base.parent / f"{base.name}.{formato.lower()}"
            fig.savefig(arquivo, dpi=dpi, facecolor="white")
            arquivos.append(arquivo)
        return arquivos

    def png(self, metrica: str, combo: tuple, largura_px: int, altura_px: int) -> Optional[bytes]:
        """PNG em memória, na resolução que cabe em largura_px x altura_px (pré-visualização)."""
        dados = self.dados(self.bloco(combo), metrica)
        if dados.empty:
            return None
        fig = self.figura(dados, metrica, combo)
        dpi = min(largura_px / (self.opcoes.largura_cm / 2.54), altura_px / (self.opcoes.altura_cm / 2.54))
        memoria = io.BytesIO()
        fig.savefig(memoria, format="png", dpi=max(40.0, min(220.0, dpi)), facecolor="white")
        return memoria.getvalue()


def gerar_graficos(plano: Plano, metricas: Sequence[str], pasta_saida: str, formatos: Sequence[str] = ("png",),
                   dpi: int = 300, log: Optional[Callable[[str], None]] = None,
                   progresso: Optional[Callable[[int, int, str], None]] = None,
                   checar_parada: Optional[Callable[[], None]] = None) -> Dict[str, object]:
    """Gera todos os gráficos do plano (cada combinação x cada métrica) e grava na pasta de saída."""
    if not metricas:
        raise ValueError("Marque ao menos uma métrica.")
    if not formatos:
        raise ValueError("Marque ao menos um formato de arquivo.")
    total = len(plano.combos) * len(metricas)
    gerados, sem_dados, arquivos = 0, 0, []
    feitos = 0
    for combo in plano.combos:
        for metrica in metricas:
            if checar_parada:
                checar_parada()
            criados = plano.salvar(metrica, combo, pasta_saida, formatos, dpi)
            feitos += 1
            if criados:
                gerados += 1
                arquivos.extend(criados)
                if log:
                    log(f"[{feitos}/{total}] {criados[0].relative_to(pasta_saida).with_suffix('')}")
            else:
                sem_dados += 1
            if progresso:
                progresso(feitos, total, f"Gráficos: {feitos}/{total}")
    return {"gerados": gerados, "sem_dados": sem_dados, "arquivos": arquivos, "total": total}


# =============================================================================
# Interface
# =============================================================================

import tkinter as tk  # noqa: E402
from tkinter import colorchooser, filedialog, messagebox, simpledialog, ttk  # noqa: E402


class ListaMarcavel(ttk.Frame):
    """Lista com caixas de seleção, barra de rolagem e botões Todos/Nenhum."""

    def __init__(self, parent, altura: int = 5, ao_mudar: Optional[Callable[[], None]] = None, largura: int = 200):
        super().__init__(parent)
        self._itens: list = []
        self._rotulos: List[str] = []
        self._marcados: set = set()
        self._ao_mudar = ao_mudar
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.tree = ttk.Treeview(self, columns=("texto",), show="", height=altura, selectmode="none")
        self.tree.column("texto", width=largura, stretch=True)
        self.tree.grid(row=0, column=0, sticky="nsew")
        barra = ttk.Scrollbar(self, orient="vertical", command=self.tree.yview)
        barra.grid(row=0, column=1, sticky="ns")
        self.tree.configure(yscrollcommand=barra.set)
        self.tree.bind("<Button-1>", self._clique)
        botoes = ttk.Frame(self)
        botoes.grid(row=0, column=2, sticky="n", padx=(4, 0))
        ttk.Button(botoes, text="Todos", width=8, command=self.marcar_todos).pack(pady=(0, 2))
        ttk.Button(botoes, text="Nenhum", width=8, command=self.desmarcar_todos).pack()

    def definir(self, itens: Sequence, rotulos: Optional[Sequence[str]] = None, marcados: Optional[Sequence] = None):
        self._itens = list(itens)
        self._rotulos = [str(r) for r in (rotulos if rotulos is not None else itens)]
        escolhidos = set(self._itens if marcados is None else marcados)
        self._marcados = {i for i, item in enumerate(self._itens) if item in escolhidos}
        self.tree.delete(*self.tree.get_children())
        for i in range(len(self._itens)):
            self.tree.insert("", "end", iid=str(i), values=(self._texto(i),))

    def _texto(self, i: int) -> str:
        return ("☑  " if i in self._marcados else "☐  ") + self._rotulos[i]

    def _redesenhar(self):
        for i in range(len(self._itens)):
            self.tree.item(str(i), values=(self._texto(i),))
        if self._ao_mudar:
            self._ao_mudar()

    def _clique(self, evento):
        linha = self.tree.identify_row(evento.y)
        if linha:
            self._marcados ^= {int(linha)}
            self._redesenhar()
        return "break"

    def itens(self) -> list:
        return list(self._itens)

    def marcados(self) -> list:
        return [self._itens[i] for i in sorted(self._marcados)]

    def marcar(self, itens: Sequence):
        escolhidos = set(itens)
        self._marcados = {i for i, item in enumerate(self._itens) if item in escolhidos}
        self._redesenhar()

    def marcar_todos(self):
        self._marcados = set(range(len(self._itens)))
        self._redesenhar()

    def desmarcar_todos(self):
        self._marcados = set()
        self._redesenhar()


class AbaGraficos:
    """Aba "10. Gráficos": tabela -> campos automáticos -> gráficos individuais e em lote."""

    MODO_INDICE = False
    TITULO = "Gráficos de linhas e pontos no padrão do ggplot2, a partir de qualquer tabela"

    def __init__(self, app, aba):
        self.app = app
        self.df: Optional[pd.DataFrame] = None
        self.estrutura = Estrutura()
        self.filtros: Dict[str, Optional[set]] = {}
        self.legenda: dict = {}
        self._valores_dim: Dict[str, list] = {}
        self._combos_previa: List[tuple] = []
        self._imagem_previa = None
        self._atualizando = False
        self.stop_event = threading.Event()
        self._montar(aba)

    # ------------------------------------------------------------------ layout
    def _montar(self, aba):
        app = self.app
        self.tabela_var, self.saida_var, self.legenda_var = tk.StringVar(), tk.StringVar(), tk.StringVar()
        self.x_var, self.serie_var = tk.StringVar(), tk.StringVar(value=NENHUMA)
        self.agregacao_var = tk.StringVar(value="Média")
        self.filtro_var = tk.StringVar()
        self.tema_var, self.posicao_legenda_var = tk.StringVar(value=TEMA_PADRAO), tk.StringVar(value=LEGENDAS[0])
        self.png_var, self.pdf_var, self.svg_var = tk.BooleanVar(value=True), tk.BooleanVar(), tk.BooleanVar()
        self.dpi_var = tk.StringVar(value="300")
        self.largura_var, self.altura_var, self.fonte_var = tk.StringVar(value="18"), tk.StringVar(value="10"), tk.StringVar(value="11")
        self.marcadores_var, self.virgula_var = tk.BooleanVar(value=True), tk.BooleanVar(value=False)
        self.lacunas_var = tk.BooleanVar(value=False)
        self.cores_arquivo_var = tk.BooleanVar(value=True)
        self.titulo_var, self.subtitulo_var = tk.StringVar(value="{painel}"), tk.StringVar()
        self.rotulo_x_var, self.rotulo_y_var = tk.StringVar(), tk.StringVar(value="{metrica}")
        self.titulo_legenda_var = tk.StringVar()
        self.previa_metrica_var, self.previa_combo_var = tk.StringVar(), tk.StringVar()

        self.frame = frame = ttk.Frame(aba, padding=10)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)
        ttk.Label(frame, text=self.TITULO, font=("Segoe UI", 12, "bold")).grid(row=0, column=0, columnspan=3, sticky="w", padx=6, pady=(0, 4))
        app._adicionar_linha_arquivo(frame, "Tabela de dados (CSV, TXT ou Parquet):", self.tabela_var,
                                     self._acao_selecionar_tabela, row=1)
        app._adicionar_linha_pasta(frame, "Pasta de saída dos gráficos:", self.saida_var,
                                   self._acao_selecionar_saida, row=2)

        corpo = ttk.Frame(frame)
        corpo.grid(row=3, column=0, columnspan=3, sticky="nsew", padx=6, pady=2)
        corpo.columnconfigure(1, weight=1)
        corpo.rowconfigure(0, weight=1)
        esquerda = ttk.Frame(corpo)
        esquerda.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
        esquerda.rowconfigure(1, weight=1)
        direita = ttk.Frame(corpo)
        direita.grid(row=0, column=1, sticky="nsew")
        direita.columnconfigure(0, weight=1)
        direita.rowconfigure(1, weight=1)

        # 1. Estrutura
        est = ttk.LabelFrame(esquerda, text="1. Estrutura do gráfico (preenchida ao ler a tabela)", padding=6)
        est.grid(row=0, column=0, sticky="ew")
        ttk.Label(est, text="Eixo X:").grid(row=0, column=0, sticky="e", padx=(0, 4))
        self.combo_x = ttk.Combobox(est, textvariable=self.x_var, state="readonly", width=15)
        self.combo_x.grid(row=0, column=1, sticky="w")
        ttk.Label(est, text="Linhas:").grid(row=0, column=2, sticky="e", padx=(10, 4))
        self.combo_serie = ttk.Combobox(est, textvariable=self.serie_var, state="readonly", width=15)
        self.combo_serie.grid(row=0, column=3, sticky="w")
        ttk.Label(est, text="Valores repetidos:").grid(row=1, column=0, columnspan=2, sticky="e", padx=(0, 4), pady=(4, 0))
        ttk.Combobox(est, textvariable=self.agregacao_var, state="readonly", width=10,
                     values=list(AGREGACOES)).grid(row=1, column=2, columnspan=2, sticky="w", pady=(4, 0))
        ttk.Label(est, text="Um gráfico para cada combinação de:").grid(row=2, column=0, columnspan=4, sticky="w", pady=(4, 0))
        self.lista_painel = ListaMarcavel(est, altura=3, ao_mudar=self._ao_mudar_estrutura, largura=230)
        self.lista_painel.grid(row=3, column=0, columnspan=4, sticky="ew")
        for combo in (self.combo_x, self.combo_serie):
            combo.bind("<<ComboboxSelected>>", lambda _e: self._ao_mudar_estrutura())

        # 2. Métricas
        met = ttk.LabelFrame(esquerda, text="2. Métricas (eixo Y) - um gráfico por métrica marcada", padding=6)
        met.grid(row=1, column=0, sticky="nsew", pady=(4, 0))
        met.columnconfigure(0, weight=1)
        met.rowconfigure(0, weight=1)
        self.lista_metricas = ListaMarcavel(met, altura=5, ao_mudar=self._ao_mudar_metricas, largura=230)
        self.lista_metricas.grid(row=0, column=0, sticky="nsew")

        # Classes: escolhe quais valores da coluna das linhas entram nos gráficos
        # (é o mesmo filtro da caixa 3 para essa coluna, só que sempre à vista).
        self.caixa_classes = ttk.LabelFrame(esquerda, text="Classes nos gráficos (valores da coluna das Linhas)", padding=6)
        self.caixa_classes.grid(row=2, column=0, sticky="ew", pady=(4, 0))
        self.caixa_classes.columnconfigure(0, weight=1)
        self.lista_classes = ListaMarcavel(self.caixa_classes, altura=4, ao_mudar=self._ao_mudar_classes, largura=230)
        self.lista_classes.grid(row=0, column=0, sticky="ew")

        # 3. Filtros
        fil = ttk.LabelFrame(esquerda, text="3. Filtros (valores de cada coluna)", padding=6)
        fil.grid(row=3, column=0, sticky="ew", pady=(4, 0))
        fil.columnconfigure(1, weight=1)
        ttk.Label(fil, text="Coluna:").grid(row=0, column=0, sticky="e", padx=(0, 4))
        self.combo_filtro = ttk.Combobox(fil, textvariable=self.filtro_var, state="readonly", width=30)
        self.combo_filtro.grid(row=0, column=1, sticky="w")
        self.combo_filtro.bind("<<ComboboxSelected>>", lambda _e: self._mostrar_filtro())
        self.lista_filtro = ListaMarcavel(fil, altura=4, ao_mudar=self._ao_mudar_filtro, largura=230)
        self.lista_filtro.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(4, 0))

        # 4. Aparência
        apa = ttk.LabelFrame(direita, text="4. Aparência", padding=6)
        apa.grid(row=0, column=0, sticky="ew")
        l0 = ttk.Frame(apa)
        l0.pack(fill="x")
        ttk.Label(l0, text="Tema:").pack(side="left")
        ttk.Combobox(l0, textvariable=self.tema_var, state="readonly", width=31, values=list(TEMAS)).pack(side="left", padx=(4, 10))
        ttk.Label(l0, text="Legenda:").pack(side="left")
        ttk.Combobox(l0, textvariable=self.posicao_legenda_var, state="readonly", width=11, values=list(LEGENDAS)).pack(side="left", padx=(4, 10))
        ttk.Label(l0, text="Formatos:").pack(side="left")
        for texto, var in (("PNG", self.png_var), ("PDF", self.pdf_var), ("SVG", self.svg_var)):
            ttk.Checkbutton(l0, text=texto, variable=var).pack(side="left", padx=2)
        ttk.Checkbutton(l0, text="Vírgula decimal", variable=self.virgula_var).pack(side="left", padx=(10, 2))
        l1 = ttk.Frame(apa)
        l1.pack(fill="x", pady=(4, 0))
        for texto, var, largura in (("Largura (cm):", self.largura_var, 5), ("Altura (cm):", self.altura_var, 5),
                                    ("Fonte (pt):", self.fonte_var, 4)):
            ttk.Label(l1, text=texto).pack(side="left")
            ttk.Entry(l1, textvariable=var, width=largura).pack(side="left", padx=(4, 10))
        ttk.Label(l1, text="DPI:").pack(side="left")
        ttk.Combobox(l1, textvariable=self.dpi_var, width=5, values=["150", "300", "600"]).pack(side="left", padx=(4, 10))
        ttk.Checkbutton(l1, text="Marcadores distintos", variable=self.marcadores_var).pack(side="left", padx=2)
        ttk.Checkbutton(l1, text="Ligar lacunas", variable=self.lacunas_var).pack(side="left", padx=2)
        l2 = ttk.Frame(apa)
        l2.pack(fill="x", pady=(4, 0))
        for texto, var, largura, cresce in (("Título:", self.titulo_var, 11, True), ("Subtítulo:", self.subtitulo_var, 9, True),
                                            ("Eixo Y:", self.rotulo_y_var, 11, True), ("Eixo X:", self.rotulo_x_var, 7, False),
                                            ("Legenda:", self.titulo_legenda_var, 7, False)):
            ttk.Label(l2, text=texto).pack(side="left")
            ttk.Entry(l2, textvariable=var, width=largura).pack(side="left", padx=(3, 8), fill="x", expand=cresce)
        l3 = ttk.Frame(apa)
        l3.pack(fill="x", pady=(4, 0))
        ttk.Label(l3, text="Cores e nomes das classes:").pack(side="left")
        ttk.Entry(l3, textvariable=self.legenda_var).pack(side="left", fill="x", expand=True, padx=4)
        ttk.Button(l3, text="Carregar .txt", command=self._acao_selecionar_legenda).pack(side="left")
        ttk.Button(l3, text="Escolher cores...", command=self.abrir_dialogo_cores).pack(side="left", padx=(4, 0))
        ttk.Checkbutton(l3, text="Usar estas cores", variable=self.cores_arquivo_var).pack(side="left", padx=(8, 0))
        ttk.Label(apa, foreground="#555555", text="Os textos aceitam {metrica}, {sigla} (nome da coluna) e {painel}; "
                  "Eixo X e Legenda vazios usam o nome da coluna.").pack(anchor="w", pady=(3, 0))

        # 5. Pré-visualização / gráfico individual
        pre = ttk.LabelFrame(direita, text="5. Gráfico individual (pré-visualização)", padding=6)
        pre.grid(row=1, column=0, sticky="nsew", pady=(4, 0))
        pre.columnconfigure(0, weight=1)
        pre.rowconfigure(1, weight=1)
        l4 = ttk.Frame(pre)
        l4.grid(row=0, column=0, sticky="ew")
        ttk.Label(l4, text="Métrica:").pack(side="left")
        self.combo_previa_metrica = ttk.Combobox(l4, textvariable=self.previa_metrica_var, state="readonly", width=16)
        self.combo_previa_metrica.pack(side="left", padx=(4, 8))
        ttk.Label(l4, text="Gráfico:").pack(side="left")
        self.combo_previa_combo = ttk.Combobox(l4, textvariable=self.previa_combo_var, state="readonly", width=22)
        self.combo_previa_combo.pack(side="left", padx=(4, 8))
        ttk.Button(l4, text="Pré-visualizar", command=self.pre_visualizar).pack(side="left", padx=2)
        ttk.Button(l4, text="Salvar este gráfico", command=self.salvar_individual).pack(side="left", padx=2)
        for combo in (self.combo_previa_metrica, self.combo_previa_combo):
            combo.bind("<<ComboboxSelected>>", lambda _e: self.pre_visualizar())
        self.rotulo_previa = ttk.Label(pre, anchor="center", foreground="#555555",
                                       text="Selecione uma tabela: os campos são criados automaticamente\n"
                                            "e o primeiro gráfico aparece aqui.")
        self.rotulo_previa.grid(row=1, column=0, sticky="nsew", pady=(4, 0))

        # 6. Índice de Estabilidade da Paisagem: caixa própria, separada das métricas
        # (na aba só do índice ela não aparece: lá a aba inteira é do índice).
        self.indice_vars = {coluna: tk.BooleanVar(value=(coluna == "IEP_u")) for coluna in NOMES_INDICE}
        self.indice_x_var = tk.StringVar(value="Raio")
        if not self.MODO_INDICE:
            ind = ttk.LabelFrame(frame, padding=4, text="6. Índice de Estabilidade da Paisagem (IEP) - por classe, calculado "
                                 "com CA, PLAND, NP, LPI e ENN_MN de todos os anos da tabela")
            ind.grid(row=4, column=0, columnspan=3, sticky="ew", padx=6, pady=(4, 0))
            for coluna, texto in (("IEP_u", "IEP"), ("V_u", "Variabilidade global"), ("V_comp", "Composição"),
                                  ("V_conf", "Configuração"), ("V_conn", "Conectividade")):
                ttk.Checkbutton(ind, text=texto, variable=self.indice_vars[coluna]).pack(side="left", padx=(2, 8))
            ttk.Label(ind, text="Eixo X:").pack(side="left", padx=(10, 4))
            ttk.Combobox(ind, textvariable=self.indice_x_var, state="readonly", width=8, values=["Raio", "Buffer"]).pack(side="left")
            self.botao_indice = ttk.Button(
                ind, text="Calcular o índice e gerar os gráficos do índice",
                command=lambda: app._executar_thread(self.botao_indice, self.botao_parar, self.stop_event, self._gerar_indice))
            self.botao_indice.pack(side="left", fill="x", expand=True, padx=(12, 2))

        # Rodapé: carregar, gerar tudo, parar, progresso e log
        rodape = ttk.Frame(frame)
        rodape.grid(row=5, column=0, columnspan=3, sticky="ew", padx=6, pady=(6, 2))
        rodape.columnconfigure(1, weight=1)
        self.botao_carregar = ttk.Button(
            rodape, text="Ler tabela e criar campos",
            command=lambda: app._executar_thread(self.botao_carregar, self.botao_parar, self.stop_event, self._ler_e_estruturar))
        self.botao_carregar.grid(row=0, column=0, padx=(0, 6))
        self.botao_gerar = ttk.Button(
            rodape, text="Gerar TODOS os gráficos selecionados e salvar na pasta de saída",
            command=lambda: app._executar_thread(self.botao_gerar, self.botao_parar, self.stop_event, self._gerar_todos))
        self.botao_gerar.grid(row=0, column=1, sticky="ew", padx=(0, 6))
        self.botao_parar = ttk.Button(rodape, text="Parar", state="disabled",
                                      command=lambda: app._solicitar_parada(self.stop_event, self.log))
        self.botao_parar.grid(row=0, column=2)
        self.rotulo_resumo = ttk.Label(frame, text="Processamento: aguardando a tabela de dados")
        self.rotulo_resumo.grid(row=6, column=0, columnspan=3, sticky="w", padx=6)
        self.barra = ttk.Progressbar(frame, mode="determinate", maximum=100)
        self.barra.grid(row=7, column=0, columnspan=3, sticky="ew", padx=6, pady=(2, 4))
        self.log = tk.Text(frame, height=3, wrap="word")
        self.log.grid(row=8, column=0, columnspan=3, sticky="nsew", padx=6, pady=(0, 4))
        barra_log = ttk.Scrollbar(frame, orient="vertical", command=self.log.yview)
        barra_log.grid(row=8, column=3, sticky="ns")
        self.log.configure(yscrollcommand=barra_log.set)
        frame.rowconfigure(3, weight=4)
        frame.rowconfigure(8, weight=1)

    def _log(self, mensagem: str):
        self.app._log(self.log, mensagem)

    # ------------------------------------------------------- seleção de arquivos
    def _acao_selecionar_tabela(self):
        caminho = filedialog.askopenfilename(
            title="Selecione a tabela de dados",
            filetypes=[("Tabelas", "*.csv *.txt *.tsv *.parquet"), ("Todos os arquivos", "*.*")])
        if not caminho:
            return
        self.tabela_var.set(caminho)
        if not self.saida_var.get():
            self.saida_var.set(str(Path(caminho).parent / "graficos"))
        self.botao_carregar.invoke()

    def _acao_selecionar_saida(self):
        pasta = self.app._selecionar_pasta("Selecione a pasta onde os gráficos serão salvos")
        if pasta:
            self.saida_var.set(pasta)

    def _acao_selecionar_legenda(self):
        caminho = filedialog.askopenfilename(
            title="Arquivo de cores e nomes das classes (valor, \"#cor\", nome)",
            filetypes=[("Texto", "*.txt *.csv"), ("Todos os arquivos", "*.*")])
        if caminho:
            self.legenda_var.set(caminho)
            self._carregar_legenda()
            self.pre_visualizar()

    def _carregar_legenda(self):
        caminho = self.legenda_var.get().strip()
        self.legenda = {}
        if caminho:
            try:
                self.legenda = ler_legenda(caminho)
                self._log(f"Cores e nomes das classes: {Path(caminho).name} ({len(self.legenda)} valores).")
            except OSError as exc:
                self._log(f"Não foi possível ler o arquivo de cores e nomes: {exc}")

    # ------------------------------------------------------------ leitura
    def _ler_e_estruturar(self):
        """Lê a tabela (em segundo plano) e cria os campos da aba a partir dela."""
        caminho = self.tabela_var.get().strip()
        if not caminho:
            raise ValueError("Selecione a tabela de dados.")
        self.log.after(0, lambda: self.app._limpar_log(self.log))
        self._log(f"Lendo a tabela: {caminho}")
        df, estrutura = self._estruturar(ler_tabela(caminho))
        legenda_auto = None
        if not self.legenda_var.get().strip():
            legenda_auto = next((str(Path(caminho).parent / n) for n in ARQUIVOS_LEGENDA
                                 if (Path(caminho).parent / n).is_file()), None)
        pronto, erros = threading.Event(), []

        def preencher():
            try:
                self._preencher_campos(df, estrutura, legenda_auto)
            except Exception as exc:  # noqa: BLE001
                erros.append(exc)
            finally:
                pronto.set()

        self.frame.after(0, preencher)
        pronto.wait(120)
        if erros:
            raise erros[0]

    def _estruturar(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, Estrutura]:
        return estruturar_tabela(df)

    def _preencher_campos(self, df: pd.DataFrame, est: Estrutura, legenda_auto: Optional[str]):
        self._atualizando = True
        try:
            self.df, self.estrutura = df, est
            self._valores_dim = {c: ordenar_valores(df[c].dropna().unique()) for c in est.dimensoes}
            filtraveis = [c for c in est.dimensoes if len(self._valores_dim[c]) >= 2]
            self.filtros = {c: None for c in filtraveis}
            for coluna, valores in est.filtros_iniciais.items():
                if coluna in self.filtros:
                    self.filtros[coluna] = set(valores)
            outras = [c for c in df.columns if c not in est.dimensoes]
            self.combo_x.configure(values=[c for c in est.dimensoes if c != est.serie] + outras)
            self.x_var.set(est.x or "")
            self.combo_serie.configure(values=[NENHUMA] + filtraveis)
            self.serie_var.set(est.serie or NENHUMA)
            self.lista_painel.definir(filtraveis, marcados=est.painel)
            self.lista_metricas.definir([], marcados=[])
            self.rotulo_x_var.set("")
            self.titulo_legenda_var.set("")
            if legenda_auto:
                self.legenda_var.set(legenda_auto)
            self._carregar_legenda()
        finally:
            self._atualizando = False
        self._log(f"Tabela lida: {len(df)} linhas x {len(df.columns)} colunas.")
        if est.derivadas:
            self._log(f"Campos criados a partir do nome do arquivo (coluna {est.origem_derivadas}): "
                      f"{', '.join(est.derivadas)}.")
        self._log(f"Estrutura sugerida -> eixo X: {est.x} | linhas: {est.serie or NENHUMA} | um gráfico por: "
                  f"{' + '.join(est.painel) or '(tabela inteira)'} | métricas numéricas: {len(est.metricas)} | "
                  f"colunas com filtro: {', '.join(str(c) for c in filtraveis) or '-'}")
        self._atualizar_dependentes(primeira_vez=True)
        self.pre_visualizar()

    # ------------------------------------------------------- estado das escolhas
    def coluna_x(self) -> str:
        return self.x_var.get()

    def coluna_serie(self) -> Optional[str]:
        serie = self.serie_var.get()
        return None if serie in ("", NENHUMA) else serie

    def colunas_painel(self) -> List[str]:
        return [c for c in self.lista_painel.marcados() if c not in (self.coluna_x(), self.coluna_serie())]

    def tabela_filtrada(self) -> pd.DataFrame:
        df, mascara = self.df, None
        for coluna, valores in self.filtros.items():
            if valores is not None and coluna in df.columns:
                atual = df[coluna].isin(list(valores))
                mascara = atual if mascara is None else (mascara & atual)
        return df if mascara is None else df[mascara]

    def _ao_mudar_estrutura(self):
        if not self._atualizando and self.df is not None:
            self._atualizar_dependentes()
            self.pre_visualizar()

    def _ao_mudar_metricas(self):
        if not self._atualizando and self.df is not None:
            self._atualizar_previa_e_resumo(self.tabela_filtrada())

    def _rotulos_filtro(self) -> List[str]:
        return [f"{c}  ({len(self._valores_dim[c]) if self.filtros[c] is None else len(self.filtros[c])}"
                f"/{len(self._valores_dim[c])})" for c in self.filtros]

    def _mostrar_filtro(self):
        indice = self.combo_filtro.current()
        if indice < 0:
            return
        coluna = list(self.filtros)[indice]
        valores = self._valores_dim[coluna]
        rotulos = [texto_valor(v) for v in valores]
        if coluna == self.coluna_serie():  # mostra também o nome da classe, se houver
            rotulos = [f"{r} - {self.legenda[r][1]}" if self.legenda.get(r, (None, None))[1] else r for r in rotulos]
        self._atualizando = True
        try:
            self.lista_filtro.definir(valores, rotulos, marcados=self.filtros[coluna])
        finally:
            self._atualizando = False

    def _ao_mudar_filtro(self):
        indice = self.combo_filtro.current()
        if self._atualizando or self.df is None or indice < 0:
            return
        coluna = list(self.filtros)[indice]
        marcados = self.lista_filtro.marcados()
        self.filtros[coluna] = None if len(marcados) == len(self._valores_dim[coluna]) else set(marcados)
        self.combo_filtro.configure(values=self._rotulos_filtro())
        self.combo_filtro.current(indice)
        self._atualizar_dependentes()
        self.pre_visualizar()

    def definir_filtro(self, coluna: str, valores: Optional[Sequence]):
        """Define o filtro de uma coluna por código (None = todos os valores)."""
        self.filtros[coluna] = None if valores is None else set(valores)
        self._atualizar_dependentes()

    def _atualizar_dependentes(self, primeira_vez: bool = False):
        """Refaz as listas que dependem das escolhas: métricas com dado, filtros e pré-visualização."""
        self._atualizando = True
        try:
            filtrada = self.tabela_filtrada()
            usadas = {self.coluna_x(), self.coluna_serie(), *self.colunas_painel()}
            metricas = [c for c in self.estrutura.metricas if c not in usadas and filtrada[c].notna().any()]
            marcadas = [c for c in self.lista_metricas.marcados() if c in metricas]
            if primeira_vez and metricas:
                marcadas = metricas[:1]
            self.lista_metricas.definir(metricas, [f"{c}  -  {nome_metrica(c)}" for c in metricas], marcados=marcadas)
            atual = max(0, self.combo_filtro.current())
            self.combo_filtro.configure(values=self._rotulos_filtro())
            if self.filtros:
                self.combo_filtro.current(min(atual, len(self.filtros) - 1))
        finally:
            self._atualizando = False
        self._mostrar_filtro()
        self._mostrar_classes()
        self._atualizar_previa_e_resumo(filtrada)

    def _mostrar_classes(self):
        """Caixa "Classes": valores da coluna das linhas, marcados conforme o filtro dessa coluna."""
        serie = self.coluna_serie()
        valores = self._valores_dim.get(serie, []) if serie else []
        rotulos = [texto_valor(v) for v in valores]
        rotulos = [f"{r} - {self.legenda[r][1]}" if self.legenda.get(r, (None, None))[1] else r for r in rotulos]
        self.caixa_classes.configure(
            text=f"Classes nos gráficos (coluna das Linhas: {serie})" if serie
            else "Classes nos gráficos (escolha a coluna das Linhas na caixa 1)")
        self._atualizando = True
        try:
            self.lista_classes.definir(valores, rotulos, marcados=self.filtros.get(serie) if serie else None)
        finally:
            self._atualizando = False

    def _ao_mudar_classes(self):
        serie = self.coluna_serie()
        if self._atualizando or self.df is None or not serie or serie not in self.filtros:
            return
        marcados = self.lista_classes.marcados()
        self.filtros[serie] = None if len(marcados) == len(self._valores_dim[serie]) else set(marcados)
        self._atualizar_dependentes()
        self.pre_visualizar()

    def classes_selecionadas(self) -> list:
        """Valores da coluna das linhas que entram nos gráficos."""
        return self.lista_classes.marcados()

    def _atualizar_previa_e_resumo(self, filtrada: pd.DataFrame):
        metricas = self.lista_metricas.marcados()
        painel = self.colunas_painel()
        if painel and len(filtrada):
            chaves = filtrada.groupby(painel, dropna=False, sort=False).indices.keys()
            combos = sorted(((k if isinstance(k, tuple) else (k,)) for k in chaves),
                            key=lambda c: [chave_natural(texto_valor(v)) for v in c])
        else:
            combos = [()] if len(filtrada) else []
        self._combos_previa = combos
        self.combo_previa_metrica.configure(values=metricas)
        if self.previa_metrica_var.get() not in metricas:
            self.previa_metrica_var.set(metricas[0] if metricas else "")
        textos = [" | ".join(texto_valor(v) for v in c) or "(tabela inteira)" for c in combos]
        self.combo_previa_combo.configure(values=textos[:3000])
        if self.previa_combo_var.get() not in textos:
            self.previa_combo_var.set(textos[0] if textos else "")
        self.rotulo_resumo.configure(
            text=f"Linhas após os filtros: {len(filtrada)} | gráficos por métrica: {len(combos)} | métricas marcadas: "
                 f"{len(metricas)} | gráficos a gerar: {len(combos) * len(metricas)}")

    # ------------------------------------------------------------- gráficos
    def opcoes(self) -> OpcoesGrafico:
        def numero(var, nome, minimo, maximo):
            try:
                valor = float(str(var.get()).replace(",", "."))
            except ValueError:
                raise ValueError(f"{nome}: informe um número.") from None
            if not minimo <= valor <= maximo:
                raise ValueError(f"{nome}: use um valor entre {minimo:g} e {maximo:g}.")
            return valor

        return OpcoesGrafico(
            tema=self.tema_var.get(), largura_cm=numero(self.largura_var, "Largura (cm)", 5, 60),
            altura_cm=numero(self.altura_var, "Altura (cm)", 4, 60), fonte_pt=numero(self.fonte_var, "Fonte (pt)", 5, 30),
            legenda=self.posicao_legenda_var.get(), marcadores_distintos=bool(self.marcadores_var.get()),
            ligar_lacunas=bool(self.lacunas_var.get()),
            virgula_decimal=bool(self.virgula_var.get()), usar_cores_do_arquivo=bool(self.cores_arquivo_var.get()),
            titulo=self.titulo_var.get(), subtitulo=self.subtitulo_var.get(), rotulo_x=self.rotulo_x_var.get(),
            rotulo_y=self.rotulo_y_var.get(), titulo_legenda=self.titulo_legenda_var.get())

    def formatos(self) -> List[str]:
        return [f for f, var in (("png", self.png_var), ("pdf", self.pdf_var), ("svg", self.svg_var)) if var.get()]

    def dpi(self) -> int:
        try:
            return int(min(1200, max(50, float(self.dpi_var.get()))))
        except ValueError:
            raise ValueError("DPI: informe um número (ex.: 300).") from None

    def plano(self) -> Plano:
        if self.df is None:
            raise ValueError("Selecione a tabela de dados primeiro.")
        return Plano(self.tabela_filtrada(), self.coluna_x(), self.coluna_serie(), self.colunas_painel(),
                     self.opcoes(), self.legenda, self.agregacao_var.get(), chaves_series=self.chaves_serie())

    # ------------------------------------------------- índice de estabilidade
    def tabela_do_indice(self) -> pd.DataFrame:
        """
        Calcula o IEP com a tabela carregada. O índice usa TODOS os buffers,
        classes e raios da tabela (a normalização depende do conjunto); só o
        filtro da coluna do eixo X (anos) entra no cálculo, definindo o período.
        """
        if self.df is None:
            raise ValueError("Selecione a tabela de dados primeiro.")
        df, x = self.df, self.coluna_x()
        if self.filtros.get(x) is not None:
            df = df[df[x].isin(list(self.filtros[x]))]
        return calcular_iep_de_tabela_de_metricas(df)

    def _gerar_indice(self):
        """Caixa 6: calcula o índice de estabilidade e gera os gráficos dele (em segundo plano)."""
        app = self.app
        self.log.after(0, lambda: app._limpar_log(self.log))
        pasta = self.saida_var.get().strip()
        if not pasta:
            raise ValueError("Selecione a pasta de saída dos gráficos.")
        variaveis = [coluna for coluna, var in self.indice_vars.items() if var.get()]
        if not variaveis:
            raise ValueError("Marque ao menos um resultado do índice (IEP, variabilidade...).")
        formatos, dpi = self.formatos() or ["png"], self.dpi()
        indice = self.tabela_do_indice()
        destino = Path(pasta) / "indice_estabilidade"
        destino.mkdir(parents=True, exist_ok=True)
        indice.to_csv(destino / "indice_estabilidade_IEP.csv", index=False, encoding="utf-8-sig")
        self._log(f"Índice calculado: {len(indice)} combinações de buffer x classe x raio | anos por combinação: "
                  f"{int(indice['n_anos'].min())} a {int(indice['n_anos'].max())} | IEP médio = {indice['IEP_u'].mean():.3f}")
        self._log(f"Tabela do índice: {destino / 'indice_estabilidade_IEP.csv'}")

        # Os filtros de buffer, raio e classe escolhem o que é desenhado (não mudam o cálculo).
        classe = next(c for c in indice.columns if c not in ("Buffer", "Raio", "n_anos")
                      and not str(c).startswith(("CV_", "NORM_", "V_", "IEP_")))
        for coluna in ("Buffer", "Raio", classe):
            if self.filtros.get(coluna) is not None and coluna in indice.columns:
                indice = indice[indice[coluna].isin(list(self.filtros[coluna]))]
        eixo = self.indice_x_var.get() if self.indice_x_var.get() in indice.columns else "Buffer"
        painel = [c for c in ("Buffer", "Raio") if c != eixo and c in indice.columns and indice[c].nunique() >= 1]
        chaves = self.chaves_serie() if classe == self.coluna_serie() else None
        plano = Plano(indice, eixo, classe, painel, self.opcoes(), self.legenda, "Média", chaves_series=chaves)
        total = len(plano.combos) * len(variaveis)
        self._log(f"Gerando {total} gráfico(s) do índice: eixo X = {eixo} | linhas = classes | um gráfico por "
                  f"{' + '.join(painel) or '(tabela inteira)'}")
        modulo = __import__(type(app).__module__)
        resultado = gerar_graficos(
            plano, variaveis, str(destino), formatos, dpi, log=self._log,
            progresso=lambda valor, maximo, msg: app._atualizar_barra_processamento(self.barra, self.rotulo_resumo, valor, maximo, msg),
            checar_parada=lambda: modulo.checar_parada(self.stop_event))
        resumo = f"Gráficos do índice gerados: {resultado['gerados']} de {resultado['total']} | sem dados: {resultado['sem_dados']}"
        self._log(resumo)
        self.log.after(0, lambda: messagebox.showinfo("Concluído", resumo + f"\n\nPasta: {destino}"))

    # ------------------------------------------------------ dicionário de cores
    def chaves_serie(self) -> List[str]:
        """Todos os valores (em texto) da coluna das linhas, na tabela sem filtros."""
        serie = self.coluna_serie()
        return [texto_valor(v) for v in self._valores_dim.get(serie, [])] if serie else []

    def cores_efetivas(self, sempre_dicionario: bool = False) -> Dict[str, str]:
        """Cor de cada valor da série: a do dicionário, se houver; senão, a da paleta do ggplot2."""
        chaves = self.chaves_serie()
        cores = dict(zip(chaves, paleta_hue(len(chaves))))
        if sempre_dicionario or self.cores_arquivo_var.get():
            for chave in chaves:
                cor = self.legenda.get(chave, (None, None))[0]
                if cor:
                    cores[chave] = cor
        return cores

    def definir_cor(self, chave: str, cor: str):
        self.legenda[chave] = (cor.upper(), self.legenda.get(chave, (None, None))[1])
        self.cores_arquivo_var.set(True)

    def abrir_dialogo_cores(self):
        if self.df is None or not self.coluna_serie():
            messagebox.showinfo("Cores das classes", "Leia uma tabela e escolha a coluna das linhas primeiro.")
            return None
        anterior = getattr(self, "dialogo_cores", None)
        if anterior is not None and anterior.winfo_exists():
            anterior.destroy()
        self.dialogo_cores = DialogoCores(self)
        return self.dialogo_cores

    def _combo_da_previa(self) -> tuple:
        indice = self.combo_previa_combo.current()
        return self._combos_previa[indice] if 0 <= indice < len(self._combos_previa) else ()

    def pre_visualizar(self):
        """Desenha na aba o gráfico escolhido em "Gráfico individual"."""
        self._imagem_previa = None
        metrica = self.previa_metrica_var.get()
        if self.df is None or not metrica or not self._combos_previa:
            self.rotulo_previa.configure(image="", text="Marque ao menos uma métrica para ver o gráfico.")
            return
        try:
            largura = max(320, self.rotulo_previa.winfo_width() - 8)
            altura = max(200, self.rotulo_previa.winfo_height() - 8)
            dados = self.plano().png(metrica, self._combo_da_previa(), largura, altura)
        except Exception as exc:  # noqa: BLE001
            self.rotulo_previa.configure(image="", text=f"Não foi possível desenhar o gráfico:\n{exc}")
            return
        if dados is None:
            self.rotulo_previa.configure(image="", text="Sem dados para esta métrica neste gráfico.")
            return
        self._imagem_previa = tk.PhotoImage(data=base64.b64encode(dados))
        self.rotulo_previa.configure(image=self._imagem_previa, text="")

    def salvar_individual(self):
        """Grava só o gráfico da pré-visualização, nos formatos marcados."""
        try:
            metrica, pasta = self.previa_metrica_var.get(), self.saida_var.get().strip()
            if not metrica:
                raise ValueError("Marque ao menos uma métrica.")
            if not pasta:
                raise ValueError("Selecione a pasta de saída dos gráficos.")
            arquivos = self.plano().salvar(metrica, self._combo_da_previa(), pasta, self.formatos() or ["png"], self.dpi())
        except Exception as exc:  # noqa: BLE001
            messagebox.showerror("Gráfico individual", str(exc))
            return
        if not arquivos:
            messagebox.showwarning("Gráfico individual", "Sem dados para esta métrica neste gráfico.")
            return
        self._log("Gráfico salvo: " + "; ".join(str(a) for a in arquivos))
        messagebox.showinfo("Gráfico individual", "Gráfico salvo em:\n" + "\n".join(str(a) for a in arquivos))

    def _gerar_todos(self):
        """Botão principal: gera todos os gráficos selecionados (em segundo plano)."""
        app = self.app
        self.log.after(0, lambda: app._limpar_log(self.log))
        pasta = self.saida_var.get().strip()
        if not pasta:
            raise ValueError("Selecione a pasta de saída dos gráficos.")
        metricas, formatos, dpi = self.lista_metricas.marcados(), self.formatos(), self.dpi()
        plano = self.plano()
        total = len(plano.combos) * len(metricas)
        self._log(f"Gerando {total} gráfico(s): {len(metricas)} métrica(s) x {len(plano.combos)} combinação(ões) de "
                  f"{' + '.join(plano.painel) or '(tabela inteira)'} | formatos: {', '.join(formatos)} | {dpi} dpi")
        self._log(f"Pasta de saída: {pasta}")
        modulo = __import__(type(app).__module__)
        resultado = gerar_graficos(
            plano, metricas, pasta, formatos, dpi, log=self._log,
            progresso=lambda valor, maximo, msg: app._atualizar_barra_processamento(self.barra, self.rotulo_resumo, valor, maximo, msg),
            checar_parada=lambda: modulo.checar_parada(self.stop_event))
        resumo = (f"Gráficos gerados: {resultado['gerados']} de {resultado['total']} | sem dados: "
                  f"{resultado['sem_dados']} | arquivos: {len(resultado['arquivos'])}")
        self._log("=" * 70)
        self._log(resumo)
        self.log.after(0, lambda: messagebox.showinfo("Concluído", resumo + f"\n\nPasta: {pasta}"))


def tabela_do_iep(df: pd.DataFrame) -> Tuple[pd.DataFrame, bool]:
    """A tabela do índice: a própria, se já tiver o IEP calculado (coluna IEP_u), ou o
    IEP calculado a partir de uma tabela geral de métricas. Devolve (tabela, calculou)."""
    if "IEP_u" in df.columns:
        return df, False
    try:
        return calcular_iep_de_tabela_de_metricas(df), True
    except Exception as exc:  # noqa: BLE001
        raise ValueError("A tabela não tem o índice calculado (coluna IEP_u) e não foi possível calculá-lo a partir "
                         f"dela: {exc}\nUse a tabela do índice ou a tabela geral de métricas (nível de classe, com "
                         "CA, PLAND, NP, LPI e ENN_MN de vários anos).") from None


def estruturar_tabela_iep(df: pd.DataFrame) -> Tuple[pd.DataFrame, Estrutura]:
    """Estrutura sugerida para a tabela do índice: eixo X = Raio, linhas = classes,
    um gráfico por Buffer e, no eixo Y, os resultados do índice (IEP primeiro)."""
    df, est = estruturar_tabela(df)
    minusculas = {str(c).lower(): c for c in df.columns}
    classe = next((minusculas[n] for n in NOMES_SERIE if n in minusculas), None)
    for coluna in ("Buffer", "Raio", classe):  # Raio e classe são categorias, mesmo sendo números
        if coluna and coluna in df.columns and coluna not in est.dimensoes:
            est.dimensoes.append(coluna)
        if coluna in est.metricas:
            est.metricas.remove(coluna)
    resultados = [c for c in NOMES_INDICE if c in df.columns]
    extras = [c for c in df.columns if str(c).startswith(("CV_", "NORM_CV_")) and c not in resultados]
    est.metricas = resultados + extras + [c for c in est.metricas if c not in resultados and c not in extras]
    est.x = "Raio" if "Raio" in df.columns else ("Buffer" if "Buffer" in df.columns else est.x)
    est.serie = classe
    est.painel = [c for c in ("Buffer", "Raio") if c in df.columns and c != est.x]
    est.filtros_iniciais = {}
    return df, est


class AbaGraficosIndice(AbaGraficos):
    """Aba "11. Gráficos do IEP": os mesmos campos da aba de gráficos, só com o índice.

    Aceita a planilha do índice já calculado (indice_estabilidade_IEP*.csv) ou a tabela
    geral de métricas, da qual o índice é calculado ao ler (e gravado na pasta de saída).
    """

    MODO_INDICE = True
    TITULO = ("Gráficos do Índice de Estabilidade da Paisagem (IEP) - planilha do índice "
              "ou tabela geral de métricas")

    def __init__(self, app, aba):
        self.iep_calculado: Optional[pd.DataFrame] = None
        super().__init__(app, aba)
        self.rotulo_previa.configure(
            text="Selecione a planilha do índice (indice_estabilidade_IEP...csv)\nou a tabela geral de métricas: "
                 "o índice é calculado e o primeiro gráfico aparece aqui.")

    def _estruturar(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, Estrutura]:
        tabela, calculou = tabela_do_iep(df)
        self.iep_calculado = tabela if calculou else None
        if calculou:
            self._log(f"A tabela é de métricas: índice calculado com todos os anos ({len(tabela)} combinações de "
                      f"buffer x classe x raio; IEP médio = {tabela['IEP_u'].mean():.3f}). Ele será gravado na pasta "
                      "de saída ao gerar os gráficos.")
        else:
            self._log("A tabela já traz o índice calculado (coluna IEP_u): usada diretamente.")
        return estruturar_tabela_iep(tabela)

    def _gerar_todos(self):
        super()._gerar_todos()
        if self.iep_calculado is not None and self.saida_var.get().strip():
            destino = Path(self.saida_var.get().strip()) / "indice_estabilidade_IEP.csv"
            self.iep_calculado.to_csv(destino, index=False, encoding="utf-8-sig")
            self._log(f"Tabela do índice gravada: {destino}")


class DialogoCores(tk.Toplevel):
    """
    Dicionário de cores e nomes das classes: uma linha por valor da coluna das
    linhas, com a cor atual. A cor é escolhida no seletor de cores do sistema
    (ou preenchida para todas a partir de uma paleta) e o dicionário pode ser
    salvo em .txt e carregado depois.
    """

    def __init__(self, aba: AbaGraficos):
        super().__init__(aba.frame)
        self.aba = aba
        self.title("Cores e nomes das classes (dicionário)")
        self.geometry("600x560")
        self.transient(aba.frame.winfo_toplevel())
        self._chaves: List[str] = []
        self._amostras: list = []
        self.paleta_var = tk.StringVar(value=list(PALETAS)[0])

        ttk.Label(self, wraplength=570, justify="left", text=(
            "Cada valor da coluna das linhas tem uma cor. Dê duplo clique em uma linha (ou use \"Escolher cor\") "
            "para abrir o seletor de cores; \"Editar nome\" muda o texto que aparece na legenda. O dicionário pode "
            "ser salvo em .txt (valor, \"#cor\", nome) e carregado em outra sessão.")).pack(anchor="w", padx=10, pady=(10, 6))
        quadro = ttk.Frame(self)
        quadro.pack(fill="both", expand=True, padx=10)
        self.tree = ttk.Treeview(quadro, columns=("valor", "cor", "nome"), show="tree headings",
                                 selectmode="browse", height=14)
        self.tree.heading("#0", text="Cor")
        self.tree.column("#0", width=64, stretch=False, anchor="center")
        for coluna, texto, largura in (("valor", "Valor", 80), ("cor", "Código da cor", 100), ("nome", "Nome na legenda", 300)):
            self.tree.heading(coluna, text=texto)
            self.tree.column(coluna, width=largura, anchor="w")
        self.tree.pack(side="left", fill="both", expand=True)
        barra = ttk.Scrollbar(quadro, orient="vertical", command=self.tree.yview)
        barra.pack(side="left", fill="y")
        self.tree.configure(yscrollcommand=barra.set)
        self.tree.bind("<Double-1>", lambda _e: self.escolher_cor())

        l1 = ttk.Frame(self)
        l1.pack(fill="x", padx=10, pady=(8, 0))
        ttk.Button(l1, text="Escolher cor...", command=self.escolher_cor).pack(side="left")
        ttk.Button(l1, text="Editar nome...", command=self.editar_nome).pack(side="left", padx=6)
        ttk.Button(l1, text="Carregar .txt", command=self.carregar).pack(side="left", padx=(20, 6))
        ttk.Button(l1, text="Salvar .txt", command=self.salvar).pack(side="left")
        ttk.Button(l1, text="Fechar", command=self.destroy).pack(side="right")
        l2 = ttk.Frame(self)
        l2.pack(fill="x", padx=10, pady=(8, 10))
        ttk.Label(l2, text="Preencher todas com a paleta:").pack(side="left")
        combo = ttk.Combobox(l2, textvariable=self.paleta_var, state="readonly", width=22, values=list(PALETAS))
        combo.pack(side="left", padx=6)
        combo.bind("<<ComboboxSelected>>", lambda _e: self._desenhar_faixa())
        ttk.Button(l2, text="Aplicar", command=self.aplicar_paleta).pack(side="left")
        self.faixa = tk.Canvas(l2, width=160, height=16, highlightthickness=1, highlightbackground="#999999")
        self.faixa.pack(side="left", padx=(10, 0))
        self._desenhar_faixa()
        self.preencher()

    def _desenhar_faixa(self):
        """Mostra o contínuo de cores da paleta escolhida."""
        self.faixa.delete("all")
        cores = cores_da_paleta(self.paleta_var.get(), 80)
        for i, cor in enumerate(cores):
            self.faixa.create_rectangle(i * 2, 0, i * 2 + 2, 16, fill=cor, outline="")

    def preencher(self, selecionar: Optional[str] = None):
        self.tree.delete(*self.tree.get_children())
        self._chaves, self._amostras = self.aba.chaves_serie(), []
        cores = self.aba.cores_efetivas(sempre_dicionario=True)
        for i, chave in enumerate(self._chaves):
            amostra = tk.PhotoImage(master=self, width=30, height=14)
            amostra.put(cores[chave], to=(0, 0, 30, 14))
            self._amostras.append(amostra)
            nome = self.aba.legenda.get(chave, (None, None))[1] or ""
            self.tree.insert("", "end", iid=str(i), image=amostra, values=(chave, cores[chave], nome))
        if selecionar in self._chaves:
            self.tree.selection_set(str(self._chaves.index(selecionar)))

    def selecionada(self) -> Optional[str]:
        selecao = self.tree.selection()
        if not selecao:
            messagebox.showinfo("Cores das classes", "Selecione uma linha da lista.", parent=self)
            return None
        return self._chaves[int(selecao[0])]

    def _aplicar(self, selecionar: Optional[str] = None):
        self.aba.cores_arquivo_var.set(True)
        self.preencher(selecionar)
        self.aba.pre_visualizar()

    def escolher_cor(self):
        chave = self.selecionada()
        if chave is None:
            return
        atual = self.aba.cores_efetivas(sempre_dicionario=True)[chave]
        cor = colorchooser.askcolor(color=atual, parent=self, title=f"Cor do valor {chave}")[1]
        if cor:
            self.aba.definir_cor(chave, cor)
            self._aplicar(chave)

    def editar_nome(self):
        chave = self.selecionada()
        if chave is None:
            return
        cor, nome = self.aba.legenda.get(chave, (None, None))
        novo = simpledialog.askstring("Nome na legenda", f"Nome do valor {chave}:", initialvalue=nome or "", parent=self)
        if novo is not None:
            self.aba.legenda[chave] = (cor, novo.strip() or None)
            self._aplicar(chave)

    def aplicar_paleta(self):
        cores = cores_da_paleta(self.paleta_var.get(), len(self._chaves))
        for chave, cor in zip(self._chaves, cores):
            self.aba.definir_cor(chave, cor)
        self._aplicar()

    def dicionario_completo(self) -> dict:
        """Todos os valores da tabela com a cor em uso, mais as demais entradas já carregadas."""
        cores = self.aba.cores_efetivas(sempre_dicionario=True)
        completo = {k: (cores[k], self.aba.legenda.get(k, (None, None))[1]) for k in self._chaves}
        completo.update({k: v for k, v in self.aba.legenda.items() if k not in completo})
        return completo

    def salvar(self):
        caminho = filedialog.asksaveasfilename(parent=self, title="Salvar dicionário de cores", defaultextension=".txt",
                                               initialfile="cores.txt", filetypes=[("Texto", "*.txt")])
        if not caminho:
            return
        try:
            salvar_legenda(caminho, self.dicionario_completo())
        except OSError as exc:
            messagebox.showerror("Erro ao salvar", str(exc), parent=self)
            return
        self.aba.legenda_var.set(caminho)
        self.aba._carregar_legenda()
        self._aplicar()
        messagebox.showinfo("Salvo", f"Dicionário de cores salvo em:\n{caminho}", parent=self)

    def carregar(self):
        caminho = filedialog.askopenfilename(parent=self, title="Carregar dicionário de cores",
                                             filetypes=[("Texto", "*.txt *.csv"), ("Todos os arquivos", "*.*")])
        if not caminho:
            return
        self.aba.legenda_var.set(caminho)
        self.aba._carregar_legenda()
        self._aplicar()
