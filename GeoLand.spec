# -*- mode: python ; coding: utf-8 -*-
"""
GeoLand.spec
=============

Arquivo de especificacao do PyInstaller para gerar o executavel Windows do
GeoLand (Grade, Buffers, Recorte, Metricas de paisagem, Reprojecao,
Fundido e Consolidar).

Uso:

    pyinstaller GeoLand.spec --noconfirm

Trata explicitamente os recursos de dados do GDAL/PROJ (rasterio e
pyogrio/geopandas), pois esses arquivos nao sao detectados automaticamente
pela analise estatica do PyInstaller, alem dos dados/binarios de
pylandstats, matplotlib, rpy2, duckdb e pyarrow. Requer o pacote
`pyinstaller-hooks-contrib` (instalado automaticamente como dependencia do
PyInstaller >= 6).

O motor R (rpy2 + landscapemetrics) e opcional: se o pacote rpy2 nao
estiver instalado no venv usado para o build, os blocos "try/except"
abaixo simplesmente pulam a coleta de seus dados e o .exe final funciona
normalmente com o motor PyLandStats (a aba de metricas apenas nao oferece
o motor R).
"""

import os

from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs, collect_submodules

block_cipher = None

# ---------------------------------------------------------------------------
# Dados e binarios exigidos em tempo de execucao
# ---------------------------------------------------------------------------

datas = []
binaries = []

# rasterio (GDAL embutido no wheel)
datas += collect_data_files("rasterio", include_py_files=False)
binaries += collect_dynamic_libs("rasterio")
try:
    import rasterio
    rasterio_gdal_data = os.path.join(os.path.dirname(rasterio.__file__), "gdal_data")
    if os.path.isdir(rasterio_gdal_data):
        datas.append((rasterio_gdal_data, "rasterio/gdal_data"))
    rasterio_proj_data = os.path.join(os.path.dirname(rasterio.__file__), "proj_data")
    if os.path.isdir(rasterio_proj_data):
        datas.append((rasterio_proj_data, "rasterio/proj_data"))
except Exception:
    pass

# pyproj (dados PROJ)
datas += collect_data_files("pyproj", include_py_files=False)
try:
    import pyproj
    proj_data_dir = pyproj.datadir.get_data_dir()
    if os.path.isdir(proj_data_dir):
        datas.append((proj_data_dir, "pyproj/proj_data"))
except Exception:
    pass

# pyogrio (GDAL embutido separadamente - usado pelo geopandas para ler/gravar .shp)
datas += collect_data_files("pyogrio", include_py_files=False)
binaries += collect_dynamic_libs("pyogrio")

# geopandas / shapely
datas += collect_data_files("geopandas", include_py_files=False)
binaries += collect_dynamic_libs("shapely")

# pylandstats (motor de metricas de paisagem em Python)
datas += collect_data_files("pylandstats", include_py_files=False)

# matplotlib (dependencia do pylandstats para plots internos)
datas += collect_data_files("matplotlib", include_py_files=False, subdir="mpl-data")

# duckdb / pyarrow (consolidacao de tabelas CSV/Parquet)
binaries += collect_dynamic_libs("duckdb")
datas += collect_data_files("pyarrow", include_py_files=False)
binaries += collect_dynamic_libs("pyarrow")

# rpy2 (motor R opcional - so entra no build se estiver instalado no venv)
try:
    import rpy2  # noqa: F401
    datas += collect_data_files("rpy2", include_py_files=False)
    binaries += collect_dynamic_libs("rpy2")
    RPY2_AVAILABLE = True
except Exception:
    RPY2_AVAILABLE = False

hidden_imports = [
    # main.py importa este modulo pelo nome (importlib, sem "import"
    # estatico), entao ele precisa ser listado explicitamente para o
    # PyInstaller incluir o bytecode dele no pacote final.
    "app_grade_buffer_shapefile_v18_pylandstats_R",
    # Aba 10 (graficos): modulo proprio, nomes das metricas e os formatos
    # de saida do matplotlib (PNG, PDF e SVG), carregados sob demanda.
    "geoland_graficos",
    "geoland_indice",
    "nomes_metricas_lsm",
    "matplotlib.backends.backend_agg",
    "matplotlib.backends.backend_pdf",
    "matplotlib.backends.backend_svg",
]
hidden_imports += collect_submodules("rasterio")
hidden_imports += collect_submodules("pyogrio")
hidden_imports += collect_submodules("shapely")
hidden_imports += collect_submodules("skimage")
hidden_imports += collect_submodules("scipy.stats")
hidden_imports += collect_submodules("pylandstats")
hidden_imports += [
    "rasterio.sample",
    "rasterio._io",
    "rasterio.vrt",
    "rasterio._features",
    "geopandas.io.arrow",
    "pandas._libs.tslibs.base",
    "duckdb",
    "pyarrow.parquet",
    "pyarrow.compute",
]

if RPY2_AVAILABLE:
    hidden_imports += collect_submodules("rpy2")
    # Modulo cffi (modo ABI) gerado na instalacao do rpy2: fica na raiz do
    # site-packages e e importado dinamicamente, entao o PyInstaller nao o ve.
    hidden_imports += ["_rinterface_cffi_abi", "_cffi_backend", "tzlocal", "jinja2"]

# R portatil (R + landscapemetrics + terra) copiado para _internal\R, para o
# motor R funcionar sem R instalado na maquina de destino. O main.py aponta
# R_HOME para essa pasta antes de carregar o rpy2.
R_PORTABLE = os.path.join(os.path.abspath("."), "R_portable")
r_tree = Tree(R_PORTABLE, prefix="R") if os.path.isdir(R_PORTABLE) else []

a = Analysis(
    ["main.py"],
    pathex=[os.path.abspath(".")],
    binaries=binaries,
    datas=datas,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="GeoLand_V1",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    icon=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    r_tree,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="GeoLand_V1",
)
