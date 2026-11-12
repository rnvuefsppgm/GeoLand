#!/usr/bin/env python3
"""
main.py
=======

Ponto de entrada do aplicativo GeoLand (Grade, Buffers, Recorte, Metricas
de paisagem, Reprojecao, Fundido e Consolidar).

Uso:

    python main.py

O arquivo original enviado pelo usuario
(app_grade_buffer_shapefile_v18_pylandstats_R.py) concentra toda a logica
(so o titulo da janela foi alterado): ele ja e uma aplicacao Tkinter completa com seu
proprio bloco "if __name__ == '__main__':" no final. Este main.py apenas
executa esse modulo como se ele tivesse sido chamado diretamente (via
runpy.run_module), o que faz o bloco original rodar normalmente - tanto
com "python main.py" quanto dentro do executavel gerado pelo PyInstaller
(por isso o modulo precisa estar listado em hiddenimports no GeoLand.spec,
ja que main.py nao faz um "import" estatico dele).
"""

from __future__ import annotations

import multiprocessing
import os
import runpy
import sys
import traceback
from pathlib import Path

MODULE_NAME = "app_grade_buffer_shapefile_v18_pylandstats_R"


def configurar_r_embutido() -> str | None:
    """Aponta o rpy2 para o R portatil distribuido junto com o aplicativo.

    Procura a pasta do R dentro do pacote do PyInstaller (_internal\\R), ao
    lado do executavel (R) ou, rodando pelo codigo-fonte, em R_portable. Se
    encontrar, define R_HOME/PATH e isola as bibliotecas de pacotes, de modo
    que o motor R (landscapemetrics) funcione sem R instalado na maquina e
    sem interferencia de um R do usuario. Se nao encontrar, nao altera nada
    e o rpy2 usa o R do sistema, se houver. Retorna o R_HOME usado ou None.
    """
    base = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
    candidatos = [
        base / "R",
        Path(sys.executable).resolve().parent / "R",
        base / "R_portable",
    ]
    for r_home in candidatos:
        if (r_home / "bin" / "x64" / "R.dll").is_file():
            break
    else:
        return None

    bin_dir = str(r_home / "bin" / "x64")
    os.environ["R_HOME"] = str(r_home)
    partes = os.environ.get("PATH", "").split(os.pathsep)
    if bin_dir not in partes:
        os.environ["PATH"] = os.pathsep.join([bin_dir] + partes)
    # Usa somente a biblioteca embutida e ignora .Rprofile/.Renviron do usuario.
    os.environ["R_LIBS_USER"] = str(r_home / "library")
    os.environ.pop("R_LIBS", None)
    os.environ.pop("R_LIBS_SITE", None)
    os.environ["R_PROFILE_USER"] = str(r_home / "etc" / "sem_Rprofile_usuario")
    os.environ["R_ENVIRON_USER"] = str(r_home / "etc" / "sem_Renviron_usuario")
    try:
        os.add_dll_directory(bin_dir)
    except (AttributeError, OSError):
        pass
    return str(r_home)


def autoteste_motores(arquivo_saida: str) -> int:
    """Diagnostico sem interface: GeoLand.exe --teste-motores relatorio.txt

    Calcula metricas de um raster minimo com os dois motores (PyLandStats e
    landscapemetrics do R embutido) e grava o resultado em `arquivo_saida`.
    """
    import importlib
    import tempfile

    linhas = ["R_HOME: %s" % os.environ.get("R_HOME")]
    ok = False
    try:
        import numpy as np
        import rasterio
        from rasterio.transform import from_origin

        app = importlib.import_module(MODULE_NAME)
        with tempfile.TemporaryDirectory() as tmp:
            tif = os.path.join(tmp, "teste_2020.tif")
            dados = np.kron(np.arange(16).reshape(4, 4) % 3 + 1, np.ones((10, 10))).astype("uint8")
            with rasterio.open(
                tif, "w", driver="GTiff", height=40, width=40, count=1, dtype="uint8",
                crs="EPSG:31984", transform=from_origin(500000, 9000000, 30, 30), nodata=0,
            ) as dst:
                dst.write(dados, 1)

            import pylandstats as pls

            df_py = pls.Landscape(tif).compute_class_metrics_df(metrics=["total_area", "number_of_patches"])
            linhas.append("PyLandStats: %d classe(s)" % len(df_py))

            funcs = {"class": ["lsm_c_ca", "lsm_c_np"], "landscape": ["lsm_l_np"]}
            _, df_r, mensagens = app._lsm_r_compute((tif, list(funcs), funcs, "8", 30.0))
            linhas += mensagens
            ro = app._R_STATE["ro"]
            linhas.append("R: %s" % ro.r("R.version.string")[0])
            linhas.append("Bibliotecas R: %s" % "; ".join(ro.r(".libPaths()")))
            ok = len(df_py) > 0 and df_r is not None and len(df_r) > 0
    except Exception:  # noqa: BLE001
        linhas.append(traceback.format_exc())
    linhas.append("RESULTADO: %s" % ("OK" if ok else "FALHA"))
    Path(arquivo_saida).write_text("\n".join(linhas), encoding="utf-8")
    return 0 if ok else 1


def main() -> int:
    # Precisa vir antes de qualquer import do rpy2 e antes de criar processos
    # filhos (que herdam estas variaveis de ambiente).
    configurar_r_embutido()
    if len(sys.argv) == 3 and sys.argv[1] == "--teste-motores":
        return autoteste_motores(sys.argv[2])
    # Necessario para o ProcessPoolExecutor (modo spawn) funcionar dentro
    # do executavel do PyInstaller; o proprio modulo original tambem chama
    # freeze_support() no seu bloco __main__, mas chamar aqui de novo e
    # inofensivo e garante o comportamento correto mesmo se o modulo for
    # atualizado no futuro.
    multiprocessing.freeze_support()
    try:
        runpy.run_module(MODULE_NAME, run_name="__main__")
        return 0
    except Exception:  # noqa: BLE001
        print("Erro fatal ao iniciar o GeoLand:", file=sys.stderr)
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
