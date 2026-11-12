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
(so o titulo da janela foi alterado): ele ja e uma aplicacao Tkinter completa.
Este main.py IMPORTA esse modulo pelo nome e abre a janela principal.

O modulo nao pode ser executado como "__main__" (runpy.run_module): os
processos paralelos (ProcessPoolExecutor, modo spawn) recebem as funcoes
por referencia "modulo.funcao" e, dentro do executavel, o "__main__" e este
main.py - o que causava "Can't pickle <function _lsm_r_init>: attribute
lookup _lsm_r_init on __main__ failed" ao calcular metricas em paralelo.
Importado pelo nome, cada processo filho reimporta o modulo e acha a funcao.
O modulo precisa estar em hiddenimports no GeoLand.spec (import dinamico).

Diagnosticos sem interface (gravam um relatorio que termina em RESULTADO):

    GeoLand_V1.exe --teste-motores relatorio.txt
    GeoLand_V1.exe --teste-metricas <pasta_rasters> <pasta_saida> [nucleos]
"""

from __future__ import annotations

import importlib
import multiprocessing
import os
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


def teste_metricas(pasta_rasters: str, pasta_saida: str, nucleos: int | None = None) -> int:
    """Teste sem interface da aba 4 (Metricas), com dados reais:

        GeoLand_V1.exe --teste-metricas <pasta_rasters> <pasta_saida> [nucleos]

    Chama as mesmas funcoes do botao da aba 4, com os dois motores, os tres
    niveis, todas as metricas (botao "Todas"), vizinhanca 8 e subpastas, em
    paralelo. Grava os CSVs e o arquivo relatorio_teste_metricas.txt em
    `pasta_saida`. Falha se algum motor nao processar todos os rasters ou se
    o motor R cair do modo paralelo para o sequencial.
    """
    import time

    saida = Path(pasta_saida)
    saida.mkdir(parents=True, exist_ok=True)
    linhas = ["R_HOME: %s" % os.environ.get("R_HOME"), "Rasters: %s" % pasta_rasters]
    ok = True
    try:
        import pandas as pd

        app = importlib.import_module(MODULE_NAME)
        nucleos = nucleos or app.numero_padrao_de_nucleos()
        total = len(app.listar_rasters_pasta(pasta_rasters, recursivo=True))
        linhas.append("Rasters encontrados: %d | nucleos: %d" % (total, nucleos))
        comuns = dict(
            pasta_rasters=pasta_rasters, pasta_saida=str(saida),
            niveis={"patch": True, "class": True, "landscape": True},
            vizinhanca="8", recursivo=True, max_workers=nucleos,
        )
        motores = [
            ("PyLandStats", app.calcular_metricas_paisagem_pylandstats, dict(
                metricas_por_nivel={"patch": list(app.PATCH_METRICS), "class": list(app.CLASS_METRICS),
                                    "landscape": list(app.LANDSCAPE_METRICS)},
                engine_label="pylandstats")),
            ("landscapemetrics (R)", app.calcular_metricas_paisagem_r, dict(
                metricas_func_por_nivel={n: [f for _, f in m] for n, m in app.R_METRICS_BY_LEVEL.items()},
                cellsize=30.0)),
        ]
        for nome, funcao, extras in motores:
            log: list = []
            inicio = time.time()
            try:
                res = funcao(log_callback=lambda m, _l=log: _l.append(str(m)), **comuns, **extras)
                geral = res.get("csv_geral")
                n_linhas = len(pd.read_csv(geral, usecols=[0])) if geral and Path(str(geral)).is_file() else 0
                n_csv = len(res.get("csvs_individuais") or [])
                problemas = [m for m in log if any(p in m for p in ("Falha", "ERRO", "Erro", "Traceback"))]
                motor_ok = (total > 0 and res.get("rasters_processados") == total and n_csv == total
                            and n_linhas > 0 and not any("Falha no modo paralelo" in m for m in log))
                linhas.append("%s: %s | rasters=%s | CSVs=%d | linhas no CSV geral=%d | avisos=%d | %.0f s"
                              % (nome, "OK" if motor_ok else "FALHA", res.get("rasters_processados"),
                                 n_csv, n_linhas, len(problemas), time.time() - inicio))
                linhas.append("  CSV geral: %s" % geral)
                linhas += ["  aviso: %s" % m[:300] for m in problemas[:5]]
                linhas += ["  modo: %s" % m for m in log if m.startswith("Executando em")]
            except Exception:  # noqa: BLE001
                motor_ok = False
                linhas.append("%s: FALHA\n%s" % (nome, traceback.format_exc()))
            ok = ok and motor_ok
    except Exception:  # noqa: BLE001
        ok = False
        linhas.append(traceback.format_exc())
    linhas.append("RESULTADO: %s" % ("OK" if ok else "FALHA"))
    (saida / "relatorio_teste_metricas.txt").write_text("\n".join(linhas), encoding="utf-8")
    return 0 if ok else 1


def main() -> int:
    # Precisa vir antes de qualquer import do rpy2 e antes de criar processos
    # filhos (que herdam estas variaveis de ambiente).
    configurar_r_embutido()
    # Necessario para o ProcessPoolExecutor (modo spawn) dentro do executavel
    # do PyInstaller: num processo filho, esta chamada executa o trabalho e
    # encerra. Vem antes de ler os argumentos de teste.
    try:
        multiprocessing.freeze_support()
    except OSError:
        if "--multiprocessing-fork" in sys.argv[1:]:
            return 1  # processo filho orfao (o processo pai ja terminou): sai sem dialogo de erro
        raise
    if len(sys.argv) == 3 and sys.argv[1] == "--teste-motores":
        return autoteste_motores(sys.argv[2])
    if len(sys.argv) in (4, 5) and sys.argv[1] == "--teste-metricas":
        return teste_metricas(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) == 5 else None)
    try:
        app = importlib.import_module(MODULE_NAME)
        app.AplicacaoGradeBuffer().mainloop()
        return 0
    except Exception:  # noqa: BLE001
        print("Erro fatal ao iniciar o GeoLand:", file=sys.stderr)
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
