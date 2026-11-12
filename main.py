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
    GeoLand_V1.exe --teste-gui <pasta_rasters_em_metros> <pasta_saida> [por_subpasta]
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


def teste_gui(pasta_rasters: str, pasta_saida: str, por_pasta: int = 3) -> int:
    """Teste das abas 6 e 4 pela propria interface, sem cliques manuais:

        GeoLand_V1.exe --teste-gui <pasta_rasters_em_metros> <pasta_saida> [rasters_por_subpasta]

    Copia uma amostra dos rasters (mantendo as subpastas), abre a janela,
    troca os dialogos por respostas automaticas e aciona os botoes:
    aba 6 - seleciona a PASTA da amostra e reprojeta para graus (EPSG 4674),
    conferindo a leitura das subpastas e a estrutura replicada na saida;
    aba 4 - "Verificar projecao (metros)" e, com "Reprojetar para metros"
    marcado, o calculo com os dois motores. Grava relatorio_teste_gui.txt.
    """
    import shutil
    import threading
    import time
    from tkinter import messagebox

    saida = Path(pasta_saida)
    saida.mkdir(parents=True, exist_ok=True)
    linhas: list = []
    falhas: list = []
    dialogos: list = []

    def conferir(nome, condicao, detalhe=""):
        linhas.append("%s: %s | %s" % (nome, "OK" if condicao else "FALHA", detalhe))
        if not condicao:
            falhas.append(nome)

    def estrutura(base):
        return sorted(str(p.relative_to(base).with_suffix("")).lower() for p in Path(base).rglob("*.tif"))

    try:
        import pandas as pd
        import rasterio

        app = importlib.import_module(MODULE_NAME)
        origem, amostra = Path(pasta_rasters), saida / "amostra"
        por_dir: dict = {}
        for r in app.listar_rasters_pasta(str(origem), recursivo=True):
            por_dir.setdefault(r.parent, []).append(r)
        n = 0
        for pasta, arquivos in por_dir.items():
            destino = amostra / pasta.relative_to(origem)
            destino.mkdir(parents=True, exist_ok=True)
            for r in arquivos[:por_pasta]:
                shutil.copy2(r, destino / r.name)
                n += 1
        linhas.append("Amostra: %d rasters em %d pasta(s)" % (n, len(por_dir)))

        def area_ha(base):
            total = 0.0
            for p in Path(base).rglob("*.tif"):
                with rasterio.open(p) as src:
                    banda = src.read(1)
                    validos = int((banda != src.nodata).sum()) if src.nodata is not None else banda.size
                    total += validos * abs(src.transform.a * src.transform.e) / 10000.0
            return total

        for nome in ("showinfo", "showwarning", "showerror"):
            setattr(messagebox, nome, lambda titulo="", msg="", _n=nome, **k: dialogos.append(
                "%s: %s" % (_n, str(msg).replace("\n", " / "))))
        janela = app.AplicacaoGradeBuffer()

        def acionar(botao, metodo, limite=7200):
            """Clica no botao e espera o processamento (em outra thread) terminar."""
            fim = threading.Event()
            original = getattr(janela, metodo)

            def embrulho():
                try:
                    original()
                finally:
                    fim.set()

            setattr(janela, metodo, embrulho)
            try:
                for _ in range(100):  # espera o botao ser reabilitado pelo processamento anterior
                    if str(botao.cget("state")) != "disabled":
                        break
                    time.sleep(0.1)
                botao.invoke()
                terminou = fim.wait(limite)
            finally:
                delattr(janela, metodo)
            time.sleep(1.0)  # deixa a janela gravar as ultimas linhas do log
            return terminou

        def roteiro():
            try:
                # Aba 6: selecionar a PASTA (dialogo trocado) e reprojetar para graus.
                graus = saida / "aba6_graus"
                janela._selecionar_pasta = lambda *a, **k: str(amostra)
                janela._acao_selecionar_rasters_reproj()
                conferir("Aba 6 - pasta e subpastas carregadas", len(janela.rasters_reproj_lista) == n,
                         "%d rasters | saida sugerida: %s | EPSG sugerido: %s" % (
                             len(janela.rasters_reproj_lista), janela.pasta_saida_reproj_var.get(),
                             janela.epsg_reproj_var.get()))
                janela.pasta_saida_reproj_var.set(str(graus))
                janela.epsg_reproj_var.set("4674")
                janela.sufixo_epsg_reproj_var.set(False)
                acionar(janela.botao_reproj, "_processar_reprojecao")
                conferir("Aba 6 - subpastas replicadas na saida", n > 0 and estrutura(graus) == estrutura(amostra),
                         "%d rasters em %d pasta(s)" % (len(estrutura(graus)),
                                                        len({p.parent for p in graus.rglob("*.tif")})))

                # Aba 4: botao de verificacao, com rasters em graus e com rasters em metros.
                csvs = saida / "aba4_csv"
                janela.pasta_saida_metricas_var.set(str(csvs))
                for pasta, esperado in ((graus, n), (amostra, 0)):
                    janela.pasta_rasters_metricas_var.set(str(pasta))
                    acionar(janela.botao_verificar_proj, "_verificar_projecao_metricas")
                    log = janela.log_metricas.get("1.0", "end")
                    conferir("Aba 4 - Verificar projecao (%s)" % pasta.name,
                             log.count("NÃO está em metros:") == esperado,
                             "; ".join(l for l in log.splitlines() if l.startswith("Rasters verificados")))

                # Aba 4: reprojetar para metros antes do calculo, com os dois motores.
                janela.pasta_rasters_metricas_var.set(str(graus))
                janela.reprojetar_metricas_var.set(True)
                janela.engine_pls_var.set(True)
                janela.engine_r_var.set(True)
                janela._atualizar_estado_metricas()
                acionar(janela.botao_metricas, "_processar_metricas")
                metros = csvs / app.PASTA_REPROJ_METRICAS
                proj = app.verificar_projecao_rasters(sorted(metros.rglob("*.tif")))
                a0, a1 = area_ha(amostra), (area_ha(metros) if metros.is_dir() else 0.0)
                conferir("Aba 4 - reprojecao para metros antes do calculo",
                         n > 0 and estrutura(metros) == estrutura(graus) and len(proj["metros"]) == n
                         and a0 > 0 and abs(a1 / a0 - 1) < 0.05,
                         "%d rasters em metros | area valida: original %.1f ha, reprojetado %.1f ha" % (
                             len(proj["metros"]), a0, a1))
                for motor in ("pylandstats", "landscapemetrics"):
                    gerais = sorted(csvs.glob("metricas_paisagem_GERAL_%s_*.csv" % motor))
                    df = pd.read_csv(gerais[-1], usecols=["source_file"]) if gerais else None
                    conferir("Aba 4 - metricas %s" % motor, df is not None and df["source_file"].nunique() == n,
                             "sem CSV geral" if df is None else "%d rasters, %d linhas" % (
                                 df["source_file"].nunique(), len(df)))
                log = janela.log_metricas.get("1.0", "end")
                conferir("Aba 4 - motor R em processos paralelos", "Falha no modo paralelo" not in log,
                         "; ".join(l for l in log.splitlines() if l.startswith(("Executando em", "Rasters prontos"))))
            except Exception:  # noqa: BLE001
                falhas.append("roteiro")
                linhas.append(traceback.format_exc())
            finally:
                janela.after(0, janela.destroy)

        janela.after(500, lambda: threading.Thread(target=roteiro, daemon=True).start())
        janela.mainloop()
        conferir("Dialogos de erro", not any(d.startswith("showerror") for d in dialogos), "%d dialogo(s)" % len(dialogos))
        linhas += ["  dialogo %s" % d[:300] for d in dialogos]
    except Exception:  # noqa: BLE001
        falhas.append("teste")
        linhas.append(traceback.format_exc())
    linhas.append("RESULTADO: %s" % ("OK" if not falhas else "FALHA"))
    (saida / "relatorio_teste_gui.txt").write_text("\n".join(linhas), encoding="utf-8")
    return 0 if not falhas else 1


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
    if len(sys.argv) in (4, 5) and sys.argv[1] == "--teste-gui":
        return teste_gui(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) == 5 else 3)
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
