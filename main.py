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
    GeoLand_V1.exe --teste-graficos <tabela.csv> <pasta_saida>
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
    """Teste sem interface da aba 6 (Metricas), com dados reais:

        GeoLand_V1.exe --teste-metricas <pasta_rasters> <pasta_saida> [nucleos]

    Chama as mesmas funcoes do botao da aba 6, com os dois motores, os tres
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
    """Teste das abas 4 e 6 pela propria interface, sem cliques manuais:

        GeoLand_V1.exe --teste-gui <pasta_rasters_em_metros> <pasta_saida> [rasters_por_subpasta]

    Copia uma amostra dos rasters (mantendo as subpastas), abre a janela,
    troca os dialogos por respostas automaticas e aciona os botoes:
    aba 4 - seleciona a PASTA da amostra e reprojeta para graus (EPSG 4674),
    conferindo a leitura das subpastas e a estrutura replicada na saida;
    aba 6 - "Verificar projecao (metros)" e, com "Reprojetar para metros"
    marcado, o calculo com os dois motores;
    aba 5 - escaneia a pasta, cria/salva/carrega regras e reclassifica,
    conferindo o resultado pixel a pixel. Grava relatorio_teste_gui.txt.
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
        import numpy as np
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
                # Aba 4: selecionar a PASTA (dialogo trocado) e reprojetar para graus.
                graus = saida / "aba6_graus"
                janela._selecionar_pasta = lambda *a, **k: str(amostra)
                janela._acao_selecionar_rasters_reproj()
                conferir("Aba 4 - pasta e subpastas carregadas", len(janela.rasters_reproj_lista) == n,
                         "%d rasters | saida sugerida: %s | EPSG sugerido: %s" % (
                             len(janela.rasters_reproj_lista), janela.pasta_saida_reproj_var.get(),
                             janela.epsg_reproj_var.get()))
                janela.pasta_saida_reproj_var.set(str(graus))
                janela.epsg_reproj_var.set("4674")
                janela.sufixo_epsg_reproj_var.set(False)
                acionar(janela.botao_reproj, "_processar_reprojecao")
                conferir("Aba 4 - subpastas replicadas na saida", n > 0 and estrutura(graus) == estrutura(amostra),
                         "%d rasters em %d pasta(s)" % (len(estrutura(graus)),
                                                        len({p.parent for p in graus.rglob("*.tif")})))

                # Aba 6: botao de verificacao, com rasters em graus e com rasters em metros.
                csvs = saida / "aba4_csv"
                janela.pasta_saida_metricas_var.set(str(csvs))
                for pasta, esperado in ((graus, n), (amostra, 0)):
                    janela.pasta_rasters_metricas_var.set(str(pasta))
                    acionar(janela.botao_verificar_proj, "_verificar_projecao_metricas")
                    log = janela.log_metricas.get("1.0", "end")
                    conferir("Aba 6 - Verificar projecao (%s)" % pasta.name,
                             log.count("NÃO está em metros:") == esperado,
                             "; ".join(l for l in log.splitlines() if l.startswith("Rasters verificados")))

                # Aba 6: reprojetar para metros antes do calculo, com os dois motores.
                janela.pasta_rasters_metricas_var.set(str(graus))
                janela.reprojetar_metricas_var.set(True)
                janela.engine_pls_var.set(True)
                janela.engine_r_var.set(True)
                janela.calcular_iep_var.set(True)
                janela._atualizar_estado_metricas()
                acionar(janela.botao_metricas, "_processar_metricas")
                metros = csvs / app.PASTA_REPROJ_METRICAS
                proj = app.verificar_projecao_rasters(sorted(metros.rglob("*.tif")))
                a0, a1 = area_ha(amostra), (area_ha(metros) if metros.is_dir() else 0.0)
                conferir("Aba 6 - reprojecao para metros antes do calculo",
                         n > 0 and estrutura(metros) == estrutura(graus) and len(proj["metros"]) == n
                         and a0 > 0 and abs(a1 / a0 - 1) < 0.05,
                         "%d rasters em metros | area valida: original %.1f ha, reprojetado %.1f ha" % (
                             len(proj["metros"]), a0, a1))
                for motor in ("pylandstats", "landscapemetrics"):
                    gerais = sorted(csvs.glob("metricas_paisagem_GERAL_%s_*.csv" % motor))
                    df = pd.read_csv(gerais[-1], usecols=["source_file"]) if gerais else None
                    conferir("Aba 6 - metricas %s" % motor, df is not None and df["source_file"].nunique() == n,
                             "sem CSV geral" if df is None else "%d rasters, %d linhas" % (
                                 df["source_file"].nunique(), len(df)))
                # Indice de estabilidade: um CSV por motor, calculado depois das metricas de todos os anos.
                ieps = {m: sorted(csvs.glob("indice_estabilidade_IEP_%s_*.csv" % m)) for m in ("pylandstats", "landscapemetrics")}
                tabelas = {m: pd.read_csv(a[-1]) for m, a in ieps.items() if a}
                conferir("Aba 6 - indice de estabilidade (IEP) dos dois motores",
                         len(tabelas) == 2 and all(len(t) > 0 and t["IEP_u"].between(0, 1).all()
                                                   and t["n_anos"].max() == len(por_dir) for t in tabelas.values()),
                         "; ".join("%s: %d combinacoes, ate %d anos, IEP medio %.3f" % (
                             m, len(t), t["n_anos"].max(), t["IEP_u"].mean()) for m, t in tabelas.items()))
                log = janela.log_metricas.get("1.0", "end")
                conferir("Aba 6 - motor R em processos paralelos", "Falha no modo paralelo" not in log,
                         "; ".join(l for l in log.splitlines() if l.startswith(("Executando em", "Rasters prontos"))))

                # Aba 5: escanear a pasta (com subpastas), criar uma regra, salvar e carregar
                # o .txt, atribuir a regra carregada a 3 arquivos e reclassificar tudo.
                from tkinter import filedialog

                reclass, txt = saida / "aba9_reclass", saida / "regra_teste.txt"
                janela.pasta_saida_reclass_var.set(str(reclass))
                janela._selecionar_pasta = lambda *a, **k: str(amostra)
                janela._acao_selecionar_pasta_entrada_reclass()
                conferir("Aba 5 - arquivos escaneados (com subpastas)", len(janela.reclass_arquivos) == n,
                         "%d arquivos" % len(janela.reclass_arquivos))
                janela.nova_regra_reclass_var.set("Teste")
                janela._acao_criar_regra_reclass()
                arvore = janela.tree_regras_reclass
                arvore.delete(*arvore.get_children())
                arvore.insert("", "end", values=(1, 100, 7))
                filedialog.asksaveasfilename = lambda **k: str(txt)
                filedialog.askopenfilename = lambda **k: str(txt)
                janela._acao_salvar_txt_reclass()
                janela._acao_carregar_txt_reclass()  # vira a regra ativa "regra_teste"
                conferir("Aba 5 - salvar e carregar regras em .txt",
                         janela.reclass_regras.get("regra_teste") == [(1.0, 100.0, 7.0)],
                         str(janela.reclass_regras.get("regra_teste")))
                # Ao carregar o .txt, a regra passa a valer para os arquivos que estavam na "Padrao".
                carregada = sum(1 for arq in janela.reclass_arquivos if arq["Regra"] == "regra_teste")
                conferir("Aba 5 - regra carregada aplicada aos arquivos", carregada == n,
                         "%d de %d arquivos com a regra carregada" % (carregada, n))
                # Volta 3 arquivos para a regra "Padrao" (mistura de regras na mesma execucao).
                janela.regra_ativa_reclass_var.set("Padrao")
                janela._acao_selecionar_regra_reclass()
                itens = janela.tree_arquivos_reclass.get_children()
                janela.tree_arquivos_reclass.selection_set(itens[:3])
                janela._acao_atribuir_selecionados_reclass()
                # Tabela de atributos antiga do ArcGIS ao lado da saida (de uma execucao anterior):
                # tem de ser apagada, senao o ArcGIS mostra as classes antigas.
                rel0 = Path(janela.reclass_arquivos[0]["Arquivo"])
                velhos = [reclass / rel0.parent / ("Reclass_" + rel0.name + s) for s in (".vat.dbf", ".aux.xml")]
                for v in velhos:
                    v.parent.mkdir(parents=True, exist_ok=True)
                    v.write_text("antigo")
                acionar(janela.botao_reclass, "_processar_reclass")
                conferir("Aba 5 - tabela de atributos antiga (.vat.dbf/.aux.xml) apagada",
                         not any(v.exists() for v in velhos), ", ".join(v.name for v in velhos if v.exists()))
                # Conferencia independente, pixel a pixel (na entrada, 0 = sem dado).
                esperado = {"regra_teste": np.full(256, 255, dtype="uint8"), "Padrao": np.full(256, 255, dtype="uint8")}
                esperado["regra_teste"][1:101] = 7
                for minimo, _maximo, novo in app.REGRAS_RECLASS_PADRAO:
                    esperado["Padrao"][int(minimo)] = novo
                diferentes, meta_ok = 0, True
                for arq in janela.reclass_arquivos:
                    rel = Path(arq["Arquivo"])
                    with rasterio.open(arq["Caminho"]) as ent, rasterio.open(
                            reclass / rel.parent / ("Reclass_" + rel.name)) as sai:
                        diferentes += int((esperado[arq["Regra"]][ent.read(1)] != sai.read(1)).sum())
                        meta_ok = (meta_ok and sai.dtypes[0] == "uint8" and sai.nodata == 255
                                   and sai.crs == ent.crs and sai.transform == ent.transform)
                com_teste = sum(1 for arq in janela.reclass_arquivos if arq["Regra"] == "regra_teste")
                conferir("Aba 5 - reclassificacao conferida pixel a pixel",
                         len(estrutura(reclass)) == n and diferentes == 0 and meta_ok and com_teste == n - 3,
                         "%d rasters (%d com a regra carregada), %d pixels diferentes, subpastas replicadas" % (
                             len(estrutura(reclass)), com_teste, diferentes))
                # Carregar o .txt ANTES de escolher a pasta: os arquivos recebem a regra do editor.
                janela.reclass_arquivos = []
                janela._acao_carregar_txt_reclass()
                janela._acao_selecionar_pasta_entrada_reclass()
                antes = sum(1 for arq in janela.reclass_arquivos if arq["Regra"] == "regra_teste")
                conferir("Aba 5 - .txt carregado antes da pasta vale para os arquivos", antes == n,
                         "%d de %d arquivos com a regra carregada" % (antes, n))
                conferir("Aba 5 - ano lido do nome com prefixo",
                         app.extrair_ano_de_recorte("Reclass_clip_id_APA_10_diam_2000_1985.tif") == "1985",
                         app.extrair_ano_de_recorte("Reclass_clip_id_APA_10_diam_2000_1985.tif"))
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


def teste_graficos(tabela: str, pasta_saida: str) -> int:
    """Teste da aba 10 (Graficos) pela propria interface, sem cliques manuais:

        GeoLand_V1.exe --teste-graficos <tabela.csv> <pasta_saida>

    Le a tabela pelos controles da aba e confere os campos criados
    automaticamente; usa filtros; edita o dicionario de cores (o seletor de
    cores e trocado por uma cor fixa), salva e recarrega o .txt; gera um lote
    pequeno em PNG, PDF e SVG; salva um grafico individual; e repete com uma
    tabela generica (ponto e virgula, virgula decimal, outros nomes de
    coluna). Grava relatorio_teste_graficos.txt em `pasta_saida`.
    """
    import re
    import threading
    import time
    from tkinter import colorchooser, filedialog, messagebox

    import pandas as pd

    saida = Path(pasta_saida)
    saida.mkdir(parents=True, exist_ok=True)
    linhas: list = []
    falhas: list = []
    dialogos: list = []

    def conferir(nome, condicao, detalhe=""):
        linhas.append("%s: %s | %s" % (nome, "OK" if condicao else "FALHA", detalhe))
        if not condicao:
            falhas.append(nome)

    try:
        app = importlib.import_module(MODULE_NAME)
        for nome in ("showinfo", "showwarning", "showerror"):
            setattr(messagebox, nome, lambda titulo="", msg="", _n=nome, **k: dialogos.append(
                "%s: %s" % (_n, str(msg).replace("\n", " / "))))
        janela = app.AplicacaoGradeBuffer()
        aba = janela.graficos

        def na_janela(funcao):
            """Executa `funcao` na thread da interface e devolve o resultado."""
            feito, caixa = threading.Event(), {}

            def rodar():
                try:
                    caixa["r"] = funcao()
                except Exception:  # noqa: BLE001
                    caixa["e"] = traceback.format_exc()
                finally:
                    feito.set()

            janela.after(0, rodar)
            feito.wait(300)
            if "e" in caixa:
                raise RuntimeError(caixa["e"])
            return caixa.get("r")

        def acionar(botao, metodo, limite=3600, alvo=None):
            """Clica no botao e espera o processamento (em outra thread) terminar."""
            alvo = alvo or aba
            fim, original = threading.Event(), getattr(alvo, metodo)

            def embrulho():
                try:
                    original()
                finally:
                    fim.set()

            setattr(alvo, metodo, embrulho)
            try:
                for _ in range(100):
                    if str(botao.cget("state")) != "disabled":
                        break
                    time.sleep(0.1)
                botao.invoke()
                fim.wait(limite)
            finally:
                delattr(alvo, metodo)
            time.sleep(1.0)

        def resumo_do_lote(alvo=None):
            achado = re.search(r"Gráficos gerados: (\d+) de (\d+) \| sem dados: (\d+)", (alvo or aba).log.get("1.0", "end"))
            return tuple(int(v) for v in achado.groups()) if achado else (0, 0, 0)

        def roteiro():
            try:
                # 1) Tabela de exemplo: leitura e campos automaticos.
                na_janela(lambda: (janela.notebook.select(janela.aba_graficos.master.master),
                                   aba.tabela_var.set(tabela), aba.saida_var.set(str(saida / "lote"))))
                acionar(aba.botao_carregar, "_ler_e_estruturar")
                est = aba.estrutura
                conferir("Leitura e campos automaticos", aba.df is not None and bool(est.x) and len(est.metricas) > 0,
                         "%d linhas | X=%s | linhas=%s | um grafico por=%s | metricas=%d | filtros=%s" % (
                             0 if aba.df is None else len(aba.df), est.x, est.serie, "+".join(est.painel),
                             len(est.metricas), ",".join(str(c) for c in aba.filtros)))
                tem_cores = (Path(tabela).parent / "cores.txt").is_file()
                conferir("Dicionario de cores lido da pasta da tabela", bool(aba.legenda) == tem_cores,
                         "%d valores" % len(aba.legenda))
                imagem = aba._imagem_previa
                conferir("Pre-visualizacao apos a leitura", imagem is not None and imagem.width() > 200,
                         "sem imagem" if imagem is None else "%dx%d px" % (imagem.width(), imagem.height()))

                # 2) Selecao: 3 metricas e, pelo filtro, 2 valores do primeiro campo de "um grafico por".
                metricas = aba.lista_metricas.itens()[:3]

                def selecionar():
                    aba.lista_metricas.marcar(metricas)
                    if est.painel:
                        aba.definir_filtro(est.painel[0], aba._valores_dim[est.painel[0]][:2])
                    aba.pdf_var.set(True)
                    aba.svg_var.set(True)
                    return len(aba._combos_previa)

                esperado = na_janela(selecionar) * len(metricas)

                # 3) Dicionario de cores: paleta, escolher uma cor, salvar e recarregar o .txt.
                dicionario = saida / "dicionario_cores.txt"
                colorchooser.askcolor = lambda **k: ((18, 52, 86), "#123456")
                filedialog.asksaveasfilename = lambda **k: str(dicionario)
                filedialog.askopenfilename = lambda **k: str(dicionario)

                def cores():
                    dialogo = aba.abrir_dialogo_cores()
                    chaves = list(dialogo._chaves)
                    dialogo.paleta_var.set("Viridis")
                    dialogo.aplicar_paleta()
                    dialogo.tree.selection_set("0")
                    dialogo.escolher_cor()
                    dialogo.salvar()
                    aba.legenda = {}  # esquece tudo e recarrega do arquivo salvo
                    dialogo.carregar()
                    usadas = aba.plano().cores
                    dialogo.destroy()
                    return chaves, usadas

                chaves, usadas = na_janela(cores)
                conferir("Dicionario de cores: escolher, salvar e carregar",
                         bool(chaves) and usadas.get(chaves[0]) == "#123456" and dicionario.is_file()
                         and len({usadas[k] for k in chaves}) == len(chaves),
                         "%d classes com cores distintas | cor escolhida para %s: %s" % (
                             len(chaves), chaves[0] if chaves else "-", usadas.get(chaves[0]) if chaves else "-"))

                # 4) Botao "Gerar TODOS os graficos".
                acionar(aba.botao_gerar, "_gerar_todos")
                gerados, total, sem_dados = resumo_do_lote()
                arquivos = [p for p in (saida / "lote").rglob("*") if p.is_file()]
                pngs = [p for p in arquivos if p.suffix == ".png"]
                validos = all(p.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n" and p.stat().st_size > 15000 for p in pngs)
                largura = int.from_bytes(pngs[0].read_bytes()[16:20], "big") if pngs else 0
                conferir("Lote de graficos (PNG, PDF e SVG)",
                         gerados > 0 and total == esperado and gerados + sem_dados == total and len(pngs) == gerados
                         and len(arquivos) == 3 * gerados and validos,
                         "%d graficos de %d previstos (%d sem dados) x 3 formatos | PNG com %d px de largura" % (
                             gerados, total, sem_dados, largura))

                # 5) Grafico individual (o da pre-visualizacao).
                def individual():
                    aba.pdf_var.set(False)
                    aba.svg_var.set(False)
                    aba.saida_var.set(str(saida / "individual"))
                    aba.salvar_individual()

                na_janela(individual)
                conferir("Grafico individual salvo", len(list((saida / "individual").rglob("*.png"))) == 1,
                         "%d arquivo(s)" % len(list((saida / "individual").rglob("*.*"))))

                # 5b) Indice de estabilidade: caixa propria, com os cinco resultados e eixo X = Raio.
                def marcar_indice():
                    for var in aba.indice_vars.values():
                        var.set(True)
                    aba.saida_var.set(str(saida / "indice"))
                    return len(aba.filtros.get(est.painel[0]) or []) if est.painel else 0

                na_janela(marcar_indice)
                acionar(aba.botao_indice, "_gerar_indice")
                pasta_indice = saida / "indice" / "indice_estabilidade"
                csv_indice = pasta_indice / "indice_estabilidade_IEP.csv"
                tabela_indice = pd.read_csv(csv_indice) if csv_indice.is_file() else None
                # Um grafico por resultado e por buffer da tabela do indice (o painel da
                # tabela de metricas pode ser Raio, quando ha um buffer so).
                buffers = int(tabela_indice["Buffer"].nunique()) if tabela_indice is not None and "Buffer" in tabela_indice else 0
                pngs_indice = list(pasta_indice.rglob("*.png")) if pasta_indice.is_dir() else []
                conferir("Indice de estabilidade: tabela e graficos",
                         tabela_indice is not None and len(tabela_indice) > 0 and tabela_indice["IEP_u"].between(0, 1).all()
                         and buffers > 0 and len(pngs_indice) == 5 * buffers,
                         "sem tabela" if tabela_indice is None else "%d combinacoes | IEP de %.3f a %.3f | %d graficos (5 resultados x %d buffers)" % (
                             len(tabela_indice), tabela_indice["IEP_u"].min(), tabela_indice["IEP_u"].max(), len(pngs_indice), buffers))

                # 5c) Caixa "Classes" da aba 10: marcar 2 classes deixa so elas nos graficos.
                serie = aba.coluna_serie()
                classes = na_janela(lambda: [c for c in aba.lista_classes.itens()
                                             if c in set(aba.tabela_filtrada()[serie].dropna())]) if serie else []

                def marcar_classes():
                    aba.lista_classes.marcar(classes[:2])
                    vistas = sorted(set(aba.tabela_filtrada()[serie].tolist()))
                    filtro = aba.filtros.get(serie)
                    aba.lista_classes.marcar_todos()
                    return vistas, filtro, aba.filtros.get(serie)

                vistas, filtro, depois = na_janela(marcar_classes) if serie and len(classes) >= 2 else ([], None, None)
                conferir("Aba 10 - caixa Classes escolhe as linhas dos graficos",
                         vistas == sorted(classes[:2]) and filtro == set(classes[:2]) and depois is None,
                         "coluna %s | %d classes | marcadas %s -> no grafico %s" % (serie, len(classes), classes[:2], vistas))

                # 5d) Aba 11 (so o indice): a partir da tabela de metricas e da planilha do indice.
                ind = janela.graficos_indice
                na_janela(lambda: (janela.notebook.select(janela.aba_graficos_indice.master.master),
                                   ind.tabela_var.set(tabela), ind.saida_var.set(str(saida / "iep_metricas"))))
                acionar(ind.botao_carregar, "_ler_e_estruturar", alvo=ind)
                e = ind.estrutura
                conferir("Aba 11 - tabela de metricas: indice calculado e campos",
                         ind.iep_calculado is not None and e.x == "Raio" and e.metricas[:1] == ["IEP_u"]
                         and ind._imagem_previa is not None,
                         "X=%s | linhas=%s | um grafico por=%s | metricas=%s" % (
                             e.x, e.serie, "+".join(e.painel), ",".join(map(str, e.metricas[:5]))))
                resultados = ("IEP_u", "V_u", "V_comp", "V_conf", "V_conn")
                na_janela(lambda: (ind.lista_classes.marcar(ind.lista_classes.itens()[:2]),
                                   ind.lista_metricas.marcar([c for c in ind.lista_metricas.itens() if c in resultados])))
                acionar(ind.botao_gerar, "_gerar_todos", alvo=ind)
                gerados, total, _sem = resumo_do_lote(ind)
                csv_iep = saida / "iep_metricas" / "indice_estabilidade_IEP.csv"
                pngs_iep = list((saida / "iep_metricas").rglob("*.png"))
                conferir("Aba 11 - lote do indice (5 resultados, 2 classes)",
                         gerados == total > 0 and len(pngs_iep) == gerados and csv_iep.is_file(),
                         "%d de %d graficos | tabela do indice gravada: %s" % (gerados, total, csv_iep.is_file()))
                if csv_iep.is_file():
                    na_janela(lambda: (ind.tabela_var.set(str(csv_iep)), ind.saida_var.set(str(saida / "iep_planilha"))))
                    acionar(ind.botao_carregar, "_ler_e_estruturar", alvo=ind)
                    conferir("Aba 11 - planilha do indice usada diretamente",
                             ind.iep_calculado is None and ind.estrutura.x == "Raio" and ind.df is not None
                             and len(ind.df) == len(pd.read_csv(csv_iep)) and ind._imagem_previa is not None,
                             "%d linhas | X=%s | linhas=%s" % (0 if ind.df is None else len(ind.df), ind.estrutura.x,
                                                              ind.estrutura.serie))

                # 6) Outra tabela: ponto e virgula, virgula decimal e nomes de coluna diferentes.
                generica = saida / "tabela_generica.csv"
                texto = ["Municipio;Ano;Cultura;Produtividade;Area_plantada"]
                for i, municipio in enumerate(("Ilheus", "Una", "Canavieiras")):
                    for j, cultura in enumerate(("Cacau", "Banana", "Mandioca")):
                        for ano in range(2010, 2024):
                            texto.append("%s;%d;%s;%s;%s" % (
                                municipio, ano, cultura, ("%.1f" % (300 + 90 * i + 60 * j + 7.5 * (ano - 2010))).replace(".", ","),
                                ("%.2f" % (5 + 3 * i + 2 * j + 0.4 * (ano - 2010))).replace(".", ",")))
                generica.write_text("\n".join(texto), encoding="utf-8")
                na_janela(lambda: (aba.tabela_var.set(str(generica)), aba.saida_var.set(str(saida / "generica")),
                                   aba.legenda_var.set("")))
                acionar(aba.botao_carregar, "_ler_e_estruturar")
                est = aba.estrutura
                conferir("Tabela generica: campos automaticos",
                         est.x == "Ano" and est.serie == "Cultura" and est.painel == ["Municipio"]
                         and est.metricas == ["Produtividade", "Area_plantada"],
                         "X=%s | linhas=%s | um grafico por=%s | metricas=%s" % (
                             est.x, est.serie, "+".join(est.painel), ",".join(est.metricas)))
                na_janela(lambda: (aba.lista_metricas.marcar_todos(), aba.tema_var.set("Branco com grade (theme_bw)"),
                                   aba.posicao_legenda_var.set("Abaixo"), aba.virgula_var.set(True)))
                acionar(aba.botao_gerar, "_gerar_todos")
                gerados, total, _sem = resumo_do_lote()
                conferir("Tabela generica: lote", gerados == total == 6
                         and len(list((saida / "generica").rglob("*.png"))) == 6, "%d de %d graficos" % (gerados, total))
            except Exception:  # noqa: BLE001
                falhas.append("roteiro")
                linhas.append(traceback.format_exc())
            finally:
                janela.after(0, janela.destroy)

        janela.after(500, lambda: threading.Thread(target=roteiro, daemon=True).start())
        janela.mainloop()
        conferir("Dialogos de erro", not any(d.startswith("showerror") for d in dialogos), "%d dialogo(s)" % len(dialogos))
        linhas += ["  dialogo %s" % d[:200] for d in dialogos if d.startswith("showerror")]
    except Exception:  # noqa: BLE001
        falhas.append("teste")
        linhas.append(traceback.format_exc())
    linhas.append("RESULTADO: %s" % ("OK" if not falhas else "FALHA"))
    (saida / "relatorio_teste_graficos.txt").write_text("\n".join(linhas), encoding="utf-8")
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
    if len(sys.argv) == 4 and sys.argv[1] == "--teste-graficos":
        return teste_graficos(sys.argv[2], sys.argv[3])
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
