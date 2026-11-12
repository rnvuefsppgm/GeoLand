# GeoLand V1

Aplicativo desktop para Windows (interface gráfica em Tkinter) que calcula
**métricas da estrutura da paisagem** em lote, a partir de rasters de uso e
cobertura da terra, com dois motores de cálculo:

- **PyLandStats** (Python)
- **landscapemetrics** (pacote do R, chamado via `rpy2`)

O executável já traz o Python, o R 4.5.2 e os pacotes `landscapemetrics` e
`terra` embutidos: **não é preciso instalar Python nem R** para usar.

## Download e instalação

1. Baixe o executável:
   **[GeoLand_V1_win64.zip](https://github.com/rnvuefsppgm/GeoLand/raw/main/executavel/GeoLand_V1_win64.zip)**
   (311 MB).
2. Extraia o `.zip` em qualquer pasta (ex.: `C:\GeoLand_V1`).
3. Abra a pasta `GeoLand_V1` e dê duplo clique em **`GeoLand_V1.exe`**.

Use o link acima: o botão **Code → Download ZIP** traz só o código-fonte,
porque o zip do executável é armazenado em Git LFS.

Requisitos: Windows 10 ou 11, 64 bits. Mantenha o `.exe` junto da pasta
`_internal` — copiar só o `.exe` não funciona. Na primeira abertura o Windows
SmartScreen pode avisar que o aplicativo não é assinado; escolha
"Mais informações" → "Executar assim mesmo".

Para conferir se os dois motores funcionam na sua máquina, abra um terminal
na pasta do executável e rode:

```bat
GeoLand_V1.exe --teste-motores relatorio.txt
```

O arquivo `relatorio.txt` deve terminar com `RESULTADO: OK`.

Para testar a aba 4 (Métricas) com os seus próprios rasters, sem abrir a
interface — os dois motores, os três níveis e todas as métricas, em paralelo:

```bat
GeoLand_V1.exe --teste-metricas C:\pasta\dos\rasters C:\pasta\de\saida
```

Os CSVs e o arquivo `relatorio_teste_metricas.txt` são gravados na pasta de
saída; o relatório também deve terminar com `RESULTADO: OK`.

## Tutorial de uso

O fluxo completo vai da criação dos pontos amostrais até a tabela final de
métricas. Cada aba é uma etapa; a pasta de saída de uma etapa serve de
entrada para a seguinte. Você pode começar em qualquer aba, desde que já
tenha os arquivos de entrada dela.

Em todas as abas: use os botões **Selecionar** para escolher os arquivos ou
pastas de entrada e a pasta de saída, ajuste as opções e inicie o
processamento. As barras de **Carregamento** e **Processamento** e o
**Log do processamento** mostram o andamento; o botão **Parar** interrompe.

### Dados de entrada

- **Polígonos da área de estudo** em shapefile (`.shp`).
- **Rasters de uso e cobertura** em GeoTIFF (`.tif`), um por ano, com valores
  inteiros representando as classes. Para métricas de área e distância
  corretas, os rasters devem estar em projeção métrica (UTM); a aba 6 faz
  essa reprojeção.

### Aba 1 — Grade de pontos

Gera uma grade regular de pontos dentro dos polígonos. Informe o shapefile
de polígonos, as **distâncias entre pontos** e a **unidade**. A saída é um
shapefile de pontos.

### Aba 2 — Buffers

Cria buffers circulares em torno dos pontos (da aba 1 ou de um ou mais
shapefiles de pontos seus). Informe os **diâmetros dos buffers** e a
**unidade**. Vários diâmetros geram um buffer por ponto para cada tamanho.

### Aba 3 — Recortar rasters

Recorta cada raster por cada feição do shapefile de buffers. Indique o
**Campo ID** e o **Campo diâmetro** do shapefile, a pasta de rasters e a
pasta de saída. Opções: `all_touched` (inclui pixels tocados pela borda) e
**Núcleos paralelos**. Os recortes são gravados em subpastas por ano.

### Aba 4 — Métricas

Calcula as métricas para todos os rasters de uma pasta.

1. Selecione a pasta de rasters (marque **Ler todas as subpastas** se os
   recortes estiverem em subpastas por ano) e a pasta de saída.
2. Marque os níveis: **Patch**, **Classe** e/ou **Paisagem**.
3. Escolha a **Vizinhança** (4 ou 8 vizinhos).
4. Em **Motores**, marque **PyLandStats**, **landscapemetrics (R)** ou ambos.
   Para o motor R, confira o **cellsize (R)** (tamanho do pixel, 30 por padrão).
5. Nas páginas "Métricas PyLandStats" e "Métricas landscapemetrics (R)",
   escolha as métricas de cada nível (botões **Todas** e **Limpar**).
6. Inicie o cálculo. São gravados um CSV por raster e um CSV geral.

### Aba 5 — Buffer + Métricas (SHP)

Junta a tabela de métricas ao shapefile de buffers pelo **Campo ID**,
gerando um shapefile com as métricas como atributos, pronto para o QGIS.

### Aba 6 — Reprojeção de rasters

Reprojeta em lote os rasters de uma pasta. Informe o **EPSG** de destino
(o botão **Sugerir EPSG** propõe um código), a **Reamostragem** e o número
de núcleos. Use antes das abas 3 e 4 se os rasters estiverem em graus.

### Aba 7 — Fundido (memória)

Faz o recorte (aba 3) e as métricas (aba 4) de uma vez, sem gravar os
recortes em disco — indicado para muitos buffers e muitos anos. Roda em
paralelo e pode **retomar de onde parou**. Pode salvar CSV por feição,
Parquet por ano, Parquet único e shapefile com métricas.

### Aba 8 — Consolidar tabelas

Junta os CSVs (ou Parquets) gerados em uma tabela única. Informe o nome base
do arquivo e os formatos de saída: **CSV único**, **Parquet único** ou
**Parquet por ano**.

## Estrutura do repositório

```
main.py                                          Ponto de entrada; configura o R embutido
app_grade_buffer_shapefile_v18_pylandstats_R.py  Interface e processamento
GeoLand.spec                                     Especificação do PyInstaller
requirements.txt                                 Dependências Python
montar_R_portatil.ps1                            Monta a pasta R_portable (R + pacotes)
build_windows.bat                                Cria o venv e gera o executável
executavel/GeoLand_V1_win64.zip                  Executável compilado (Git LFS)
```

O zip contém a pasta `GeoLand_V1` com o `GeoLand_V1.exe` e a subpasta
`_internal` (Python, bibliotecas e o R portátil), cerca de 760 MB extraídos.

## Executar pelo código-fonte

Requer Python 3.11 (64 bits). O motor R, nesse modo, usa a pasta
`R_portable` se ela existir; senão, o R instalado no sistema.

```bat
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
python main.py
```

## Gerar o executável

Requer Python 3.11 (64 bits), R 4.x (64 bits) instalado e internet.

```bat
build_windows.bat
```

O script cria o venv, monta a pasta `R_portable` (cópia do R instalado, com
`landscapemetrics` e `terra` baixados do CRAN), instala as dependências e
roda o PyInstaller. O resultado fica em `dist\GeoLand_V1\GeoLand_V1.exe`.

## Solução de problemas

| Sintoma | Solução |
|---|---|
| "Failed to load Python DLL" ao abrir | O `.exe` foi separado da pasta `_internal`. Use a pasta `GeoLand_V1` inteira. |
| Antivírus ou SmartScreen bloqueia | O executável não é assinado; libere a execução. |
| Motor R falha | Rode `--teste-motores` e veja o relatório; confirme que `_internal\R` existe. |
| "Can't pickle <function _lsm_r_init>" ou "WinError 6" ao calcular métricas | Erro da primeira compilação da V1, já corrigido: baixe o zip de novo. |
| Áreas ou distâncias estranhas | Rasters em graus: reprojete para UTM na aba 6. |
| Demora para abrir | Normal na primeira execução (o antivírus examina os arquivos). |

## Componentes de terceiros

O pacote distribuído inclui o R, `landscapemetrics`, `terra` e `rpy2`,
licenciados sob a GPL, além de bibliotecas Python de código aberto
(PyLandStats, GeoPandas, Rasterio e outras). Os códigos-fonte desses
componentes estão em https://cran.r-project.org e https://pypi.org.
