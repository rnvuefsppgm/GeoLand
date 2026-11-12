@echo off
REM ============================================================
REM build_windows.bat
REM Cria o ambiente virtual, instala dependencias e gera o
REM executavel AppGrade.exe usando PyInstaller (Windows 10/11, 64 bits).
REM Nao depende de Conda nem de internet apos a instalacao inicial.
REM ============================================================

setlocal

echo ============================================================
echo  AppGrade - Grade, Buffers, Recorte, Metricas, Reprojecao,
echo  Fundido e Consolidar
echo  Build para Windows
echo ============================================================

echo.
echo [1/5] Criando ambiente virtual em .\venv com Python 3.11 ...

REM Removemos um venv antigo/quebrado (ex.: criado com a versao errada do
REM Python) para evitar reaproveitar um ambiente inconsistente.
if exist venv (
    echo Removendo ambiente virtual existente...
    rmdir /s /q venv
)

REM Tentamos primeiro o Python Launcher pedindo explicitamente a versao
REM 3.11, pois em maquinas com varias versoes instaladas (3.12, 3.13, 3.14...)
REM o comando "python" isolado pode apontar para uma versao mais nova sem
REM wheels pre-compiladas para numpy/scipy/geopandas/rasterio, forcando
REM compilacao do zero (o que costuma falhar sem compilador C/C++).
where py >nul 2>nul
if not errorlevel 1 (
    py -3.11 -m venv venv 2>nul
)

if not exist venv\Scripts\python.exe (
    echo [AVISO] Nao foi possivel criar o venv com "py -3.11". Tentando com
    echo         o Python padrao do PATH ^(comando "python"^)...
    where python >nul 2>nul
    if errorlevel 1 (
        echo [ERRO] Nenhum Python encontrado ^(nem "py -3.11" nem "python"^).
        echo        Instale o Python 3.11 ^(64 bits^) em
        echo        https://www.python.org/downloads/ marcando "Add python.exe to PATH".
        exit /b 1
    )
    python -m venv venv
)

if not exist venv\Scripts\python.exe (
    echo [ERRO] Falha ao criar o ambiente virtual.
    exit /b 1
)

echo.
echo Verificando a versao do Python usada no venv...
venv\Scripts\python.exe --version
venv\Scripts\python.exe -c "import sys; exit(0 if sys.version_info[:2]==(3,11) else 1)"
if errorlevel 1 (
    echo.
    echo [AVISO] O venv NAO foi criado com Python 3.11. Isso pode causar
    echo         falhas ao instalar numpy/scipy/geopandas/rasterio
    echo         ^(compilacao do zero em vez de usar wheels prontas^).
    echo         Instale o Python 3.11 ^(64 bits^) de https://www.python.org/downloads/,
    echo         apague a pasta venv e rode este script novamente.
    echo         Continuando mesmo assim...
)

echo.
echo [2/5] Ativando ambiente virtual...
call venv\Scripts\activate.bat

echo.
echo Preparando o R portatil ^(R + landscapemetrics + terra^) em .\R_portable ...
if not exist R_portable\bin\x64\R.dll powershell -NoProfile -ExecutionPolicy Bypass -File montar_R_portatil.ps1
if not exist R_portable\bin\x64\R.dll goto sem_r_portatil
REM O rpy2 precisa encontrar o R durante a instalacao e o build.
set "R_HOME=%CD%\R_portable"
set "PATH=%CD%\R_portable\bin\x64;%PATH%"
goto r_portatil_ok
:sem_r_portatil
echo [AVISO] R portatil nao montado - o .exe sera gerado SEM o R embutido.
:r_portatil_ok

echo.
echo [3/5] Atualizando pip e instalando dependencias...
python -m pip install --upgrade pip
pip install -r requirements.txt
if errorlevel 1 (
    echo [AVISO] Falha ao instalar TODAS as dependencias. Isso normalmente
    echo         acontece por causa do rpy2 ^(motor R, opcional^) em maquinas
    echo         sem R instalado ou sem compilador C/C++. Veja a secao
    echo         "Motor R" no README.md.
    echo.
    echo         Tentando novamente SEM o rpy2 para nao bloquear o build...
    findstr /v /b "rpy2" requirements.txt > requirements.sem_r.txt
    pip install -r requirements.sem_r.txt
    if errorlevel 1 (
        echo [ERRO] Falha ao instalar dependencias mesmo sem o rpy2. Veja a
        echo        secao de solucao de problemas no README.md.
        exit /b 1
    )
    del requirements.sem_r.txt
)

echo.
echo [4/5] Verificando importacao dos modulos principais...
python -c "import geopandas, rasterio, pyproj, shapely, pylandstats; print('OK: dependencias principais importadas com sucesso')"
if errorlevel 1 (
    echo [ERRO] Falha ao importar dependencias principais. Corrija antes de
    echo        gerar o executavel.
    exit /b 1
)
python -c "import rpy2" 2>nul
if errorlevel 1 (
    echo [INFO] rpy2 nao disponivel neste venv - o .exe sera gerado SEM o
    echo        motor R ^(landscapemetrics^). O motor PyLandStats continua
    echo        funcionando normalmente.
) else (
    echo [INFO] rpy2 disponivel - o .exe incluira o motor R ^(landscapemetrics^).
    echo        Com a pasta R_portable presente, o R vai embutido no pacote e
    echo        nao precisa estar instalado na maquina de destino.
)

echo.
echo [5/5] Gerando executavel com PyInstaller (GeoLand.spec)...
pyinstaller GeoLand.spec --noconfirm
if errorlevel 1 (
    echo [ERRO] Falha ao gerar o executavel.
    exit /b 1
)

echo.
echo ============================================================
echo  Build concluido. Executavel disponivel em: dist\GeoLand_V1\GeoLand_V1.exe
echo  Distribua a pasta dist\GeoLand_V1 INTEIRA - o executavel depende dos
echo  arquivos auxiliares nela (GDAL/PROJ, R portatil, etc.).
echo ============================================================

endlocal
