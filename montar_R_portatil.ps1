<#
Monta a pasta R_portable (R + landscapemetrics + terra) que o GeoLand.spec
embute no executavel, para o motor R funcionar sem R instalado.

Requer R 4.x (64 bits) instalado nesta maquina e internet (CRAN).

Uso:
    powershell -File montar_R_portatil.ps1
    powershell -File montar_R_portatil.ps1 -RHome "C:\Program Files\R\R-4.5.2"
#>
param(
    [string]$RHome,
    [string]$Destino = (Join-Path $PSScriptRoot 'R_portable')
)
$ErrorActionPreference = 'Stop'

# Sem -RHome, usa o R que estiver no PATH.
if (-not $RHome) {
    $rscript = Get-Command Rscript -ErrorAction SilentlyContinue
    if ($rscript) {
        $RHome = Split-Path $rscript.Source
        while ($RHome -and -not (Test-Path (Join-Path $RHome 'bin\x64\R.dll'))) {
            $RHome = Split-Path $RHome
        }
    }
}
if (-not $RHome -or -not (Test-Path (Join-Path $RHome 'bin\x64\R.dll'))) {
    throw "R (64 bits) nao encontrado. Informe -RHome, ex.: -RHome 'C:\Program Files\R\R-4.5.2'"
}

# 1) R base, sem documentacao, testes, fontes, Tcl/Tk e desinstalador.
robocopy $RHome $Destino /E /NFL /NDL /NJH /NJS /NP /XD doc tests src include Tcl /XF unins000.dat unins000.exe | Out-Null
if ($LASTEXITCODE -ge 8) { throw "Falha ao copiar o R (robocopy $LASTEXITCODE)." }

# 2) Pacotes instalados somente na biblioteca do R portatil.
$env:R_LIBS_USER = Join-Path $Destino 'library'
$env:R_LIBS = $null
$env:R_LIBS_SITE = $null
$codigo = "install.packages(c('landscapemetrics','terra'), lib=.Library, repos='https://cloud.r-project.org', type='binary', quiet=TRUE); ok <- all(sapply(c('landscapemetrics','terra','tibble'), requireNamespace, quietly=TRUE)); quit(status=if (ok) 0 else 1)"
& (Join-Path $Destino 'bin\x64\Rscript.exe') --vanilla -e $codigo
if ($LASTEXITCODE -ne 0) { throw 'Falha ao instalar landscapemetrics/terra no R portatil.' }

Write-Host "R portatil pronto em: $Destino"
