# Setup of the automatic annotation's engine, once per machine: a Python
# environment of its own (.venv, beside this script) with torch, the CUDA build
# when an NVIDIA GPU is present, the other packages pinned in
# tools\requirements_auto_annotation.txt, and a self-test that loads the bundled
# weights.
#
#   cd D:\sep_histology\code
#   .\registration\auto_annotation\setup.ps1
#
# Running it again is safe. Needs Python 3.12 on PATH (the environment was frozen
# with Anaconda's 3.12.7) and internet access for the packages (about 2.5 GB with
# CUDA torch).

$ErrorActionPreference = 'Stop'

# this script sits in the engine's folder, next to weights\, two levels below the
# code root, whose tools\ holds the pinned packages
$pydir = Split-Path -Parent $MyInvocation.MyCommand.Path
$code_root = Split-Path -Parent (Split-Path -Parent $pydir)
$pinned = Join-Path $code_root 'tools\requirements_auto_annotation.txt'
$venv  = Join-Path $pydir '.venv'
$py    = Join-Path $venv 'Scripts\python.exe'

if (-not (Test-Path $py)) {
    Write-Host "creating venv at $venv"
    python -m venv $venv
}
& $py -m pip install --quiet --upgrade pip

# the CUDA build when nvidia-smi finds a GPU: the registration runs about 500
# optimisation steps per batch of sections, minutes a brain on a GPU and much
# longer on the CPU, whose build still works
$gpu = $false
try { & nvidia-smi -L 2>$null | Out-Null; $gpu = ($LASTEXITCODE -eq 0) } catch { $gpu = $false }
if ($gpu) {
    Write-Host "NVIDIA GPU found: installing the CUDA build of torch (about 2.5 GB)"
    & $py -m pip install --quiet --index-url https://download.pytorch.org/whl/cu124 torch==2.6.0
} else {
    Write-Host "no NVIDIA GPU: installing the CPU build of torch (a brain will take a while)"
    & $py -m pip install --quiet --index-url https://download.pytorch.org/whl/cpu torch==2.6.0
}

# the other pinned packages, without the torch line: its +cu124 build would fail on
# a machine without CUDA, and torch is in place already
$packages = Get-Content $pinned | Where-Object {
    $_ -match '^[A-Za-z0-9_.-]+==' -and $_ -notmatch '^torch=='
}
& $py -m pip install --quiet @packages

# the weights are tracked with the code; stop if either is missing
foreach ($w in 'landmark.pt', 'matcher.pt') {
    if (-not (Test-Path (Join-Path $pydir "weights\$w"))) { throw "weights\$w is missing" }
}

# self-test: load both models with this environment's torch
Write-Host "running a self-test..."
$test = @'
import numpy as np, os, sys
sys.path.insert(0, r'__PYDIR__')
import core
w = os.path.join(r'__PYDIR__', 'weights')
lnet, mnet = core.load_models(os.path.join(w, 'landmark.pt'), os.path.join(w, 'matcher.pt'))
print('auto_annotation: ok, models loaded on', core.device())
'@ -replace '__PYDIR__', $pydir
& $py -c $test

Write-Host "done. MATLAB finds this interpreter automatically via auto_annotate.m"
