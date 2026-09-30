# One-time setup for auto_annotation: a private Python environment with torch
# (the CUDA build when an NVIDIA GPU is present), the packages, and a self-test
# on the bundled weights.
#
#   cd D:\sep_histology\code
#   .\registration\auto_annotation\setup.ps1
#
# Re-running is safe. Needs a Python 3.10+ on PATH (Anaconda's is fine) and
# internet access for the packages (~2.5 GB with CUDA torch).

$ErrorActionPreference = 'Stop'
# this script sits in the engine's folder, next to requirements.txt and weights\
$pydir = Split-Path -Parent $MyInvocation.MyCommand.Path
$venv  = Join-Path $pydir '.venv'
$py    = Join-Path $venv 'Scripts\python.exe'

if (-not (Test-Path $py)) {
    Write-Host "creating venv at $venv"
    python -m venv $venv
}
& $py -m pip install --quiet --upgrade pip

# The registration runs ~150 optimisation steps per section; on a GPU a whole
# brain takes minutes, on the CPU much longer. The CPU build still works.
$gpu = $false
try { & nvidia-smi -L 2>$null | Out-Null; $gpu = ($LASTEXITCODE -eq 0) } catch { $gpu = $false }
if ($gpu) {
    Write-Host "NVIDIA GPU found: installing the CUDA build of torch (about 2.5 GB)"
    & $py -m pip install --quiet --index-url https://download.pytorch.org/whl/cu124 torch==2.6.0
} else {
    Write-Host "no NVIDIA GPU: installing the CPU build of torch (a brain will take a while)"
    & $py -m pip install --quiet --index-url https://download.pytorch.org/whl/cpu torch==2.6.0
}
& $py -m pip install --quiet -r (Join-Path $pydir 'requirements.txt')

foreach ($w in 'landmark.pt', 'matcher.pt') {
    if (-not (Test-Path (Join-Path $pydir "weights\$w"))) { throw "weights\$w is missing" }
}

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
