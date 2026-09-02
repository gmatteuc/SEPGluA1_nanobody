# One-time setup for landmark_refine: a private Python environment with CPU torch
# and kornia, the LoFTR weights downloaded, and a self-test.
#
#   cd D:\sep_histology\code
#   .\setup_landmark_refine.ps1
#
# Re-running is safe. Needs a Python 3.10+ on PATH (Anaconda's is fine) and
# internet access for the ~300 MB of packages and the 45 MB of weights.

$ErrorActionPreference = 'Stop'
$here  = Split-Path -Parent $MyInvocation.MyCommand.Path
$pydir = Join-Path $here 'landmark_refine'
$venv  = Join-Path $pydir '.venv'
$py    = Join-Path $venv 'Scripts\python.exe'

if (-not (Test-Path $py)) {
    Write-Host "creating venv at $venv"
    python -m venv $venv
}
& $py -m pip install --quiet --upgrade pip
& $py -m pip install --quiet --index-url https://download.pytorch.org/whl/cpu torch
& $py -m pip install --quiet -r (Join-Path $pydir 'requirements.txt')

Write-Host "downloading LoFTR weights and running a self-test..."
$test = @'
import numpy as np, sys
sys.path.insert(0, r'__PYDIR__')
import core
rng = np.random.default_rng(0)
a = rng.random((200, 280)).astype(np.float32)
core.match_multiscale(a, a, sigmas=(0,))          # forces the weight download
print('landmark_refine: ok, LoFTR loaded')
'@ -replace '__PYDIR__', $pydir
& $py -c $test

Write-Host "done. MATLAB will find this interpreter automatically via landmark_refine.m"
