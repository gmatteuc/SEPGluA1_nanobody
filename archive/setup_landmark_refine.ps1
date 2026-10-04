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

# An NVIDIA GPU makes the matching ~10x faster; use the CUDA build of torch when
# one is present, the CPU build otherwise. Either works, the wrapper does not care.
$gpu = $false
try { & nvidia-smi -L 2>$null | Out-Null; $gpu = ($LASTEXITCODE -eq 0) } catch { $gpu = $false }
if ($gpu) {
    Write-Host "NVIDIA GPU found: installing the CUDA build of torch (about 2.5 GB)"
    & $py -m pip install --quiet --index-url https://download.pytorch.org/whl/cu124 torch==2.6.0
} else {
    Write-Host "no NVIDIA GPU: installing the CPU build of torch"
    & $py -m pip install --quiet --index-url https://download.pytorch.org/whl/cpu torch==2.6.0
}
& $py -m pip install --quiet -r (Join-Path $pydir 'requirements.txt')

# The LoFTR weights ship with the package (weights\loftr_outdoor.ckpt, 45 MB),
# so no download is needed; core.py only fetches them if that file is missing.
$w = Join-Path $pydir 'weights\loftr_outdoor.ckpt'
if (Test-Path $w) { Write-Host "weights: bundled file found" }
else { Write-Host "weights: bundled file MISSING, they will be downloaded on first use" }

Write-Host "running a self-test..."
$test = @'
import numpy as np, sys
sys.path.insert(0, r'__PYDIR__')
import core
rng = np.random.default_rng(0)
a = rng.random((200, 280)).astype(np.float32)
core.match_multiscale(a, a, sigmas=(0,))
print('landmark_refine: ok, LoFTR loaded on', core.device())
'@ -replace '__PYDIR__', $pydir
& $py -c $test

Write-Host "done. MATLAB will find this interpreter automatically via landmark_refine.m"
