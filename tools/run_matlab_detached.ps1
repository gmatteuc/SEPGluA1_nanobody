# Runs a MATLAB pipeline script in a process that survives the calling session
# closing, and appends everything it prints to a log.
#
# Used for the long unattended stages: copying raw data off the share, running
# extraction over a whole cohort, and so on.
#
#   powershell -File run_matlab_detached.ps1 -Script P0_copy_raw_data
#
# The runner itself lives in the repo so it does not get lost when the log and
# scratch files in the data folder are cleared out.

param(
    [Parameter(Mandatory = $true)][string]$Script,
    [string]$CodeDir = "D:\sep_histology\code",
    [string]$LogDir  = "D:\sep_histology\data\young",
    [string]$MatlabExe = "C:\Program Files\MATLAB\R2024b\bin\matlab.exe"
)

$log = Join-Path $LogDir ("_{0}.log" -f $Script)

"" | Out-File -FilePath $log -Encoding utf8 -Append
"=== $Script started $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ===" |
    Out-File -FilePath $log -Encoding utf8 -Append

& $MatlabExe -batch "cd('$CodeDir'); $Script" 2>&1 |
    Out-File -FilePath $log -Encoding utf8 -Append

"=== $Script finished $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') exit=$LASTEXITCODE ===" |
    Out-File -FilePath $log -Encoding utf8 -Append
