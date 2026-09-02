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
    # One script, or several separated by commas, run in order.
    # Note this has to be a plain string rather than a string[]: powershell -File
    # binds "a,b" as a single element, so the splitting is done here instead.
    [Parameter(Mandatory = $true)][string]$Script,
    [int]$WaitForPid = 0,          # optional: hold until this process exits first
    [int]$MaxWaitMinutes = 480,
    [string]$CodeDir = "D:\sep_histology\code",
    [string]$LogDir  = "D:\sep_histology\data\young",
    [string]$MatlabExe = "C:\Program Files\MATLAB\R2024b\bin\matlab.exe"
)

# @() around the pipeline: with a single script name the pipeline returns a bare
# string rather than an array, and $scripts[0] then indexes its first character,
# so every single-stage run logged to _P.log instead of _<script>.log.
$scripts = @($Script -split ',' | ForEach-Object { $_.Trim() } | Where-Object { $_ -ne '' })
$log = Join-Path $LogDir ("_{0}.log" -f $scripts[0])

# Chain onto an earlier stage when asked, so a long copy can hand over to
# processing without anyone being at the keyboard.
if ($WaitForPid -gt 0) {
    "waiting for PID $WaitForPid to finish ($(Get-Date -Format 'HH:mm:ss'))" |
        Out-File -FilePath $log -Encoding utf8 -Append
    $waited = 0
    while ((Get-Process -Id $WaitForPid -ErrorAction SilentlyContinue) -and ($waited -lt $MaxWaitMinutes)) {
        Start-Sleep -Seconds 60
        $waited++
    }
    if (Get-Process -Id $WaitForPid -ErrorAction SilentlyContinue) {
        "PID $WaitForPid still running after $MaxWaitMinutes min - not starting $Script" |
            Out-File -FilePath $log -Encoding utf8 -Append
        exit 1
    }
    "PID $WaitForPid finished after ~$waited min" | Out-File -FilePath $log -Encoding utf8 -Append
}

# Run the scripts one after another. If one fails we stop rather than feeding a
# later stage with half-finished input.
foreach ($s in $scripts) {
    "" | Out-File -FilePath $log -Encoding utf8 -Append
    "=== $s started $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ===" |
        Out-File -FilePath $log -Encoding utf8 -Append

    & $MatlabExe -batch "cd('$CodeDir'); $s" 2>&1 |
        Out-File -FilePath $log -Encoding utf8 -Append

    $code = $LASTEXITCODE
    "=== $s finished $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') exit=$code ===" |
        Out-File -FilePath $log -Encoding utf8 -Append

    if ($code -ne 0) {
        "$s did not exit cleanly - stopping, later stages not started" |
            Out-File -FilePath $log -Encoding utf8 -Append
        exit 1
    }
}
