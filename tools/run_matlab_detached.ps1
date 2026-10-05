# Runs a MATLAB pipeline script in a process that survives the calling session
# closing, and appends everything it prints to a log.
#
# Used for the long unattended stages: copying raw data off the share, running
# extraction over a whole cohort, and so on, and for the refactor's checks,
# which run the old and the new code on copies of the data
# (docs/history/REFACTOR_PLAN.md, Verification design).
#
#   powershell -NonInteractive -File run_matlab_detached.ps1 -Script run_collect_by_group -CodeDir D:\sep_histology\code -DataRoot D:\sep_histology\data -LogDir D:\sep_histology\data\young
#
# -NonInteractive makes a forgotten parameter an error instead of a prompt
# that a detached launch would wait on forever.
#
# The code folder, the data root and the log folder are always given: a check
# must never fall back on the production folders by default. The data root is
# passed to MATLAB as SEP_DATA_ROOT, which get_paths honours, so the same code
# can run on the production data or on a check tree's copy. Each session
# starts from the default path, adds the code with sep_setup_paths, and prints
# the code root, which get_paths it resolved, the data root and the commit
# before running anything.
#
# A stage can be a script name or a call, such as
# "addpath('tools'); sep_run_driver_copy('run_collect_by_group', 'G:\sep_refactor\check\data')".
# Use single quotes inside it: double quotes do not survive the hand-over to
# matlab.exe.
#
# The runner itself lives in the repo so it does not get lost when the log and
# scratch files in the data folder are cleared out.

param(
    # one script, or several separated by commas, run in order; a plain string, not a
    # string[], which powershell -File would fill with "a,b" as one element
    [Parameter(Mandatory = $true)][string]$Script,
    [Parameter(Mandatory = $true)][string]$CodeDir,
    [Parameter(Mandatory = $true)][string]$DataRoot,
    [Parameter(Mandatory = $true)][string]$LogDir,

    # optional: hold the start until this process has ended
    [int]$WaitForPid = 0,
    [int]$MaxWaitMinutes = 480,
    [string]$MatlabExe = "C:\Program Files\MATLAB\R2024b\bin\matlab.exe"
)

# the folders a check must not touch
$ProductionCode = "D:\sep_histology\code"
$ProductionData = "D:\sep_histology\data"
$SnapshotPrefix = "G:\sep_histology_snapshot"

# Absolute, with backslashes and no trailing separator, so the checks below see one
# spelling of each folder; a relative path is taken from the current location.
function Get-CanonicalPath([string]$Path) {
    $full = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($Path)
    $full = [System.IO.Path]::GetFullPath($full)
    if ($full.Length -gt 3) {
        $full = $full.TrimEnd('\')
    }
    return $full
}

# Whether a path is the folder or inside it, ignoring case as Windows does.
function Test-Inside([string]$Path, [string]$Folder) {
    return ($Path -ieq $Folder) -or
        $Path.StartsWith($Folder + '\', [System.StringComparison]::OrdinalIgnoreCase)
}

# Refusals go to the console only: the log folder may be the reason for them.
function Stop-Refused([string]$Message) {
    Write-Host "refused: $Message"
    exit 1
}

# Split on the commas between stages, not on those inside a call's brackets
# or quotes, so a stage can be a call with several arguments.
function Split-Stages([string]$Text) {
    $stages = @()
    $current = ""
    $depth = 0
    $quoted = $false
    foreach ($c in $Text.ToCharArray()) {
        # a quote opens or closes a text; a bracket outside a text changes the depth
        if ($c -eq "'") {
            $quoted = -not $quoted
        } elseif (-not $quoted -and "([{".Contains([string]$c)) {
            $depth++
        } elseif (-not $quoted -and ")]}".Contains([string]$c)) {
            $depth--
        }

        # a comma outside brackets and texts ends a stage
        if ($c -eq ',' -and $depth -eq 0 -and -not $quoted) {
            $stages += $current
            $current = ""
        } else {
            $current += $c
        }
    }
    $stages += $current
    return @($stages | ForEach-Object { $_.Trim() } | Where-Object { $_ -ne '' })
}

# one spelling of each folder
$CodeDir = Get-CanonicalPath $CodeDir
$DataRoot = Get-CanonicalPath $DataRoot
$LogDir = Get-CanonicalPath $LogDir

# only the production code runs on the production data (get_paths insists too), a
# check's log stays out of it, and the snapshot on G: is never run on or written to
if ((Test-Inside $DataRoot $ProductionData) -and -not ($CodeDir -ieq $ProductionCode)) {
    Stop-Refused "the code in $CodeDir is a copy, but the data root $DataRoot is the production data. Point it at the copy's own check tree."
}
if ($CodeDir.StartsWith($SnapshotPrefix, [System.StringComparison]::OrdinalIgnoreCase)) {
    Stop-Refused "the code folder $CodeDir is inside the snapshot on G:, which is a backup, never run."
}
if ((Test-Inside $LogDir $ProductionData) -and -not ($DataRoot -ieq $ProductionData)) {
    Stop-Refused "the log folder $LogDir is in the production data, but the data root is $DataRoot. Put the log in the check tree."
}
if ($DataRoot.StartsWith($SnapshotPrefix, [System.StringComparison]::OrdinalIgnoreCase)) {
    Stop-Refused "the data root $DataRoot is inside the snapshot on G:, which is a backup."
}
if ($LogDir.StartsWith($SnapshotPrefix, [System.StringComparison]::OrdinalIgnoreCase)) {
    Stop-Refused "the log folder $LogDir is inside the snapshot on G:, which is a backup."
}
if (-not (Test-Path -LiteralPath $DataRoot -PathType Container)) {
    Stop-Refused "the data root $DataRoot does not exist."
}
if ($Script.Contains('"')) {
    Stop-Refused "the script contains a double quote, which matlab.exe would not receive intact. Use single quotes."
}

# the stages; @() keeps a single stage an array, since the function then returns a
# bare string, and $scripts[0] would be its first character
$scripts = @(Split-Stages $Script)
if ($scripts.Count -eq 0) {
    Stop-Refused "no script given."
}
if (-not (Test-Path -LiteralPath $LogDir -PathType Container)) {
    New-Item -ItemType Directory -Path $LogDir | Out-Null
}

# the log is named after the first stage: a run('<file>') stage after its file,
# anything else by its letters, digits, _ and -, at most 60 of them
$logName = $scripts[0]
if ($logName -match "^run\('([^']+)'\)$") {
    $logName = [System.IO.Path]::GetFileNameWithoutExtension($Matches[1])
}
$logName = ($logName -replace '[^\w\-]+', '_').Trim('_')
if ($logName.Length -gt 60) {
    $logName = $logName.Substring(0, 60)
}
$log = Join-Path $LogDir ("_{0}.log" -f $logName)

# the commit the code is at, and its uncommitted changes, so the log says which code
# ran; git refuses a folder on exFAT (the check trees on G:), which records no owner,
# unless it is marked safe: here for these two calls only, never in a git configuration
$safe = "safe.directory=$($CodeDir.Replace('\', '/'))"
$commit = (& git -c $safe -C $CodeDir rev-parse HEAD 2>$null)
if ($LASTEXITCODE -ne 0 -or -not $commit) {
    $commit = "unknown (not a git checkout)"
} else {
    $changed = @(& git -c $safe -C $CodeDir status --porcelain --untracked-files=no 2>$null)
    if ($changed.Count -gt 0) {
        $commit = "$commit with $($changed.Count) uncommitted changes"
    }
}

# the log's header: when, which code, which data, which commit, which stages
"" | Out-File -FilePath $log -Encoding utf8 -Append
"=== run_matlab_detached $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ===" |
    Out-File -FilePath $log -Encoding utf8 -Append
"code     $CodeDir" | Out-File -FilePath $log -Encoding utf8 -Append
"data     $DataRoot" | Out-File -FilePath $log -Encoding utf8 -Append
"commit   $commit" | Out-File -FilePath $log -Encoding utf8 -Append
"scripts  $($scripts -join ' | ')" | Out-File -FilePath $log -Encoding utf8 -Append

# the data root reaches MATLAB through the environment; the previous value is put
# back at the end, for a run from an open PowerShell session
$previousDataRoot = $env:SEP_DATA_ROOT
$env:SEP_DATA_ROOT = $DataRoot
try {
    # wait for an earlier process when asked, so a long copy can hand over to
    # processing unattended; one check a minute, at most MaxWaitMinutes
    if ($WaitForPid -gt 0) {
        "waiting for PID $WaitForPid to finish ($(Get-Date -Format 'HH:mm:ss'))" |
            Out-File -FilePath $log -Encoding utf8 -Append
        $waited = 0
        while ((Get-Process -Id $WaitForPid -ErrorAction SilentlyContinue) -and
                ($waited -lt $MaxWaitMinutes)) {
            Start-Sleep -Seconds 60
            $waited++
        }
        if (Get-Process -Id $WaitForPid -ErrorAction SilentlyContinue) {
            "PID $WaitForPid still running after $MaxWaitMinutes min - not starting $Script" |
                Out-File -FilePath $log -Encoding utf8 -Append
            exit 1
        }
        "PID $WaitForPid finished after ~$waited min" |
            Out-File -FilePath $log -Encoding utf8 -Append
    }

    # every session starts from the default path, so nothing a previous run added
    # can shadow the code under test, and says what it resolved
    $codeDirMatlab = $CodeDir -replace "'", "''"
    $setup = "restoredefaultpath; cd('$codeDirMatlab'); sep_setup_paths; " +
        "fprintf('code root  %s\n', fileparts(which('sep_setup_paths'))); " +
        "fprintf('get_paths  %s\n', which('get_paths')); " +
        "fprintf('data root  %s\n', get_paths().data); " +
        "fprintf('commit     %s\n', '$commit');"

    # run the stages in order; a failing stage stops the chain, so no later stage
    # reads half-finished input
    foreach ($s in $scripts) {
        "" | Out-File -FilePath $log -Encoding utf8 -Append
        "=== $s started $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ===" |
            Out-File -FilePath $log -Encoding utf8 -Append

        # a fresh MATLAB for each stage, everything it prints into the log
        & $MatlabExe -batch "$setup $s" 2>&1 |
            Out-File -FilePath $log -Encoding utf8 -Append

        $code = $LASTEXITCODE
        "=== $s finished $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') exit=$code ===" |
            Out-File -FilePath $log -Encoding utf8 -Append

        # stop the chain on a failure
        if ($code -ne 0) {
            "$s did not exit cleanly - stopping, later stages not started" |
                Out-File -FilePath $log -Encoding utf8 -Append
            exit 1
        }
    }
} finally {
    $env:SEP_DATA_ROOT = $previousDataRoot
}
