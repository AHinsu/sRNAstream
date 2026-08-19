# =============================================================================
# run_query_windows.ps1
# -----------------------------------------------------------------------------
# Run an SQL query against the RNAcentral public PostgreSQL database from
# Windows PowerShell.
#
# Supports two backends:
#   1. Docker (Windows Docker Desktop) – recommended if Docker is installed
#   2. Local psql client (if PostgreSQL client tools are installed natively)
#
# Prerequisites:
#   - Docker Desktop for Windows  OR  psql installed (see README)
#   - A .env file in the referenceDownload\ folder, OR set the $env:
#     variables listed below before running this script.
#
# Usage examples:
#   # List RNA types available for Bos taurus (taxid 9913):
#   .\run_query_windows.ps1 -Taxid 9913 -ListTypes
#
#   # Generate IDs file for a single type:
#   .\run_query_windows.ps1 -Taxid 9913 -RnaType miRNA
#
#   # Generate IDs file for multiple types (combined into one file):
#   .\run_query_windows.ps1 -Taxid 9913 -RnaType miRNA,tRNA,rRNA
#
#   # Generate IDs file and specify output path:
#   .\run_query_windows.ps1 -Taxid 9913 -RnaType tRNA -OutFile C:\data\bos_tRNA.ids.txt
#
#   # Skip automatic metadata TSV:
#   .\run_query_windows.ps1 -Taxid 9913 -RnaType miRNA -NoMetadata
#
#   # Force Docker even if local psql is present:
#   .\run_query_windows.ps1 -Taxid 9913 -RnaType rRNA -UseDocker
# =============================================================================

[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [int]$Taxid,

    # One or more RNA types; comma-separated string, e.g. "miRNA,tRNA,rRNA"
    [string]$RnaType = "",

    [string]$OutFile = "",

    # Output path for the metadata TSV (default: <taxid>_<types>.metadata.tsv)
    [string]$MetaOutFile = "",

    # Output path for the --list-types summary (default: <taxid>_rna_types.txt)
    [string]$TypesOutFile = "",

    [switch]$ListTypes,

    # Skip generating the sequence-metadata TSV
    [switch]$NoMetadata,

    [switch]$UseDocker
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# ---------------------------------------------------------------------------
# Resolve paths
# ---------------------------------------------------------------------------
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Definition
$RepoRoot  = Split-Path -Parent $ScriptDir
$SqlDir    = Join-Path $RepoRoot "sql"
$EnvFile   = Join-Path $RepoRoot ".env"

# ---------------------------------------------------------------------------
# Load .env file if it exists (key=value format, one per line)
# ---------------------------------------------------------------------------
if (Test-Path $EnvFile) {
    Get-Content $EnvFile | ForEach-Object {
        if ($_ -match '^\s*([^#][^=]+?)\s*=\s*(.*)\s*$') {
            [System.Environment]::SetEnvironmentVariable($Matches[1].Trim(), $Matches[2].Trim(), 'Process')
        }
    }
    Write-Host "INFO: Loaded connection parameters from $EnvFile"
}

# ---------------------------------------------------------------------------
# Connection parameters
# Override by setting environment variables or entries in .env.
# Current values: https://rnacentral.org/help/public-database
# ---------------------------------------------------------------------------
$DbHost = if ($env:RNACENTRAL_HOST)     { $env:RNACENTRAL_HOST }     else { "hh-pgsql-public.ebi.ac.uk" }
$DbPort = if ($env:RNACENTRAL_PORT)     { $env:RNACENTRAL_PORT }     else { "5432" }
$DbName = if ($env:RNACENTRAL_DB)       { $env:RNACENTRAL_DB }       else { "pfmegrnargs" }
$DbUser = if ($env:RNACENTRAL_USER)     { $env:RNACENTRAL_USER }     else { "reader" }
$DbPass = if ($env:RNACENTRAL_PASSWORD) { $env:RNACENTRAL_PASSWORD } else { "NWDMCE5xdipIjRrp" }

# ---------------------------------------------------------------------------
# Parse RNA types – split comma-separated value into an array
# ---------------------------------------------------------------------------
$RnaTypesArray = @()
if (-not [string]::IsNullOrWhiteSpace($RnaType)) {
    $RnaTypesArray = $RnaType -split ',' | ForEach-Object { $_.Trim() } | Where-Object { $_ -ne "" }
}

# ---------------------------------------------------------------------------
# Validate arguments
# ---------------------------------------------------------------------------
if (-not $ListTypes -and $RnaTypesArray.Count -eq 0) {
    Write-Error "You must supply -RnaType (e.g. miRNA or miRNA,tRNA,rRNA) unless you use -ListTypes."
    exit 1
}

# ---------------------------------------------------------------------------
# Helper: build the SQL IN-list string, e.g. 'miRNA','tRNA','rRNA'
# ---------------------------------------------------------------------------
function Build-TypesIn {
    ($RnaTypesArray | ForEach-Object { "'$_'" }) -join ","
}

# ---------------------------------------------------------------------------
# Build a label string for file names, e.g. miRNA_tRNA_rRNA
# ---------------------------------------------------------------------------
function Types-Label {
    $RnaTypesArray -join "_"
}

# ---------------------------------------------------------------------------
# Default output files
# ---------------------------------------------------------------------------
if ([string]::IsNullOrWhiteSpace($OutFile) -and -not $ListTypes) {
    $OutFile = Join-Path (Get-Location) "${Taxid}_$(Types-Label).ids.txt"
}
if ([string]::IsNullOrWhiteSpace($MetaOutFile) -and -not $ListTypes) {
    $MetaOutFile = Join-Path (Get-Location) "${Taxid}_$(Types-Label).metadata.tsv"
}
if ([string]::IsNullOrWhiteSpace($TypesOutFile)) {
    $TypesOutFile = Join-Path (Get-Location) "${Taxid}_rna_types.txt"
}

# ---------------------------------------------------------------------------
# Determine backend: Docker or local psql
# ---------------------------------------------------------------------------
$UseDockerBackend = $UseDocker.IsPresent

if (-not $UseDockerBackend) {
    $psqlPath = Get-Command psql -ErrorAction SilentlyContinue
    if (-not $psqlPath) {
        $dockerPath = Get-Command docker -ErrorAction SilentlyContinue
        if ($dockerPath) {
            $UseDockerBackend = $true
            Write-Host "INFO: psql not found – falling back to Docker (postgres:16)."
        } else {
            Write-Error @"
Neither 'psql' nor 'docker' was found.
Options:
  1. Install Docker Desktop: https://www.docker.com/products/docker-desktop/
  2. Install PostgreSQL client tools (which include psql):
     https://www.postgresql.org/download/windows/
     During installation, select 'Command Line Tools'.
     Then add the bin directory to your PATH, e.g.:
       C:\Program Files\PostgreSQL\16\bin
"@
            exit 1
        }
    }
}

# ---------------------------------------------------------------------------
# Helper: build and run the SQL command
# ---------------------------------------------------------------------------
$ConnStr = "host=${DbHost} port=${DbPort} dbname=${DbName} user=${DbUser} password=${DbPass} sslmode=require"

function Invoke-PsqlQuery {
    param(
        [string]$SqlContent,
        [bool]$TuplesOnly    = $false,
        [bool]$UseTsvFlags   = $false   # -A -F tab -q  (for metadata TSV)
    )

    # Write SQL to a temp file (avoids quoting nightmares)
    $TmpSql = [System.IO.Path]::GetTempFileName() -replace '\.tmp$', '.sql'
    Set-Content -Path $TmpSql -Value $SqlContent -Encoding UTF8

    try {
        if ($UseDockerBackend) {
            $SqlBytes = [System.Text.Encoding]::UTF8.GetBytes($SqlContent)
            $PsqlArgs = @(
                "run", "--rm", "-i",
                "-e", "PGPASSWORD=$DbPass",
                "postgres:16", "psql",
                $ConnStr
            )
            if ($TuplesOnly)  { $PsqlArgs += @("-t", "-A") }
            if ($UseTsvFlags) { $PsqlArgs += @("-A", "-F", "`t", "-q") }

            $result = $SqlBytes | & docker @PsqlArgs
        } else {
            $env:PGPASSWORD = $DbPass
            $PsqlConnStr = "host=${DbHost} port=${DbPort} dbname=${DbName} user=${DbUser} sslmode=require"
            $PsqlArgs = @("-d", $PsqlConnStr, "-f", $TmpSql)
            if ($TuplesOnly)  { $PsqlArgs += @("-t", "-A") }
            if ($UseTsvFlags) { $PsqlArgs += @("-A", "-F", "`t", "-q") }
            $result = & psql @PsqlArgs
            Remove-Item Env:PGPASSWORD -ErrorAction SilentlyContinue
        }
        return $result
    } finally {
        Remove-Item $TmpSql -Force -ErrorAction SilentlyContinue
    }
}

# ===========================================================================
# --list-types: show RNA types for the taxid, save summary to file
# ===========================================================================
if ($ListTypes) {
    Write-Host "==> Listing RNA types for taxid $Taxid ..."
    $SqlTemplate = Get-Content (Join-Path $SqlDir "list_types_for_taxid.sql") -Raw
    $SqlContent  = $SqlTemplate -replace ':TAXID', $Taxid
    $output = Invoke-PsqlQuery -SqlContent $SqlContent -TuplesOnly $false

    # Print to terminal
    $output | Write-Host

    # Save to file
    [System.IO.File]::WriteAllLines($TypesOutFile, $output, [System.Text.UTF8Encoding]::new($false))
    Write-Host "==> RNA type summary saved to: $TypesOutFile"
    exit 0
}

# ===========================================================================
# Fetch sequence IDs
# ===========================================================================
$TypesIn      = Build-TypesIn
$TypesDisplay = $RnaTypesArray -join ", "

Write-Host "==> Querying IDs for taxid=$Taxid, type(s): $TypesDisplay ..."
$SqlTemplate = Get-Content (Join-Path $SqlDir "ids_for_taxid_and_type.sql") -Raw
$SqlContent  = $SqlTemplate -replace ':TAXID', $Taxid
$SqlContent  = $SqlContent  -replace ':RNA_TYPES', $TypesIn

$output = Invoke-PsqlQuery -SqlContent $SqlContent -TuplesOnly $true

# Write output – ensure Unix line endings (LF only) for seqkit compatibility
$lines = $output | Where-Object { $_ -ne "" }
[System.IO.File]::WriteAllLines($OutFile, $lines, [System.Text.UTF8Encoding]::new($false))

$count = $lines.Count
Write-Host "==> Done. $count IDs written to: $OutFile"

# ===========================================================================
# Fetch sequence metadata TSV – skipped when -NoMetadata is specified
# ===========================================================================
if (-not $NoMetadata) {
    Write-Host "==> Querying metadata for taxid=$Taxid, type(s): $TypesDisplay ..."
    $MetaTemplate = Get-Content (Join-Path $SqlDir "sequence_metadata_for_taxid_and_type.sql") -Raw
    $MetaSql      = $MetaTemplate -replace ':TAXID', $Taxid
    $MetaSql      = $MetaSql      -replace ':RNA_TYPES', $TypesIn

    try {
        $metaOutput = Invoke-PsqlQuery -SqlContent $MetaSql -UseTsvFlags $true
        [System.IO.File]::WriteAllLines($MetaOutFile, $metaOutput, [System.Text.UTF8Encoding]::new($false))
        Write-Host "==> Metadata ($($metaOutput.Count) lines incl. header) written to: $MetaOutFile"
    } catch {
        Write-Warning "Metadata query failed – skipping metadata file. Error: $_"
    }
}
