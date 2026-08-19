# =============================================================================
# extract_fasta.ps1
# -----------------------------------------------------------------------------
# Extract sequences from an RNAcentral FASTA file for a given RNA type using
# a pre-generated IDs file, using the seqkit tool (Windows binary).
#
# Requirements:
#   seqkit for Windows – download the latest release binary from:
#     https://github.com/shenwei356/seqkit/releases
#   Place seqkit.exe somewhere on your PATH (e.g. C:\tools\seqkit.exe) or
#   supply the full path via -SeqkitPath.
#
# Usage:
#   .\extract_fasta.ps1 `
#       -FastaFile rnacentral_species_specific_ids.fasta.gz `
#       -IdsFile   9913_miRNA.ids.txt `
#       -OutFile   bos_taurus_9913_miRNA.fasta
#
#   # Loop over multiple types:
#   foreach ($type in @("miRNA","tRNA","rRNA","snRNA","snoRNA")) {
#       .\extract_fasta.ps1 `
#           -FastaFile rnacentral_species_specific_ids.fasta.gz `
#           -IdsFile   "9913_${type}.ids.txt" `
#           -OutFile   "bos_taurus_9913_${type}.fasta"
#   }
# =============================================================================

[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$FastaFile,

    [Parameter(Mandatory = $true)]
    [string]$IdsFile,

    [string]$OutFile = "output.fasta",

    [string]$SeqkitPath = "seqkit"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
if (-not (Test-Path $FastaFile)) {
    Write-Error "FASTA file not found: $FastaFile"
    exit 1
}
if (-not (Test-Path $IdsFile)) {
    Write-Error "IDs file not found: $IdsFile"
    exit 1
}

# ---------------------------------------------------------------------------
# Check seqkit
# ---------------------------------------------------------------------------
$seqkitCmd = Get-Command $SeqkitPath -ErrorAction SilentlyContinue
if (-not $seqkitCmd) {
    Write-Error @"
seqkit not found at '$SeqkitPath'.
Download the Windows binary from:
  https://github.com/shenwei356/seqkit/releases
Unzip seqkit.exe and either:
  - Place it in a folder that is on your PATH, or
  - Pass its full path via -SeqkitPath, e.g.:
      -SeqkitPath 'C:\tools\seqkit.exe'
"@
    exit 1
}

# ---------------------------------------------------------------------------
# Count IDs
# ---------------------------------------------------------------------------
$idLines = Get-Content $IdsFile | Where-Object { $_ -match '\S' }
$idCount = $idLines.Count
Write-Host "==> Extracting $idCount sequences from $FastaFile ..."
Write-Host "    IDs file : $IdsFile"
Write-Host "    Output   : $OutFile"

# ---------------------------------------------------------------------------
# Run seqkit grep
# ---------------------------------------------------------------------------
& $SeqkitPath grep -f $IdsFile $FastaFile | Set-Content -Path $OutFile -Encoding UTF8

# ---------------------------------------------------------------------------
# Verify output
# ---------------------------------------------------------------------------
$outCount = (Select-String -Path $OutFile -Pattern '^>' -SimpleMatch).Count
Write-Host "==> Done. $outCount sequences written to: $OutFile"

if ($idCount -ne $outCount) {
    Write-Warning "Expected $idCount sequences but found $outCount."
    Write-Warning "Some IDs may not be present in the FASTA file."
    Write-Warning "Check that you are using the correct FASTA release."
}
