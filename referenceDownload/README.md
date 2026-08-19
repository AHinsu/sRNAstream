# RNAcentral Reference Download Guide

A **detailed, beginner-friendly** step-by-step guide for downloading and
extracting RNAcentral reference sequences filtered by **species (NCBI taxid)**
and **RNA type** into separate FASTA files.

> **Quick summary of the approach:**
> 1. Query the **RNAcentral public PostgreSQL database** to get a list of
>    sequence IDs matching your species + RNA type.
> 2. **Download** the large RNAcentral species-specific FASTA file from FTP.
> 3. **Extract** only the sequences you want locally using `seqkit grep`.
>
> This avoids having to host the entire RNAcentral database yourself.

---

## Table of Contents

1. [Folder structure](#folder-structure)
2. [Overview of the approach](#overview-of-the-approach)
3. [Prerequisites and installation](#prerequisites-and-installation)
   - [Windows (native PowerShell)](#windows-native-powershell)
   - [Linux / Ubuntu (including WSL)](#linux--ubuntu-including-wsl)
4. [RNAcentral connection parameters](#rnacentral-connection-parameters)
5. [Quick start: Bos taurus (taxid 9913)](#quick-start-bos-taurus-taxid-9913)
   - [Step 1 – List available RNA types](#step-1--list-available-rna-types)
   - [Step 2 – Generate IDs files for each type](#step-2--generate-ids-files-for-each-type)
   - [Step 3 – Download the RNAcentral FASTA](#step-3--download-the-rnacentral-fasta)
   - [Step 4 – Extract per-type FASTA files](#step-4--extract-per-type-fasta-files)
6. [Multiple RNA types in one run](#multiple-rna-types-in-one-run)
7. [Sequence metadata table](#sequence-metadata-table)
8. [Adapting to other species](#adapting-to-other-species)
9. [Adapting to other RNA types](#adapting-to-other-rna-types)
10. [Credentials and the .env file](#credentials-and-the-env-file)
11. [Troubleshooting](#troubleshooting)
12. [Script reference](#script-reference)

---

## Folder structure

```
referenceDownload/
├── README.md                                      ← this file
├── .gitignore                                     ← ignores large downloads and .env
├── .env.example                                   ← template for DB credentials
├── sql/
│   ├── list_types_for_taxid.sql                   ← lists RNA types for a taxid
│   ├── ids_for_taxid_and_type.sql                 ← gets ID list for one or more types
│   └── sequence_metadata_for_taxid_and_type.sql   ← gets metadata table for sequences
└── scripts/
    ├── run_query_linux.sh         ← bash: run SQL query, save IDs + metadata
    ├── run_query_windows.ps1      ← PowerShell: same
    ├── extract_fasta.sh           ← bash: extract sequences with seqkit
    └── extract_fasta.ps1          ← PowerShell: same
```

---

## Overview of the approach

RNAcentral aggregates non-coding RNA sequences from dozens of expert databases
(Ensembl, RefSeq, Rfam, miRBase, etc.) and assigns each unique sequence a
stable **URS** identifier. Species-specific identifiers combine a URS with an
NCBI taxid, e.g. `URS0000000001_9913` for Bos taurus.

**Why this approach?**

| Option | Pros | Cons |
|--------|------|------|
| RNAcentral public Postgres + FTP FASTA | Fast, handles millions of sequences, fully reproducible | Need psql or Docker |
| REST API | No install needed | Slow for large datasets, pagination limits |
| Full Postgres dump | Complete local copy | Very large (100+ GB), complex to set up |

We use the **public Postgres + FTP FASTA** approach because it is the method
RNAcentral itself recommends on its
[Public database help page](https://rnacentral.org/help/public-database).

---

## Prerequisites and installation

### Windows (native PowerShell)

You will need **one** of the following to run PostgreSQL queries:

- **Option A (recommended): Docker Desktop** — already installed on your system.
  Docker lets you run `psql` inside a temporary container without installing
  PostgreSQL natively.
- **Option B: Native psql** — install PostgreSQL client tools, which include
  the `psql` command-line tool.

You will also need **seqkit** to extract sequences from the FASTA file.

#### Install Docker Desktop (if not already installed)

1. Download Docker Desktop from <https://www.docker.com/products/docker-desktop/>.
2. Run the installer and follow the prompts (Windows Subsystem for Linux 2
   integration is optional — Docker Desktop runs on its own Hyper-V VM by
   default).
3. Start Docker Desktop from the Start menu.
4. Open PowerShell and verify it works:
   ```powershell
   docker run --rm hello-world
   ```
   You should see a "Hello from Docker!" message.

> **WSL note:** If your WSL Ubuntu installation cannot access the Docker daemon
> directly (common when Docker Desktop is set to use the Windows-native Hyper-V
> backend rather than the WSL 2 backend), run all Docker commands from a
> **native Windows PowerShell** window, not from inside WSL. All scripts in
> this folder are designed to work in both environments. See
> [Troubleshooting](#troubleshooting) for more details.

#### Option B: Install native psql (PostgreSQL client only)

1. Go to <https://www.postgresql.org/download/windows/> and download the
   EDB installer.
2. Run the installer. On the "Select Components" screen, tick only
   **Command Line Tools** (you do not need the server itself).
3. Add the PostgreSQL `bin` directory to your `PATH`. The default location is
   something like `C:\Program Files\PostgreSQL\16\bin`. To add it:
   - Open **Settings → System → About → Advanced system settings**.
   - Click **Environment Variables**.
   - Under **User variables**, select `Path` and click **Edit**.
   - Add the path (e.g. `C:\Program Files\PostgreSQL\16\bin`).
4. Open a new PowerShell window and verify:
   ```powershell
   psql --version
   ```

#### Install seqkit (Windows)

1. Go to <https://github.com/shenwei356/seqkit/releases>.
2. Download the latest `seqkit_windows_amd64.exe.tar.gz` (or the `.zip` if
   available).
3. Extract the archive. You will get a file called `seqkit.exe`.
4. Place `seqkit.exe` in a folder on your `PATH`, for example `C:\tools\`.
   (Create that folder if it does not exist and add it to the `Path`
   environment variable as shown in the psql step above.)
5. Verify in PowerShell:
   ```powershell
   seqkit version
   ```

#### Download curl (PowerShell already has Invoke-WebRequest; curl is a built-in alias)

PowerShell's `Invoke-WebRequest` (alias `curl`) works for downloading files.
No extra installation is required.

---

### Linux / Ubuntu (including WSL)

#### Install psql (PostgreSQL client)

```bash
sudo apt update
sudo apt install -y postgresql-client
psql --version
```

#### Install seqkit

**Option 1: via conda/mamba (recommended if you use conda)**

```bash
conda install -c bioconda seqkit
# or
mamba install -c bioconda seqkit
```

**Option 2: via apt (Ubuntu 22.04+)**

```bash
sudo apt install -y seqkit
seqkit version
```

**Option 3: download binary directly**

```bash
# Check https://github.com/shenwei356/seqkit/releases for the latest version
# Replace X.Y.Z with the current release number shown on that page
SEQKIT_VERSION="X.Y.Z"
wget "https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz"
tar -xzf seqkit_linux_amd64.tar.gz
sudo mv seqkit /usr/local/bin/
seqkit version
```

#### Install Docker (optional – only needed if you prefer Docker over native psql)

```bash
# Ubuntu / Debian
sudo apt update
sudo apt install -y docker.io
sudo usermod -aG docker "$USER"
# Log out and back in, then verify:
docker run --rm hello-world
```

> **WSL note:** If you are inside WSL and Docker is installed on the *Windows*
> side only (not the WSL 2 backend), Docker commands will **not** work from
> inside WSL. In that case, either:
> - Use **native psql** inside WSL (recommended — just `apt install
>   postgresql-client`), or
> - Run the PowerShell scripts from a Windows PowerShell window instead.

---

## RNAcentral connection parameters

RNAcentral operates a **free, public PostgreSQL database** that anyone can
query. The connection details are published on their help page:

> **👉 Always obtain the latest credentials from:**
> **<https://rnacentral.org/help/public-database>**

As of writing, the connection parameters are:

| Parameter | Value |
|-----------|-------|
| Host      | `hh-pgsql-public.ebi.ac.uk` |
| Port      | `5432` |
| Database  | `pfmegrnargs` |
| Username  | `reader` |
| Password  | `NWDMCE5xdipIjRrp` |

> **Note:** These are *public, read-only* credentials explicitly published by
> RNAcentral for community use. They are not secret. If the connection fails,
> visit the link above for the current values.

### Storing credentials in a .env file

The wrapper scripts automatically load a `.env` file placed in the
`referenceDownload/` directory. This keeps connection values in one place.
A template is provided:

```
# referenceDownload/.env
# Copy this to .env and fill in current values from:
# https://rnacentral.org/help/public-database

RNACENTRAL_HOST=hh-pgsql-public.ebi.ac.uk
RNACENTRAL_PORT=5432
RNACENTRAL_DB=pfmegrnargs
RNACENTRAL_USER=reader
RNACENTRAL_PASSWORD=NWDMCE5xdipIjRrp
```

The `.env` file is listed in `.gitignore` so it will never be committed to
the repository.

To create it:

```bash
# Linux / WSL
cp referenceDownload/.env.example referenceDownload/.env
# Edit the file if the values have changed since this guide was written:
nano referenceDownload/.env
```

```powershell
# Windows PowerShell
Copy-Item referenceDownload\.env.example referenceDownload\.env
notepad referenceDownload\.env
```

---

## Quick start: Bos taurus (taxid 9913)

Work from the `referenceDownload/` directory for the examples below.

```bash
# Linux / WSL
cd referenceDownload/

# Windows PowerShell
Set-Location referenceDownload\
```

### Step 1 – List available RNA types

Before downloading anything, confirm which RNA types exist for taxid 9913.

**Linux / WSL:**

```bash
bash scripts/run_query_linux.sh --taxid 9913 --list-types
```

**Windows PowerShell:**

```powershell
.\scripts\run_query_windows.ps1 -Taxid 9913 -ListTypes
```

Expected output (truncated example):

```
    rna_type    |   n
----------------+-------
 miRNA          |  1234
 rRNA           |   456
 tRNA           |   789
 snRNA          |    98
 snoRNA         |   321
 lncRNA         |  5432
 ...
(N rows)
```

This confirms the exact spelling of each `rna_type` value. **Use these exact
strings** in the next step.

> **New behaviour:** The RNA-type summary is automatically saved to
> `<taxid>_rna_types.txt` (e.g. `9913_rna_types.txt`) in your working
> directory in addition to being printed to the terminal.  Use `--types-out`
> (bash) or `-TypesOutFile` (PowerShell) to override the output path.

### Step 2 – Generate IDs files for each type

Run one command per RNA type. The script queries the RNAcentral Postgres
database and writes a plain text file with one RNAcentral ID per line.
It also automatically generates a companion **metadata TSV** (see
[Sequence metadata table](#sequence-metadata-table)).

**Linux / WSL:**

```bash
bash scripts/run_query_linux.sh --taxid 9913 --type miRNA --out 9913_miRNA.ids.txt
bash scripts/run_query_linux.sh --taxid 9913 --type tRNA  --out 9913_tRNA.ids.txt
bash scripts/run_query_linux.sh --taxid 9913 --type rRNA  --out 9913_rRNA.ids.txt
```

Or loop over all types at once:

```bash
for TYPE in miRNA tRNA rRNA snRNA snoRNA; do
    bash scripts/run_query_linux.sh \
        --taxid 9913 \
        --type  "$TYPE" \
        --out   "9913_${TYPE}.ids.txt"
done
```

**Windows PowerShell:**

```powershell
.\scripts\run_query_windows.ps1 -Taxid 9913 -RnaType miRNA -OutFile 9913_miRNA.ids.txt
.\scripts\run_query_windows.ps1 -Taxid 9913 -RnaType tRNA  -OutFile 9913_tRNA.ids.txt
.\scripts\run_query_windows.ps1 -Taxid 9913 -RnaType rRNA  -OutFile 9913_rRNA.ids.txt
```

Or as a loop:

```powershell
foreach ($type in @("miRNA","tRNA","rRNA","snRNA","snoRNA")) {
    .\scripts\run_query_windows.ps1 `
        -Taxid   9913 `
        -RnaType $type `
        -OutFile "9913_${type}.ids.txt"
}
```

Each output file looks like:

```
URS0000000001_9913
URS0000000042_9913
URS0000001234_9913
...
```

> **Metadata file:** Each run also writes a tab-separated metadata file
> named `<taxid>_<type>.metadata.tsv` (e.g. `9913_miRNA.metadata.tsv`) by
> default.  Add `--no-metadata` (bash) or `-NoMetadata` (PowerShell) to skip
> it.  See [Sequence metadata table](#sequence-metadata-table) for details.

### Step 3 – Download the RNAcentral FASTA

RNAcentral publishes a **species-specific IDs FASTA** file on their FTP
server. This file contains sequences for every species and is identified by
`URS_TAXID` in its headers — which matches the IDs you generated above.

**FTP documentation:** <https://rnacentral.org/help/ftp>

The file is located at:

```
ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/<release>/sequences/rnacentral_species_specific_ids.fasta.gz
```

To find the latest release number, visit:
<https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/>

> **Important for reproducibility:** Pin the release number so your results
> match a specific version of RNAcentral. Replace `<release>` with the version
> number, e.g. `25.0`.

**Linux / WSL (using wget or curl):**

```bash
# Set the release you want to use (check FTP for latest):
RELEASE="25.0"

wget "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/${RELEASE}/sequences/rnacentral_species_specific_ids.fasta.gz"
```

Or with curl:

```bash
RELEASE="25.0"
curl -O "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/${RELEASE}/sequences/rnacentral_species_specific_ids.fasta.gz"
```

**Windows PowerShell:**

```powershell
$release = "25.0"
$url = "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/$release/sequences/rnacentral_species_specific_ids.fasta.gz"
Invoke-WebRequest -Uri $url -OutFile "rnacentral_species_specific_ids.fasta.gz"
```

> **File size warning:** This file is several gigabytes compressed. Ensure you
> have at least 20 GB of free disk space. The download can take tens of minutes
> on a typical home connection.

### Step 4 – Extract per-type FASTA files

Now use `seqkit grep` to pull out only the sequences whose IDs appear in each
IDs file.

**Linux / WSL:**

```bash
FASTA="rnacentral_species_specific_ids.fasta.gz"

bash scripts/extract_fasta.sh --fasta "$FASTA" --ids 9913_miRNA.ids.txt --out bos_taurus_9913_miRNA.fasta
bash scripts/extract_fasta.sh --fasta "$FASTA" --ids 9913_tRNA.ids.txt  --out bos_taurus_9913_tRNA.fasta
bash scripts/extract_fasta.sh --fasta "$FASTA" --ids 9913_rRNA.ids.txt  --out bos_taurus_9913_rRNA.fasta
```

Or as a loop:

```bash
FASTA="rnacentral_species_specific_ids.fasta.gz"
for TYPE in miRNA tRNA rRNA snRNA snoRNA; do
    bash scripts/extract_fasta.sh \
        --fasta "$FASTA" \
        --ids   "9913_${TYPE}.ids.txt" \
        --out   "bos_taurus_9913_${TYPE}.fasta"
done
```

**Windows PowerShell:**

```powershell
$fasta = "rnacentral_species_specific_ids.fasta.gz"

.\scripts\extract_fasta.ps1 -FastaFile $fasta -IdsFile 9913_miRNA.ids.txt -OutFile bos_taurus_9913_miRNA.fasta
.\scripts\extract_fasta.ps1 -FastaFile $fasta -IdsFile 9913_tRNA.ids.txt  -OutFile bos_taurus_9913_tRNA.fasta
.\scripts\extract_fasta.ps1 -FastaFile $fasta -IdsFile 9913_rRNA.ids.txt  -OutFile bos_taurus_9913_rRNA.fasta
```

Or as a loop:

```powershell
$fasta = "rnacentral_species_specific_ids.fasta.gz"
foreach ($type in @("miRNA","tRNA","rRNA","snRNA","snoRNA")) {
    .\scripts\extract_fasta.ps1 `
        -FastaFile $fasta `
        -IdsFile   "9913_${type}.ids.txt" `
        -OutFile   "bos_taurus_9913_${type}.fasta"
}
```

Each output FASTA file will contain only the sequences for that RNA type in
Bos taurus, ready to use as a reference in tools like miRDeep2, STAR, etc.

---

## Multiple RNA types in one run

You can fetch IDs for **several RNA types at once** and have them combined
into a single output file.  Pass `--type` multiple times **or** supply a
comma-separated list in a single `--type` value.

**Linux / WSL – multiple `--type` flags:**

```bash
bash scripts/run_query_linux.sh \
    --taxid 9913 \
    --type miRNA \
    --type tRNA  \
    --type rRNA
# Output: 9913_miRNA_tRNA_rRNA.ids.txt
#         9913_miRNA_tRNA_rRNA.metadata.tsv
```

**Linux / WSL – comma-separated (equivalent):**

```bash
bash scripts/run_query_linux.sh \
    --taxid 9913 \
    --type  miRNA,tRNA,rRNA
```

**Windows PowerShell – comma-separated string:**

```powershell
.\scripts\run_query_windows.ps1 -Taxid 9913 -RnaType "miRNA,tRNA,rRNA"
# Output: 9913_miRNA_tRNA_rRNA.ids.txt
#         9913_miRNA_tRNA_rRNA.metadata.tsv
```

The resulting IDs file can be passed directly to `extract_fasta.sh` /
`extract_fasta.ps1` to obtain a combined FASTA containing all requested
types, which is useful when building a multi-type reference index.

> **Tip:** When running multiple types, the default output file names include
> all type labels joined by underscores (e.g. `9913_miRNA_tRNA_rRNA.ids.txt`).
> Override with `--out` / `-OutFile` to use a custom path.

---

## Sequence metadata table

Every time IDs are fetched (Step 2 above), the scripts also run a second
query that retrieves **descriptive metadata** for each sequence and writes it
as a tab-separated (TSV) file.  This companion table is intended to help you
track and manage the sequences you have downloaded.

### Columns in the metadata TSV

| Column | Description |
|--------|-------------|
| `rnacentral_id` | Species-specific RNAcentral identifier (`URS…_<taxid>`) |
| `rna_type` | RNA type (miRNA, tRNA, rRNA, …) |
| `description` | Human-readable sequence description from RNAcentral |
| `gene` | Gene symbol from the source annotation (if available) |
| `gene_synonym` | Alternate gene symbols (if available) |
| `product` | RNA product name from the source annotation |
| `locus_tag` | Locus tag from the source annotation |
| `optional_id` | Optional external identifier supplied by the source database |
| `source_database` | Name of the contributing expert database (miRBase, Ensembl, …) |
| `external_accession` | Accession of the sequence in the source database |

> **Note:** Number of exons information is stored in the
> `rnc_genome_coordinates` table of the RNAcentral schema and is not
> included in the default metadata query to keep it fast.  Contact the
> RNAcentral team or consult their schema documentation if you need
> exon-level detail.

### Default file names

| Scenario | IDs file | Metadata file |
|----------|----------|---------------|
| Single type (miRNA) | `9913_miRNA.ids.txt` | `9913_miRNA.metadata.tsv` |
| Multiple types (miRNA + tRNA) | `9913_miRNA_tRNA.ids.txt` | `9913_miRNA_tRNA.metadata.tsv` |

Override either path with `--out` / `--meta-out` (bash) or `-OutFile` /
`-MetaOutFile` (PowerShell).

### Skipping the metadata query

If you do not need the metadata file (e.g. you just want the IDs quickly),
pass `--no-metadata` (bash) or `-NoMetadata` (PowerShell):

```bash
bash scripts/run_query_linux.sh --taxid 9913 --type miRNA --no-metadata
```

```powershell
.\scripts\run_query_windows.ps1 -Taxid 9913 -RnaType miRNA -NoMetadata
```

---

## Adapting to other species

All you need to change is the **NCBI taxonomy ID**. Look up taxids at:
<https://www.ncbi.nlm.nih.gov/taxonomy>

| Species | NCBI taxid |
|---------|-----------|
| Bos taurus (cattle) | 9913 |
| Homo sapiens (human) | 9606 |
| Mus musculus (mouse) | 10090 |
| Gallus gallus (chicken) | 9031 |
| Sus scrofa (pig) | 9823 |
| Ovis aries (sheep) | 9940 |

**Linux example for human (taxid 9606):**

```bash
# 1. List types
bash scripts/run_query_linux.sh --taxid 9606 --list-types

# 2. Get miRNA IDs
bash scripts/run_query_linux.sh --taxid 9606 --type miRNA --out 9606_miRNA.ids.txt

# 3. Extract (same FASTA file – it contains all species)
bash scripts/extract_fasta.sh \
    --fasta rnacentral_species_specific_ids.fasta.gz \
    --ids   9606_miRNA.ids.txt \
    --out   human_9606_miRNA.fasta
```

**PowerShell example for chicken (taxid 9031):**

```powershell
.\scripts\run_query_windows.ps1 -Taxid 9031 -ListTypes
.\scripts\run_query_windows.ps1 -Taxid 9031 -RnaType miRNA -OutFile 9031_miRNA.ids.txt
.\scripts\extract_fasta.ps1 `
    -FastaFile rnacentral_species_specific_ids.fasta.gz `
    -IdsFile   9031_miRNA.ids.txt `
    -OutFile   chicken_9031_miRNA.fasta
```

> **One FASTA file serves all species.** You only need to download
> `rnacentral_species_specific_ids.fasta.gz` once per release. Then run
> different ID queries for different species and extract from the same file.

---

## Adapting to other RNA types

Run `--list-types` first to see which RNA type labels are available for your
taxid. Common values include:

| RNAcentral rna_type | Description |
|---------------------|-------------|
| `miRNA`  | micro RNA |
| `tRNA`   | transfer RNA |
| `rRNA`   | ribosomal RNA |
| `snRNA`  | small nuclear RNA |
| `snoRNA` | small nucleolar RNA |
| `lncRNA` | long non-coding RNA |
| `scaRNA` | small Cajal body-specific RNA |
| `piRNA`  | Piwi-interacting RNA |
| `precursor_RNA` | precursor sequences (e.g. pre-miRNA) |
| `other`  | RNA types not falling into the above categories |

Replace any type name in the commands above. Spelling is **case-sensitive**.

---

## Credentials and the .env file

The scripts read credentials from environment variables. The easiest way to
set them is via a `.env` file.

1. Copy the example:

   ```bash
   # Linux
   cp .env.example .env
   ```

   ```powershell
   # Windows
   Copy-Item .env.example .env
   ```

2. Edit `.env` to match the current values from
   <https://rnacentral.org/help/public-database>.

3. The `.env` file is listed in `.gitignore` and will **never** be committed.

You can also export environment variables directly (useful in CI/CD):

```bash
export RNACENTRAL_HOST=hh-pgsql-public.ebi.ac.uk
export RNACENTRAL_PORT=5432
export RNACENTRAL_DB=pfmegrnargs
export RNACENTRAL_USER=reader
export RNACENTRAL_PASSWORD=NWDMCE5xdipIjRrp
```

```powershell
$env:RNACENTRAL_HOST     = "hh-pgsql-public.ebi.ac.uk"
$env:RNACENTRAL_PORT     = "5432"
$env:RNACENTRAL_DB       = "pfmegrnargs"
$env:RNACENTRAL_USER     = "reader"
$env:RNACENTRAL_PASSWORD = "NWDMCE5xdipIjRrp"
```

---

## Troubleshooting

### SSL / TLS connection errors to PostgreSQL

RNAcentral's public database requires SSL (`sslmode=require`). If you see an
error like `SSL connection error` or `certificate verify failed`:

- **Check your system clock.** TLS certificates are time-sensitive. If your
  system clock is significantly off, certificate verification can fail.
- **Try disabling certificate verification** (not recommended for production,
  but useful for testing):
  ```bash
  PGSSLMODE=verify-ca psql ...   # stricter
  PGSSLMODE=require  psql ...    # what we use (verifies encryption only)
  ```
- **Corporate proxies / firewalls** sometimes block outbound PostgreSQL
  connections (port 5432). Ask your IT department to allow outbound TCP to
  `hh-pgsql-public.ebi.ac.uk:5432`.

### Firewall blocking port 5432

Many corporate or university networks block non-standard ports. Options:

1. Use a personal/home network.
2. Use a VPN that allows port 5432.
3. Ask your network administrator to whitelist `hh-pgsql-public.ebi.ac.uk:5432`.

### Docker not available inside WSL

If you get `Cannot connect to the Docker daemon` inside WSL, this means Docker
Desktop is using the Windows Hyper-V backend and is not bridged into WSL.

**Solutions:**

1. **(Recommended) Install psql natively inside WSL** — this is the simplest:
   ```bash
   sudo apt install -y postgresql-client
   ```
   Then run `run_query_linux.sh` normally. It will use the native `psql`.

2. **Run Docker commands from Windows PowerShell** — open a native PowerShell
   window (not WSL) and use `run_query_windows.ps1` instead.

3. **Enable WSL 2 integration in Docker Desktop:**
   - Open Docker Desktop → Settings → Resources → WSL Integration.
   - Enable integration for your WSL distro.
   - Restart Docker Desktop and your WSL terminal.

### Windows line endings (CRLF) in IDs files

`seqkit grep` reads the IDs file line-by-line. If the file has Windows-style
line endings (`\r\n` — CRLF), trailing `\r` characters get appended to each ID
and no sequences are found.

The scripts attempt to write Unix line endings (`\n` — LF only), but if you
edit an IDs file with Notepad or another Windows editor, line endings may be
converted.

Fix with `dos2unix` (Linux):

```bash
sudo apt install -y dos2unix
dos2unix 9913_miRNA.ids.txt
```

Or in PowerShell (convert to Unix line endings):

```powershell
(Get-Content 9913_miRNA.ids.txt) | Set-Content -NoNewline -Path 9913_miRNA.ids.txt
```

Or with sed (Linux):

```bash
sed -i 's/\r//' 9913_miRNA.ids.txt
```

### seqkit finds 0 sequences

Possible causes:

1. **Wrong FASTA file:** Make sure you downloaded
   `rnacentral_species_specific_ids.fasta.gz`, not the plain
   `rnacentral.fasta.gz` (which uses URS-only headers without the taxid
   suffix).
2. **Release mismatch:** The IDs in your `.ids.txt` file came from a different
   release than the FASTA you downloaded. Regenerate the IDs file or download
   the matching FASTA release.
3. **Line endings in the IDs file** — see above.
4. **Empty IDs file:** The query returned no results. Verify by running
   `--list-types` and confirming the type exists.

### psql: could not translate host name / name resolution failure

Your machine cannot resolve the hostname `hh-pgsql-public.ebi.ac.uk`. This
can happen in WSL when DNS is not configured correctly.

Fix (WSL):

```bash
# Check current DNS
cat /etc/resolv.conf

# If it points to 127.0.0.53 or is empty, try adding Google DNS:
echo "nameserver 8.8.8.8" | sudo tee /etc/resolv.conf
```

For a permanent fix on WSL 2, create `/etc/wsl.conf`:

```ini
[network]
generateResolvConf = false
```

Then set `/etc/resolv.conf` to your preferred DNS server and restart WSL.

### Performance tips

- **Compress the IDs file:** Not necessary — seqkit handles plain text IDs
  efficiently even for millions of entries.
- **Run multiple `seqkit grep` in parallel:** Safe because each reads the same
  compressed FASTA in parallel (seqkit is thread-safe for reading).
- **Keep the FASTA local:** Don't stream it from a network path. Download it
  once and reuse it for all extractions.
- **Use `--threads` / `-j` with seqkit** for faster extraction on multi-core
  machines:
  ```bash
  seqkit grep -j 4 -f ids.txt rnacentral_species_specific_ids.fasta.gz > out.fasta
  ```

---

## Script reference

| Script | Platform | Purpose |
|--------|----------|---------|
| `scripts/run_query_linux.sh` | Linux / WSL Bash | Run SQL queries against RNAcentral Postgres; save IDs and metadata to files |
| `scripts/run_query_windows.ps1` | Windows PowerShell | Same as above |
| `scripts/extract_fasta.sh` | Linux / WSL Bash | Extract sequences using seqkit |
| `scripts/extract_fasta.ps1` | Windows PowerShell | Same as above |
| `sql/list_types_for_taxid.sql` | Any (psql) | List RNA types + counts for a taxid |
| `sql/ids_for_taxid_and_type.sql` | Any (psql) | Get IDs for a taxid + one or more RNA types |
| `sql/sequence_metadata_for_taxid_and_type.sql` | Any (psql) | Get metadata table for a taxid + one or more RNA types |

### Key options for `run_query_linux.sh`

| Option | Description |
|--------|-------------|
| `--taxid NUM` | NCBI taxonomy ID (required) |
| `--type STRING` | RNA type (required unless `--list-types`). Repeat or comma-separate for multiple types |
| `--out FILE` | Output IDs file (default: `<taxid>_<type>.ids.txt`) |
| `--meta-out FILE` | Output metadata TSV (default: `<taxid>_<type>.metadata.tsv`) |
| `--no-metadata` | Skip generating the metadata TSV |
| `--list-types` | List RNA types; saves summary to `<taxid>_rna_types.txt` |
| `--types-out FILE` | Override path for `--list-types` output |
| `--use-docker` | Force Docker even if local psql is available |

### Key options for `run_query_windows.ps1`

| Option | Description |
|--------|-------------|
| `-Taxid INT` | NCBI taxonomy ID (required) |
| `-RnaType STRING` | RNA type(s); comma-separated for multiple types (required unless `-ListTypes`) |
| `-OutFile PATH` | Output IDs file |
| `-MetaOutFile PATH` | Output metadata TSV |
| `-NoMetadata` | Skip generating the metadata TSV |
| `-ListTypes` | List RNA types; saves summary to `<taxid>_rna_types.txt` |
| `-TypesOutFile PATH` | Override path for `-ListTypes` output |
| `-UseDocker` | Force Docker even if local psql is available |

All bash scripts accept `--help` to print usage.

---

## Further resources

- RNAcentral Public Database: <https://rnacentral.org/help/public-database>
- RNAcentral FTP: <https://rnacentral.org/help/ftp>
- seqkit documentation: <https://bioinf.shenwei.me/seqkit/>
- NCBI Taxonomy Browser: <https://www.ncbi.nlm.nih.gov/taxonomy>
- RNAcentral REST API (for small queries): <https://rnacentral.org/api>
