#!/usr/bin/env bash

set -e

PFAM_DIR="./pfam_db"
FTP_URL="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"

echo "Environment Check"

if ! command -v hmmpress &> /dev/null; then
    echo "ERROR: 'hmmpress' could not be found."
    echo "hmmpress indexes the HMM database so domain scans start faster."
    echo "It ships with HMMER:"
    echo "  - conda:                  conda install -c bioconda hmmer"
    echo "  - macOS:                  brew install hmmer"
    echo "  - Linux (Ubuntu/Debian):  sudo apt-get install hmmer"
    echo "  - Linux (Arch):           sudo pacman -S hmmer"
    exit 1
fi

if command -v curl &> /dev/null; then
    DOWNLOAD_CMD="curl -fLO"
elif command -v wget &> /dev/null; then
    DOWNLOAD_CMD="wget -c"
else
    echo "ERROR: Neither 'curl' nor 'wget' was found. Please install one of them."
    exit 1
fi

echo "Creating database directory: ${PFAM_DIR}"
mkdir -p "${PFAM_DIR}"
cd "${PFAM_DIR}"

echo "Downloading Pfam files"
$DOWNLOAD_CMD "${FTP_URL}/Pfam-A.hmm.gz"
$DOWNLOAD_CMD "${FTP_URL}/Pfam-A.clans.tsv.gz"

echo "Extracting the HMM database"
gunzip -f Pfam-A.hmm.gz

echo "Indexing Pfam-A database for HMMER (hmmpress)"
hmmpress -f Pfam-A.hmm

echo "Verification"
ls -lh Pfam-A.hmm* Pfam-A.clans.tsv.gz

echo "DONE"
echo
echo "Use with FlaGs2:"
echo "  --domains --hmmdb ${PFAM_DIR}/Pfam-A.hmm"
echo "  --clans ${PFAM_DIR}/Pfam-A.clans.tsv.gz   (optional clan colouring)"
