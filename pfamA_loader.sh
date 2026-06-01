#!/usr/bin/env bash

set -e

PFAM_DIR="./pfam_db"
FTP_URL="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"

echo "Environment Check"

if ! command -v hmmpress &> /dev/null; then
    echo "ERROR: 'hmmpress' could not be found."
    echo "Please install HMMER before running this script."
    echo "  - macOS: brew install hmmer"
    echo "  - Linux (Ubuntu/Debian): sudo apt-get install hmmer"
    echo "  - Linux (Arch): yay -S hmmer"
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
$DOWNLOAD_CMD "${FTP_URL}/Pfam-A.hmm.dat.gz"
$DOWNLOAD_CMD "${FTP_URL}/active_site.dat.gz"

echo "Extracting compressed files"
gunzip -f Pfam-A.hmm.gz
gunzip -f Pfam-A.hmm.dat.gz
gunzip -f active_site.dat.gz

echo "Assembling and indexing Pfam-A database for HMMER"
hmmpress -f Pfam-A.hmm

echo "Verification"
ls -lh Pfam-A.hmm*

echo "DONE"
