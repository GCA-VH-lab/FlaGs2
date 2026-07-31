#!/usr/bin/env bash

set -e

PFAM_DIR="./pfam_db"
FTP_URL="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"

echo "Environment Check"

if ! python3 -c "import pyhmmer" &> /dev/null; then
    echo "ERROR: pyhmmer is not available to python3."
    echo "pyhmmer is a FlaGs2 dependency and is used here to index (press) the"
    echo "HMM database. Install it into the environment you run FlaGs2 with:"
    echo "  - conda:  conda install -c bioconda pyhmmer"
    echo "  - pip:    pip install pyhmmer"
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

fetch() {
    local name="$1"
    if [[ -s "${name}" ]]; then
        echo "  ${name} already present, skipping download"
    else
        echo "  downloading ${name}"
        $DOWNLOAD_CMD "${FTP_URL}/${name}"
    fi
}

echo "Downloading Pfam files"
if [[ -s "Pfam-A.hmm" ]]; then
    echo "  Pfam-A.hmm already extracted, skipping download"
else
    fetch "Pfam-A.hmm.gz"
    echo "Extracting the HMM database"
    gunzip -f Pfam-A.hmm.gz
fi
fetch "Pfam-A.clans.tsv.gz"

echo "Indexing Pfam-A database with pyhmmer"

rm -f Pfam-A.hmm.h3f Pfam-A.hmm.h3i Pfam-A.hmm.h3m Pfam-A.hmm.h3p

python3 - "Pfam-A.hmm" << 'PYEOF'
import sys, pyhmmer
path = sys.argv[1]
with pyhmmer.plan7.HMMFile(path) as hf:
    n = pyhmmer.hmmpress(hf, path)
print("   pressed {} HMMs".format(n))
PYEOF

echo "Verification"
ls -lh Pfam-A.hmm* Pfam-A.clans.tsv.gz

echo "DONE"
echo
echo "Use with FlaGs2:"
echo "  --domains --hmmdb ${PFAM_DIR}/Pfam-A.hmm"
echo "  --clans ${PFAM_DIR}/Pfam-A.clans.tsv.gz (optional clan colouring)"