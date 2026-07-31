#!/usr/bin/env bash
set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

info()  { printf "${GREEN}[INFO]${NC}  %s\n" "$*"; }
warn()  { printf "${YELLOW}[WARN]${NC}  %s\n" "$*"; }
error() { printf "${RED}[ERROR]${NC} %s\n" "$*" >&2; }
say()   { printf '%b\n' "$*"; }

ask() {
    printf "\n"
    for line in "$@"; do printf "  %s\n" "${line}"; done
    printf "  [y/N] "
    read -r REPLY
    case "${REPLY}" in
        [Yy]|[Yy][Ee][Ss]) return 0 ;;
        *) return 1 ;;
    esac
}

if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    error "Do not source this script. Run it directly: bash build.sh"
    return 1
fi

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${THIS_DIR}/environment.yml"
PFAM_SCRIPT="${THIS_DIR}/pfamA_loader.sh"
ENV_NAME="FlaGs2"

if [[ "$(uname)" == "Darwin" ]]; then
    unset DYLD_LIBRARY_PATH || true
fi


info "Checking for Conda..."

if ! command -v conda &>/dev/null; then
    error "Conda executable not found in PATH."
    error "Please install Miniconda or Anaconda and re-run this script."
    error "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

info "Found: $(conda --version 2>&1)"


if [[ ! -f "${ENV_FILE}" ]]; then
    error "Environment file not found: ${ENV_FILE}"
    error "Expected it next to this script: ${THIS_DIR}/environment.yml"
    exit 1
fi

info "Using environment file: ${ENV_FILE}"
info "Detected platform: $(uname -s) / $(uname -m)"


if conda env list | grep -qE "^${ENV_NAME}[[:space:]]"; then
    warn "Environment '${ENV_NAME}' already exists."
    if ask "Remove and recreate it?"; then
        info "Removing existing environment '${ENV_NAME}'..."
        conda env remove --name "${ENV_NAME}" --yes
    else
        info "Keeping existing environment. Skipping creation."
        info "Activate it with:  conda activate ${ENV_NAME}"
        exit 0
    fi
fi


info "Creating Conda environment '${ENV_NAME}' — this may take a few minutes..."

if ! conda env create --name "${ENV_NAME}" --file "${ENV_FILE}"; then
    error "Conda environment creation failed."
    error "Possible fixes:"
    error "  • Check your internet connection."
    error "  • Ensure the bioconda channel is reachable."
    error "  • Run:  conda clean --all  and retry."
    exit 1
fi

info "Environment '${ENV_NAME}' created successfully."


info "Verifying the installation..."

if ! conda run --name "${ENV_NAME}" python -c "import Bio, requests, pyhmmer" &>/dev/null; then
    error "Post-install check failed: core Python packages not importable in '${ENV_NAME}'."
    error "The environment may be incomplete. Remove and re-run:"
    error "  conda env remove --name ${ENV_NAME} && bash build.sh"
    exit 1
fi
info "Core packages OK (Bio, requests, pyhmmer) — default pipeline is ready."

for tool in mafft VeryFastTree; do
    if conda run --name "${ENV_NAME}" command -v "${tool}" &>/dev/null; then
        info "Found tree tool: ${tool}"
    else
        warn "Tree tool '${tool}' not found — --tree / --tree_order will be skipped."
    fi
done

if conda run --name "${ENV_NAME}" python -c "import biolib" &>/dev/null; then
    info "Found pybiolib — --tmhmm / --signalp available."
else
    warn "pybiolib not importable — --tmhmm / --signalp will be skipped."
    warn "  Install with:  conda run --name ${ENV_NAME} pip install pybiolib"
fi

if conda run --name "${ENV_NAME}" python -c "import sismis" &>/dev/null; then
    info "Found sismis — --sismis available."
elif ask "Install sismis now? It is only needed for secretion-system" \
          "detection (the --sismis option) and pulls in gecco," \
          "scikit-learn, scipy, polars and pyrodigal."; then
    info "Installing sismis..."
    if ! conda run --no-capture-output --name "${ENV_NAME}" pip install sismis; then
        error "sismis installation failed. Check the output above."
        error "You can retry at any time:"
        error "  conda run --name ${ENV_NAME} pip install sismis"
        exit 1
    fi
    if ! conda run --name "${ENV_NAME}" python -c "import sismis, pyhmmer" &>/dev/null; then
        error "sismis installed but the environment is now inconsistent"
        error "(sismis or pyhmmer no longer imports). Remove and re-run:"
        error "  conda env remove --name ${ENV_NAME} && bash build.sh"
        exit 1
    fi
    info "sismis installed successfully."
else
    info "Skipping sismis. Install it later with:"
    info "  conda run --name ${ENV_NAME} pip install sismis"
fi

PFAM_HMM="${THIS_DIR}/pfam_db/Pfam-A.hmm"

if [[ -f "${PFAM_HMM}" && -f "${PFAM_HMM}.h3m" ]]; then
    info "Found an indexed Pfam-A database — --domains available."
    info "  --domains --hmmdb pfam_db/Pfam-A.hmm"
elif ask "Install the Pfam-A database now? It is only needed for domain" \
         "annotation (the --domains option). Download is large (~1.5 GB)."; then
    if [[ ! -f "${PFAM_SCRIPT}" ]]; then
        error "pfamA_loader.sh not found at: ${PFAM_SCRIPT}"
        error "Place pfamA_loader.sh next to build.sh and retry."
        exit 1
    fi
    if [[ ! -x "${PFAM_SCRIPT}" ]]; then
        info "Making pfamA_loader.sh executable..."
        chmod +x "${PFAM_SCRIPT}"
    fi
    info "Running pfamA_loader.sh..."
    if ! conda run --no-capture-output --name "${ENV_NAME}" bash "${PFAM_SCRIPT}"; then
        error "pfamA_loader.sh encountered an error. Check the output above."
        error "You can re-run it manually at any time:"
        error "  conda run --name ${ENV_NAME} bash pfamA_loader.sh"
        exit 1
    fi
    info "Pfam-A database installed successfully."
else
    info "Skipping Pfam-A installation."
    info "You can install it later by running:"
    info "  conda run --name ${ENV_NAME} bash pfamA_loader.sh"
fi


printf "\n"
info "Installation complete."
printf "\n"
printf "  Activate the environment:\n"
say "    ${GREEN}conda activate ${ENV_NAME}${NC}"
printf "\n"
printf "  Run FlaGs2 (default — neighbours figure + data tables):\n"
say "    ${GREEN}python FlaGs2.py -i input.txt -u you@example.com -O myrun${NC}"
printf "\n"
printf "  Input list — one query per line, either form, freely mixed:\n"
printf "    WP_047256880.1                            protein only, genome found via NCBI\n"
say "    WP_047256880.1${YELLOW}<TAB>${NC}GCF_000001765.3      protein + NCBI assembly"
say "    MGYG000454827_00001${YELLOW}<TAB>${NC}MGYG000454827   protein + MGnify genome"
say "    ${YELLOW}MGnify genomes need the exact locus tag from that genome's annotation.${NC}"
printf "\n"
printf "  Figures and analysis:\n"
say "    Tree:             ${GREEN}--tree${NC} / ${GREEN}--tree_order${NC}   (mafft + VeryFastTree)"
say "    Domains:          ${GREEN}--domains --hmmdb pfam_db/Pfam-A.hmm${NC}"
say "    Clan colouring:   ${GREEN}--clans pfam_db/Pfam-A.clans.tsv.gz${NC}"
say "    Membrane/signal:  ${GREEN}--tmhmm${NC} / ${GREEN}--signalp${NC}       (pybiolib, uploads sequences)"
say "    Secretion:        ${GREEN}--sismis${NC}                  (sismis)"
say "    Cluster RNAs:     ${GREEN}--cluster_rna${NC}"
printf "\n"
printf "  Other useful options:\n"
say "    Local genomes:    ${GREEN}--use_local DIR${NC}   search a .gff/.faa directory first"
say "    NCBI API key:     ${GREEN}-api KEY${NC}          raises the download rate cap"
say "    Progress:         ${GREEN}-vb${NC}               per-stage progress and timings"
printf "\n"
say "  Deactivate:         ${GREEN}conda deactivate${NC}"
printf "\n"

exit 0