#!/usr/bin/env bash
set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

info()  { printf "${GREEN}[INFO]${NC}  %s\n" "$*"; }
warn()  { printf "${YELLOW}[WARN]${NC}  %s\n" "$*"; }
error() { printf "${RED}[ERROR]${NC} %s\n" "$*" >&2; }

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
    printf "  Remove and recreate it? [y/N] "
    read -r REPLY
    case "${REPLY}" in
        [Yy]|[Yy][Ee][Ss])
            info "Removing existing environment '${ENV_NAME}'..."
            conda env remove --name "${ENV_NAME}" --yes
            ;;
        *)
            info "Keeping existing environment. Skipping creation."
            info "Activate it with:  conda activate ${ENV_NAME}"
            exit 0
            ;;
    esac
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
info "Core packages OK (Bio, requests, pyhmmer)."

for tool in mafft VeryFastTree; do
    if conda run --name "${ENV_NAME}" command -v "${tool}" &>/dev/null; then
        info "Found tree tool: ${tool}"
    else
        warn "Tree tool '${tool}' not found in the environment."
        warn "  The default pipeline still works; only --tree / --tree_order need it."
    fi
done


printf "\n"
printf "  Install the Pfam-A database now? It is only needed for domain\n"
printf "  annotation (the --domains option). Download is large (~1.5 GB).\n"
printf "  [y/N] "
read -r PFAM_REPLY

case "${PFAM_REPLY}" in
    [Yy]|[Yy][Ee][Ss])
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
        if ! bash "${PFAM_SCRIPT}"; then
            error "pfamA_loader.sh encountered an error. Check the output above."
            error "You can re-run it manually at any time:  bash pfamA_loader.sh"
            exit 1
        fi
        info "Pfam-A database installed successfully."
        ;;
    *)
        info "Skipping Pfam-A installation."
        info "You can install it later by running:  bash pfamA_loader.sh"
        ;;
esac

printf "\n"
info "Installation complete."
printf "\n"
printf "  Activate the environment:\n"
printf "    ${GREEN}conda activate ${ENV_NAME}${NC}\n"
printf "\n"
printf "  Run FlaGs2 (default — neighbours figure + data tables):\n"
printf "    ${GREEN}python FlaGs2.py -i input.txt -u you@example.com -O myrun${NC}\n"
printf "\n"
printf "  With a tree:        ${GREEN}--tree${NC}        (needs mafft + VeryFastTree)\n"
printf "  With domains:       ${GREEN}--domains --hmmdb pfam_db/Pfam-A.hmm${NC}\n"
printf "  Clan colouring:     ${GREEN}--clans pfam_db/Pfam-A.clans.tsv.gz${NC}\n"
printf "\n"
printf "  Deactivate:         ${GREEN}conda deactivate${NC}\n"
printf "\n"

exit 0
