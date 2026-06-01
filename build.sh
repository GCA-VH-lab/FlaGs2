#!/usr/bin/env bash

# FlaGs2 — Environment installer
# Supports: Linux (linux-64) and macOS (osx-64)

set -euo pipefail

# Output colorizer
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Colour

info()    { printf "${GREEN}[INFO]${NC}  %s\n" "$*"; }
warn()    { printf "${YELLOW}[WARN]${NC}  %s\n" "$*"; }
error()   { printf "${RED}[ERROR]${NC} %s\n" "$*" >&2; }

# Check for run with bash / sh, not sourced

if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    error "Do not source this script. Run it directly: bash build.sh"
    return 1
fi

# Path resolver 

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_PATH="${THIS_DIR}/env"
MAC_FILE="${ENV_PATH}/eFlaGs2.txt"
LIN_FILE="${ENV_PATH}/elinFlaGs2.txt"
PFAM_SCRIPT="${THIS_DIR}/pfamA_loader.sh"
ENV_NAME="eFlaGs2"


# macOS only: unset DYLD_LIBRARY_PATH to avoid library conflicts

if [[ "$(uname)" == "Darwin" ]]; then
    unset DYLD_LIBRARY_PATH
fi

#  1. Check for Conda

info "Checking for Conda..."

if ! command -v conda &>/dev/null; then
    error "Conda executable not found in PATH."
    error "Please install Miniconda or Anaconda and re-run this script."
    error "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

CONDA_VERSION="$(conda --version 2>&1)"
info "Found: ${CONDA_VERSION}"

#  2. Detect OS and select the right package list

OS="$(uname -s)"

case "${OS}" in
    Linux)
        PLATFORM="linux-64"
        ENV_FILE="${LIN_FILE}"
        ;;
    Darwin)
        PLATFORM="osx-64"
        ENV_FILE="${MAC_FILE}"
        ;;
    *)
        error "Unsupported OS: ${OS}. Only Linux and macOS are supported."
        exit 1
        ;;
esac

info "Detected platform: ${PLATFORM}"

#  3. Check for environment spec file

if [[ ! -f "${ENV_FILE}" ]]; then
    error "Package list not found: ${ENV_FILE}"
    error "Expected directory layout:"
    error "  $(dirname "${ENV_FILE}")/"
    error "    ├── eFlaGs2.txt      (macOS)"
    error "    └── elinFlaGs2.txt   (Linux)"
    exit 1
fi

info "Using package list: ${ENV_FILE}"

#  4. Remove stale environment if it already exists

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
            ;;
    esac
fi

#  5. Create the Conda environment

info "Creating Conda environment '${ENV_NAME}' — this may take a few minutes..."

if ! conda create --name "${ENV_NAME}" --file "${ENV_FILE}" --yes; then
    error "Conda environment creation failed."
    error "Possible fixes:"
    error "  • Check your internet connection."
    error "  • Run:  conda clean --all  and retry."
    error "  • Verify the package URLs in ${ENV_FILE} are still valid."
    exit 1
fi

info "Environment '${ENV_NAME}' created successfully."

#  6. macOS — verify critical binaries landed inside the env

if [[ "${OS}" == "Darwin" ]]; then
    CONDA_ENV_BIN="$(conda run --name "${ENV_NAME}" conda info --base 2>/dev/null || true)"
    # Light smoke-test: check Python is reachable in the new env
    if ! conda run --name "${ENV_NAME}" python --version &>/dev/null; then
        error "Post-install check failed: Python not found in '${ENV_NAME}'."
        error "The environment may be incomplete. Try removing it and re-running:"
        error "  conda env remove --name ${ENV_NAME} && bash build.sh"
        exit 1
    fi
    info "macOS post-install check passed."
fi

#  7. Optional: download Pfam-A database

printf "\n"
printf "  Do you want to install the Pfam-A database for domain annotation?\n"
printf "  (Required for HMMER-based domain annotation — download is ~500 MB)\n"
printf "  [y/N] "
read -r PFAM_REPLY

case "${PFAM_REPLY}" in
    [Yy]|[Yy][Ee][Ss])
        if [[ ! -f "${PFAM_SCRIPT}" ]]; then
            error "pfamA_loader.sh not found at: ${PFAM_SCRIPT}"
            error "Please place pfamA_loader.sh in the same directory as build.sh and retry."
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

#  8. Done — print usage instructions

printf "\n"
info "Installation complete."
printf "\n"
printf "  To activate the FlaGs2 environment:\n"
printf "    ${GREEN}conda activate ${ENV_NAME}${NC}\n"
printf "\n"
printf "  To deactivate:\n"
printf "    ${GREEN}conda deactivate${NC}\n"
printf "\n"

exit 0