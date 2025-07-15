#!/bin/bash

# exit when any command fails
set -e

MANAGERS=(dnf apt vcpkg brew)
MANAGER=""
LIST=0
VERBOSE=0
DRYRUN=0
SUDO=""

PKGS_OPTIONS=()
PACKAGES=()

# for apt install mkl
APT_LIST_FILE="/etc/apt/sources.list.d/oneAPI.list"
ONEAIP_KEYRINGS_URL="https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB"
ONEAIP_KEYRINGS_FILE="/usr/share/keyrings/oneapi-archive-keyring.gpg"

# Parse Command line
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -v|--verbose)
      VERBOSE=1
      shift
      ;;
    -d|--dry-run)
      DRYRUN=1
      shift
      ;;
    -l|--list)
      LIST=1
      shift
      ;;
    -m|--package-manager)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        MANAGER=($2)
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -h|--help)
      echo "$0 [-vh] [-m package-manager-list]"
      echo "  -m, --package-manager:     preferred package manager order (default: \"${MANAGERS[*]}\")"
      echo "  -v, --verbose:             verbose output"
      echo "  -d, --dry-run:             print actions, but do not execute"
      echo "  -l, --list:                just list the packages to install"
      echo "  -h, --help:                this help message"
      exit 0
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$1"
      shift
      ;;
  esac
done

# Make lower case
PARAMS=$(echo "$PARAMS" | tr '[:upper:]' '[:lower:]' | tr -s "[:blank:]")

# Find an available package manager from the preferred list
# if one has not already been selected manually.
if [ -z "$MANAGER" ]
then
  for m in ${MANAGERS[@]}
  do
      if [ -x "$(command -v $m)" ]; then
          MANAGER="$m"
          break
      fi
  done
fi

# If no package manager is found, exit
if [ -z "$MANAGER" ]
then
      echo "Error: No preferred package managers from list [${MANAGERS[*]}] found. Use -m to select manually." >&2
      exit 1
fi
if ((VERBOSE > 0)); then echo "Using \"$MANAGER\" package manager (select another using -m)"; fi

if [[ "$MANAGER" == "apt" ]]; then
    SUDO="sudo"
    PKGS_OPTIONS+=(install --no-install-suggests --no-install-recommends)
    if ((DRYRUN > 0));  then PKGS_OPTIONS+=(--dry-run); SUDO=""; fi
    PACKAGES+=(build-essential cmake intel-oneapi-mkl-devel libboost-all-dev libgsl-dev libeigen3-dev zlib1g-dev libzstd-dev libsqlite3-dev)
elif [[ "$MANAGER" == "vcpkg" ]]; then
    PKGS_OPTIONS+=(install --triplet=x64-windows)
    if ((DRYRUN > 0));  then PKGS_OPTIONS+=(--dry-run); fi
    PACKAGES+=(gsl eigen3 zlib zstd sqlite3 boost-algorithm boost-lexical-cast boost-random boost-math boost-crc)
elif [[ "$MANAGER" == "dnf" ]]; then
    # TODO:
    SUDO="sudo"
    PKGS_OPTIONS+=(install)
    PACKAGES+=()
    if ((DRYRUN > 0));  then
        MANAGER="echo $MANAGER"
        SUDO=""
    fi
elif [[ "$MANAGER" == "brew" ]]; then
    PKGS_OPTIONS+=(install)
    if ((VERBOSE > 0)); then PKGS_OPTIONS+=(--verbose); fi
    PACKAGES+=(cmake llvm openblas boost gsl eigen zlib zstd sqlite )
    # Brew doesn't have a dryrun option
    if ((DRYRUN > 0));  then
        MANAGER="echo $MANAGER"
    fi
else
    echo "Error: Don't know how to use \"$MANAGER\", please fix the script." >&2
    exit 1
fi

if ((LIST > 0)); then
    echo "${PACKAGES[*]}"
    exit 0
fi

if ((VERBOSE > 0)); then echo "Requesting install of: ${PACKAGES[*]}"; fi



# Add oneAPI apt repository if using apt
if [[ "$MANAGER" == "apt" ]]; then
  if [ ! -f "$APT_LIST_FILE" ]; then
    wget -O- $ONEAIP_KEYRINGS_URL | gpg --dearmor | sudo tee "$ONEAIP_KEYRINGS_FILE"
    echo "deb [signed-by=$ONEAIP_KEYRINGS_FILE] https://apt.repos.intel.com/oneapi all main" | sudo tee "$APT_LIST_FILE" > /dev/null
    sudo apt update
  fi
fi

# Install
$SUDO $MANAGER ${PKGS_OPTIONS[*]} ${PACKAGES[*]}
