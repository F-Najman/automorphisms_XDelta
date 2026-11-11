#!/bin/bash
# Prints the versions of the programs used (MAGMA can be set, default magma).
MAGMA=${MAGMA:-magma}
echo "Magma $(echo 'print GetVersion(); quit;' | "$MAGMA" -b 2>/dev/null | tr ' ' '.' | tr -d '\n')"
sage --version
python3 --version
gcc --version | head -1
uname -srm
