#!/usr/bin/env bash
# Rend les scripts exécutables et propose de les ajouter au PATH.
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
chmod +x "$DIR"/*.py "$DIR"/*.sh
echo "Scripts rendus exécutables."
echo "Ajoute ceci à ton ~/.bashrc pour les appeler de partout :"
echo "  export PATH=\"\$PATH:$DIR/scripts:$DIR/bash\""
