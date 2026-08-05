# Sets $OVERLEAF: the folder the .tex tables are copied into.
# Source this after cd'ing to the Filtering folder:  . ./overleaf_dir.sh
#
# Same policy as overleaf_dir.m. Overleaf is the destination ONLY on a machine
# that already has the authors' sync folder; everyone else gets Filtering/output
# so a replicator never writes into a Dropbox path they don't have.
#
#   1. $SCRAMBLING_QUANTFIGS if set (created if absent)
#   2. $HOME/Dropbox/Apps/Overleaf/<project>/quantfigs, only if it exists
#   3. <this folder>/output

_ovl_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -n "${SCRAMBLING_QUANTFIGS:-}" ]; then
    OVERLEAF="$SCRAMBLING_QUANTFIGS"
else
    OVERLEAF="$HOME/Dropbox/Apps/Overleaf/ScramblingDollarsLiquidity_NewVersion_Restud/quantfigs"
    [ -d "$OVERLEAF" ] || OVERLEAF="$_ovl_here/output"
fi

mkdir -p "$OVERLEAF"
export OVERLEAF
unset _ovl_here
