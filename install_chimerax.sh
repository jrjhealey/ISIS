#!/bin/bash
#
# ISIS ChimeraX plugin installer.
#
# Installs both halves - the isis-epitope prediction library into ChimeraX's own
# Python, and the ChimeraX-ISIS bundle into ChimeraX - then VERIFIES the result
# by running `isis doctor` and reporting what ChimeraX actually sees.
#
# That verification step is the point. The previous version of this script could
# fail at every stage and still print "Installation complete!", because it piped
# pip through `while read` and then tested $? - which is the exit status of the
# while loop, not of pip. Silent failure is the whole reason an install appears
# to "do nothing".

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUNDLE_DIR="$SCRIPT_DIR/src/isis_chimerax"

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC} $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; }

CHIMERAX_PYTHON=""
CHIMERAX_BIN=""
CHIMERAX_SITE_PACKAGES=""

find_chimerax() {
    # An explicit override is the only reliable answer on an unusual layout.
    if [[ -n "${CHIMERAX_PYTHON_OVERRIDE:-}" ]]; then
        CHIMERAX_PYTHON="$CHIMERAX_PYTHON_OVERRIDE"
        CHIMERAX_BIN="${CHIMERAX_BIN_OVERRIDE:-$(command -v chimerax || true)}"
        return
    fi

    local prefixes=()
    if [[ "$OSTYPE" == darwin* ]]; then
        # Newest app bundle first, so a current install wins over an old one.
        while IFS= read -r app; do
            prefixes+=("$app/Contents")
        done < <(ls -d /Applications/ChimeraX*.app 2>/dev/null | sort -rV)
    elif [[ "$OSTYPE" == linux* ]]; then
        # The Debian/Ubuntu .deb installs to /usr/lib/ucsf-chimerax, which the
        # previous version of this script never looked for - so on a stock
        # apt/dpkg install it simply failed to find ChimeraX.
        while IFS= read -r d; do
            prefixes+=("$d")
        done < <(ls -d /usr/lib/ucsf-chimerax* /opt/UCSF/ChimeraX* /usr/lib/chimerax* \
                       /opt/chimerax* /usr/local/chimerax* "$HOME"/.local/chimerax* \
                       2>/dev/null | sort -rV)
    else
        while IFS= read -r d; do
            prefixes+=("$d")
        done < <(ls -d "/c/Program Files/ChimeraX"* "/mnt/c/Program Files/ChimeraX"* \
                       2>/dev/null | sort -rV)
    fi

    for prefix in "${prefixes[@]:-}"; do
        [[ -d "$prefix" ]] || continue
        # Glob for the real interpreter rather than guessing a version, and skip
        # python*-config, which is executable and easy to pick by mistake.
        for py in "$prefix"/bin/python3.[0-9] "$prefix"/bin/python3.[0-9][0-9] \
                  "$prefix"/bin/python3 "$prefix"/bin/python.exe; do
            [[ -x "$py" && "$py" != *-config ]] || continue
            if "$py" -c 'import sys' >/dev/null 2>&1; then
                CHIMERAX_PYTHON="$py"
                break
            fi
        done
        if [[ -n "$CHIMERAX_PYTHON" ]]; then
            for b in "$prefix/MacOS/ChimeraX" "$prefix/bin/ChimeraX" \
                     "$prefix/bin/chimerax" "$prefix/bin/ChimeraX.exe"; do
                [[ -x "$b" ]] && { CHIMERAX_BIN="$b"; break; }
            done
            [[ -z "$CHIMERAX_BIN" ]] && CHIMERAX_BIN="$(command -v chimerax || true)"
            return
        fi
    done
}

main() {
    echo
    echo "=================================="
    echo "  ISIS ChimeraX plugin installer"
    echo "=================================="
    echo

    info "Looking for ChimeraX..."
    find_chimerax

    # Ask ChimeraX itself for the interpreter and site-packages it actually uses,
    # and prefer that over anything guessed from the directory layout.
    #
    # This is the fix for the worst variant of "the installer ran and nothing
    # changed". On macOS, globbing the app bundle finds
    #   .../ChimeraX.app/Contents/bin/python3.11
    # but ChimeraX runs out of
    #   .../Contents/Library/Frameworks/Python.framework/Versions/3.11/...
    # If those do not share a site-packages, pip installs into a directory
    # ChimeraX never imports from: pip reports success, and nothing changes.
    # ChimeraX's own `pip` command cannot be used instead - it only accepts
    # package names, not paths or wheels - so the interpreter has to be right.
    if [[ -n "${CHIMERAX_BIN:-}" && -x "${CHIMERAX_BIN:-}" ]]; then
        local probe="${TMPDIR:-/tmp}/isis_probe_$$.py"
        cat > "$probe" <<'PYPROBE'
import os, sys, sysconfig
def run(session):
    py = os.path.join(sys.prefix, "bin",
                      "python%d.%d" % sys.version_info[:2])
    if not os.path.exists(py):
        py = os.path.join(sys.prefix, "bin", "python3")
    print("ISIS_PY:" + (py if os.path.exists(py) else ""))
    print("ISIS_SP:" + (sysconfig.get_paths().get("purelib") or ""))
run(session)
PYPROBE
        local probe_out authoritative_py
        probe_out="$("$CHIMERAX_BIN" --nogui --exit --script "$probe" 2>/dev/null)"
        rm -f "$probe"
        authoritative_py="$(sed -n 's/^ISIS_PY://p' <<< "$probe_out" | head -1)"
        CHIMERAX_SITE_PACKAGES="$(sed -n 's/^ISIS_SP://p' <<< "$probe_out" | head -1)"
        if [[ -n "$authoritative_py" && -x "$authoritative_py" ]]; then
            if [[ "$authoritative_py" != "$CHIMERAX_PYTHON" ]]; then
                warn "Detected python was $CHIMERAX_PYTHON"
                warn "but ChimeraX itself runs $authoritative_py - using ChimeraX's."
            fi
            CHIMERAX_PYTHON="$authoritative_py"
        fi
    fi

    if [[ -z "$CHIMERAX_PYTHON" || ! -x "$CHIMERAX_PYTHON" ]]; then
        error "Could not find ChimeraX's Python interpreter."
        echo
        echo "Set it explicitly and re-run, for example:"
        echo "  CHIMERAX_PYTHON_OVERRIDE=/Applications/ChimeraX.app/Contents/bin/python3.11 $0"
        echo "  CHIMERAX_PYTHON_OVERRIDE=/usr/lib/ucsf-chimerax/bin/python3.11 $0"
        exit 1
    fi

    info "ChimeraX Python : $CHIMERAX_PYTHON"
    info "ChimeraX binary : ${CHIMERAX_BIN:-<not found>}"
    [[ -n "$CHIMERAX_SITE_PACKAGES" ]] && \
        info "Site-packages   : $CHIMERAX_SITE_PACKAGES"
    # Confirm the chosen interpreter writes where ChimeraX reads. If it does not,
    # every later step would succeed while changing nothing ChimeraX can see.
    if [[ -n "$CHIMERAX_SITE_PACKAGES" ]]; then
        local pip_target
        pip_target="$("$CHIMERAX_PYTHON" -c \
            'import sysconfig;print(sysconfig.get_paths().get("purelib",""))' 2>/dev/null)"
        if [[ -n "$pip_target" && "$pip_target" != "$CHIMERAX_SITE_PACKAGES" ]]; then
            error "The interpreter found does NOT install where ChimeraX reads:"
            error "  interpreter installs to : $pip_target"
            error "  ChimeraX imports from   : $CHIMERAX_SITE_PACKAGES"
            error "Anything installed with it would be invisible to ChimeraX."
            error "Set CHIMERAX_PYTHON_OVERRIDE to an interpreter whose"
            error "site-packages matches, and re-run."
            exit 1
        fi
    fi
    echo

    # ---- 1. the prediction library, WITH the extras it needs -----------------
    # [ml] pulls scikit-learn and joblib, without which the MHC models cannot be
    # loaded; [plot] pulls matplotlib for `isis plot`. Installing bare - as this
    # script used to - yields a working-looking install with T-cell prediction
    # unusable.
    info "Installing the isis-epitope library into ChimeraX's Python..."
    if ! "$CHIMERAX_PYTHON" -m pip install --upgrade "${SCRIPT_DIR}[ml,plot]"; then
        error "Failed to install the isis-epitope library."
        error "ChimeraX would load the plugin but every prediction method would"
        error "report itself unavailable. Fix the error above and re-run."
        exit 1
    fi
    info "Library installed."
    echo

    # ---- 2. the ChimeraX bundle ---------------------------------------------
    if [[ ! -d "$BUNDLE_DIR" ]]; then
        error "Bundle directory not found: $BUNDLE_DIR"
        exit 1
    fi

    if [[ -z "$CHIMERAX_BIN" || ! -x "$CHIMERAX_BIN" ]]; then
        warn "Could not find the ChimeraX executable, so the bundle was not installed."
        warn "Run this inside ChimeraX:"
        warn "    devel install \"$BUNDLE_DIR\""
        exit 1
    fi

    # Clear stale build artefacts first: a previously built wheel gets reused and
    # can carry old metadata, silently reinstalling the very version you are
    # trying to replace.
    rm -rf "$BUNDLE_DIR/build" "$BUNDLE_DIR/dist" "$BUNDLE_DIR"/*.egg-info 2>/dev/null || true

    # Remove EVERY existing copy of the bundle first.
    #
    # ChimeraX searches its user data directory (on macOS
    # ~/Library/Application Support/ChimeraX/<ver>/lib/pythonX.Y/site-packages)
    # BEFORE the application directory. A bundle put there by an earlier
    # `devel install` therefore shadows a newer one installed into the app: the
    # new files are on disk but never imported, so new commands come back as
    # "Unknown command" and the listing stays stale, while the installer reports
    # success. Ask ChimeraX itself where the copies are - it knows its own
    # sys.path - then delete them all.
    info "Removing any previously installed copies of the bundle..."
    local finder="${TMPDIR:-/tmp}/isis_find_copies_$$.py"
    cat > "$finder" <<'PYFIND'
import os, sys
def run(session):
    seen = []
    for p in sys.path:
        cand = os.path.join(p, "chimerax", "isis")
        if os.path.isdir(cand):
            real = os.path.realpath(cand)
            if real not in seen:
                seen.append(real)
                print("ISIS_COPY:" + real)
    try:
        from chimerax import app_dirs
        for root, dirs, _f in os.walk(app_dirs.user_data_dir):
            if os.path.basename(root) == "isis" and \
               os.path.basename(os.path.dirname(root)) == "chimerax":
                real = os.path.realpath(root)
                if real not in seen:
                    seen.append(real)
                    print("ISIS_COPY:" + real)
    except Exception:
        pass
run(session)
PYFIND
    local copies
    copies="$("$CHIMERAX_BIN" --nogui --exit --script "$finder" 2>/dev/null \
              | sed -n 's/^ISIS_COPY://p' | sort -u)"
    rm -f "$finder"
    if [[ -n "$copies" ]]; then
        while IFS= read -r c; do
            [[ -n "$c" && -d "$c" ]] || continue
            info "  removing $c"
            rm -rf "$c"
            # take the dist-info/egg-info with it, or pip believes it is present
            rm -rf "$(dirname "$(dirname "$c")")"/chimerax_isis-*.dist-info 2>/dev/null || true
            rm -rf "$(dirname "$(dirname "$c")")"/ChimeraX_ISIS-*.egg-info 2>/dev/null || true
        done <<< "$copies"
    else
        info "  none found"
    fi

    # `devel install` REFUSES to replace a bundle already installed at the same
    # version - it prints "already installed with the same version" and leaves the
    # old code running. That is the single most confusing failure here: the
    # installer reports success and nothing changes. So build the wheel, then
    # install it with `reinstall true`, which does replace it.
    info "Building the ChimeraX bundle..."
    "$CHIMERAX_BIN" --nogui --exit --cmd "devel build \"$BUNDLE_DIR\"" 2>&1 \
        | sed 's/^/    /'

    local wheel
    wheel="$(ls -t "$BUNDLE_DIR"/dist/*.whl 2>/dev/null | head -1)"
    if [[ -z "$wheel" ]]; then
        error "Bundle build produced no wheel in $BUNDLE_DIR/dist"
        exit 1
    fi
    # Install with pip --force-reinstall via ChimeraX's own interpreter.
    #
    # Neither `devel install` nor `toolshed install ... reinstall true` will
    # replace a bundle already present at the same version - both stop at
    # "already installed with the same version", leaving the old code running
    # while reporting success. Verified by planting a marker in the source,
    # rebuilding at an unchanged version, and checking whether it reached the
    # installed bundle: it did not for either of those, and does for this.
    # ChimeraX still registers the bundle correctly afterwards (`toolshed list`
    # shows the new metadata), so nothing is lost by going through pip.
    info "Installing bundle: $(basename "$wheel")"
    if ! "$CHIMERAX_PYTHON" -m pip install --force-reinstall --no-deps "$wheel"; then
        error "Bundle installation failed."
        exit 1
    fi
    echo

    # ---- 3. verify, and say what is actually true ---------------------------
    info "Verifying the installation (isis doctor)..."
    echo
    local report
    report="$("$CHIMERAX_BIN" --nogui --exit --cmd "isis doctor" 2>&1)"
    echo "$report" | grep -vE '^(INFO:|WARNING:|STATUS:)$' | sed 's/^/    /'
    echo

    if echo "$report" | grep -q "copies of the bundle are installed"; then
        error "More than one copy of the bundle is installed and an old one is"
        error "shadowing the new one. The report above lists them. Remove the"
        error "extras by hand and re-run, or delete the ChimeraX user bundle dir."
        exit 1
    fi

    if echo "$report" | grep -q "Everything required is present"; then
        info "Installation verified."
        echo "$report" | sed -n 's/.*Bundle loaded *: *//p' | head -1 \
            | while read -r loaded; do info "Bundle in use: $loaded"; done
        echo
        echo "Try it:"
        echo "    open 1ubq"
        echo "    isis list"
        echo "    isis bcell consensus #1"
        echo "    isis plot #1"
        echo
        exit 0
    fi

    error "Verification FAILED - ChimeraX cannot see a usable installation."
    error "The report above names the problem. Common causes:"
    error "  * ChimeraX was already running: restart it, then run 'isis doctor'."
    error "  * More than one ChimeraX installed. This script used:"
    error "      $CHIMERAX_PYTHON"
    error "    Set CHIMERAX_PYTHON_OVERRIDE to the one you actually launch."
    exit 1
}

uninstall() {
    find_chimerax
    info "Uninstalling ISIS..."
    if [[ -x "$CHIMERAX_PYTHON" ]]; then
        "$CHIMERAX_PYTHON" -m pip uninstall -y isis-epitope || true
    fi
    if [[ -n "${CHIMERAX_BIN:-}" && -x "${CHIMERAX_BIN:-}" ]]; then
        "$CHIMERAX_BIN" --nogui --exit --cmd "toolshed uninstall ChimeraX-ISIS" || true
    fi
    info "Done. Restart ChimeraX."
}

case "${1:-}" in
    --uninstall|-u) uninstall ;;
    --help|-h)
        echo "Usage: $0 [--uninstall] [--help]"
        echo
        echo "Installs the isis-epitope library into ChimeraX's Python and the"
        echo "ChimeraX-ISIS bundle into ChimeraX, then verifies with 'isis doctor'."
        echo
        echo "Environment:"
        echo "  CHIMERAX_PYTHON_OVERRIDE   path to ChimeraX's python (skips detection)"
        echo "  CHIMERAX_BIN_OVERRIDE      path to the ChimeraX executable"
        ;;
    *) main ;;
esac
