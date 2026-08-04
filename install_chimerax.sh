#!/bin/bash
#
# ISIS ChimeraX Plugin Installer
# Installs both the core library and ChimeraX bundle
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
echo_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
echo_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Find ChimeraX installation
find_chimerax() {
    local chimerax_python=""
    local chimerax_app=""

    # macOS locations
    if [[ "$OSTYPE" == "darwin"* ]]; then
        for app in /Applications/ChimeraX*.app; do
            if [[ -d "$app" ]]; then
                chimerax_app="$app"
                # Find Python in the bundle
                for py in "$app"/Contents/bin/python* "$app"/Contents/Library/Frameworks/Python.framework/Versions/*/bin/python*; do
                    if [[ -x "$py" && ! -L "$py" ]]; then
                        chimerax_python="$py"
                        break
                    fi
                done
                break
            fi
        done
        chimerax_bin="$chimerax_app/Contents/MacOS/ChimeraX"

    # Linux locations
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        for path in /usr/lib/chimerax /opt/chimerax /usr/local/chimerax ~/.local/chimerax; do
            if [[ -d "$path" ]]; then
                chimerax_python="$path/bin/python3"
                chimerax_bin="$path/bin/ChimeraX"
                chimerax_app="$path"
                break
            fi
        done

    # Windows (Git Bash / WSL)
    else
        for path in "/c/Program Files/ChimeraX"* "/mnt/c/Program Files/ChimeraX"*; do
            if [[ -d "$path" ]]; then
                chimerax_python="$path/bin/python.exe"
                chimerax_bin="$path/bin/ChimeraX.exe"
                chimerax_app="$path"
                break
            fi
        done
    fi

    echo "$chimerax_python|$chimerax_bin|$chimerax_app"
}

# Main installation
main() {
    echo ""
    echo "=================================="
    echo "  ISIS ChimeraX Plugin Installer"
    echo "=================================="
    echo ""

    # Find ChimeraX
    echo_info "Looking for ChimeraX..."
    IFS='|' read -r CHIMERAX_PYTHON CHIMERAX_BIN CHIMERAX_APP <<< "$(find_chimerax)"

    if [[ -z "$CHIMERAX_PYTHON" || ! -x "$CHIMERAX_PYTHON" ]]; then
        echo_error "Could not find ChimeraX Python interpreter."
        echo ""
        echo "Please specify the path to ChimeraX's Python:"
        echo "  macOS:   /Applications/ChimeraX-1.10.app/Contents/bin/python3.11"
        echo "  Linux:   /usr/lib/chimerax/bin/python3"
        echo "  Windows: C:\\Program Files\\ChimeraX\\bin\\python.exe"
        echo ""
        read -p "ChimeraX Python path: " CHIMERAX_PYTHON

        if [[ ! -x "$CHIMERAX_PYTHON" ]]; then
            echo_error "Invalid path: $CHIMERAX_PYTHON"
            exit 1
        fi
    fi

    echo_info "Found ChimeraX Python: $CHIMERAX_PYTHON"
    echo_info "Found ChimeraX App: $CHIMERAX_APP"
    echo ""

    # Step 1: Install core ISIS library
    echo_info "Installing ISIS core library into ChimeraX Python..."
    "$CHIMERAX_PYTHON" -m pip install --upgrade "$SCRIPT_DIR" 2>&1 | while read line; do
        echo "    $line"
    done

    if [[ $? -ne 0 ]]; then
        echo_error "Failed to install ISIS core library"
        exit 1
    fi
    echo_info "Core library installed successfully"
    echo ""

    # Step 2: Install ChimeraX bundle
    echo_info "Installing ChimeraX plugin bundle..."

    BUNDLE_DIR="$SCRIPT_DIR/src/isis_chimerax"

    if [[ ! -d "$BUNDLE_DIR" ]]; then
        echo_error "Bundle directory not found: $BUNDLE_DIR"
        exit 1
    fi

    # Use ChimeraX to install the bundle
    if [[ -x "$CHIMERAX_BIN" ]]; then
        echo_info "Running: ChimeraX --nogui --cmd 'devel install $BUNDLE_DIR; exit'"
        "$CHIMERAX_BIN" --nogui --exit --cmd "devel install \"$BUNDLE_DIR\"" 2>&1 | while read line; do
            echo "    $line"
        done
    else
        echo_warn "Could not find ChimeraX executable for automatic bundle install."
        echo_warn "Please run this command manually in ChimeraX:"
        echo ""
        echo "    devel install $BUNDLE_DIR"
        echo ""
    fi

    echo ""
    echo "=================================="
    echo_info "Installation complete!"
    echo "=================================="
    echo ""
    echo "Next steps:"
    echo "  1. Restart ChimeraX"
    echo "  2. Open a structure: open 1ubq"
    echo "  3. Run prediction:   isis predict #1"
    echo "  4. Color by scores:  isis color #1"
    echo ""
    echo "For help: isis list"
    echo ""
}

# Uninstall option
uninstall() {
    echo_info "Uninstalling ISIS..."

    IFS='|' read -r CHIMERAX_PYTHON CHIMERAX_BIN CHIMERAX_APP <<< "$(find_chimerax)"

    if [[ -x "$CHIMERAX_PYTHON" ]]; then
        "$CHIMERAX_PYTHON" -m pip uninstall -y isis-epitope 2>/dev/null || true
    fi

    if [[ -x "$CHIMERAX_BIN" ]]; then
        "$CHIMERAX_BIN" --nogui --exit --cmd "toolshed uninstall ChimeraX-ISIS" 2>/dev/null || true
    fi

    echo_info "Uninstall complete. Restart ChimeraX."
}

# Parse arguments
case "${1:-}" in
    --uninstall|-u)
        uninstall
        ;;
    --help|-h)
        echo "Usage: $0 [OPTIONS]"
        echo ""
        echo "Options:"
        echo "  --uninstall, -u    Uninstall ISIS from ChimeraX"
        echo "  --help, -h         Show this help"
        echo ""
        ;;
    *)
        main
        ;;
esac
