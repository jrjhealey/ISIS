# ISIS - In Silico Immunogenicity Studies

B-cell epitope prediction using established amino acid property scales. Predict which regions of a protein are likely to be recognized by antibodies.

## Features

- **5 prediction methods** based on peer-reviewed amino acid property scales
- **Consensus analysis** combining all methods for higher confidence
- **Python API** for scripting and pipelines
- **Command-line interface** for batch processing
- **ChimeraX plugin** for 3D visualization on protein structures

---

## Installation

### Core Library

```bash
git clone https://github.com/jrjhealey/ISIS.git
cd ISIS
pip install -e .
```

### ChimeraX Plugin

#### Easy Install (Recommended)

```bash
# macOS/Linux - run from the ISIS directory:
./install_chimerax.sh

# Or use Python directly with ChimeraX's interpreter:
/Applications/ChimeraX-1.10.app/Contents/bin/python3.11 install_chimerax.py
```

The installer will:
1. Find your ChimeraX installation
2. Install the ISIS core library into ChimeraX's Python
3. Install the ChimeraX plugin bundle

**Restart ChimeraX after installation.**

#### Manual Install

If the automatic installer doesn't work:

1. **Install the core library into ChimeraX's Python:**

```bash
# macOS
/Applications/ChimeraX-1.10.app/Contents/bin/python3.11 -m pip install /path/to/ISIS

# Linux
/usr/lib/chimerax/bin/python3 -m pip install /path/to/ISIS

# Windows
"C:\Program Files\ChimeraX\bin\python.exe" -m pip install C:\path\to\ISIS
```

2. **Install the ChimeraX bundle** (in ChimeraX command line):

```
devel install /path/to/ISIS/src/isis_chimerax
```

3. **Restart ChimeraX**

#### Uninstall

```bash
./install_chimerax.sh --uninstall
```

Or manually in ChimeraX:
```
toolshed uninstall ChimeraX-ISIS
```

---

## Quick Start

### Python

```python
from isis import predict, predict_all

result = predict("MKTAYIAKQRQISFVKSHFSRQLE", method="emini")
for ep in result.epitopes:
    print(f"{ep.start}-{ep.end}: {ep.sequence}")
```

### Command Line

```bash
isis predict sequence.fasta
isis list-methods
```

### ChimeraX

```
open 1ubq
isis consensus #1
```

---

## Prediction Methods

| Method | Property | Window | Threshold | Use Case |
|--------|----------|--------|-----------|----------|
| `emini` | Surface accessibility | 6 | 1.0 | Surface-exposed regions |
| `parker` | Hydrophilicity | 7 | average | Hydrophilic/antigenic regions |
| `chou-fasman` | Beta-turn propensity | 7 | average | Loop regions |
| `kolaskar-tongaonkar` | Antigenicity | 7 | 1.0 | ~75% accuracy on known epitopes |
| `karplus-schulz` | Flexibility | 7 | 1.0 | Flexible, accessible regions |

---

## ChimeraX Commands

### `isis predict` - Single Method Prediction

Run epitope prediction using one method. Stores scores as residue attributes.

**Syntax:**
```
isis predict <structures> [method <name>] [window <int>] [threshold <float>]
```

**Examples:**
```
# Open a structure first
open 1ubq

# Default prediction (Emini method)
isis predict #1

# Specify method
isis predict #1 method parker

# Custom window size
isis predict #1 method emini window 9

# Custom threshold (stricter)
isis predict #1 method kolaskar-tongaonkar threshold 1.2

# Multiple structures
isis predict #1-3 method emini
```

---

### `isis consensus` - Multi-Method Consensus

Run ALL methods and create a consensus score (0-5) showing how many methods agree each position is epitopic. Automatically colors the structure.

**Syntax:**
```
isis consensus <structures> [min_methods <int>] [min_length <int>]
```

**Parameters:**
- `min_methods`: Minimum methods that must agree (default: 2)
- `min_length`: Minimum epitope length in amino acids (default: 6)

**Examples:**
```
open 1ubq

# Basic consensus (all defaults)
isis consensus #1

# Require 3+ methods to agree
isis consensus #1 min_methods 3

# Require longer epitopes (8+ aa)
isis consensus #1 min_length 8

# Strict consensus
isis consensus #1 min_methods 4 min_length 8
```

**Output:**
- Colors structure: white(0) → yellow(1-2) → orange(3-4) → red(5)
- Logs consensus epitopes with sequences and average method agreement
- If no consensus found with default thresholds, automatically loosens thresholds

---

### `isis color` - Color by Scores

Color structure using a gradient based on prediction scores.

**Syntax:**
```
isis color <structures> [method <name>] [palette <colors>]
```

**Palette format:** `low:mid:high` or `low:high`

**Examples:**
```
# First run a prediction
isis predict #1 method emini

# Default coloring (white → yellow → red)
isis color #1

# Specify which method's scores to use
isis color #1 method emini

# Custom color palette
isis color #1 palette blue:white:red
isis color #1 palette white:red
isis color #1 palette cyan:magenta
isis color #1 method parker palette green:yellow:red
```

---

### `isis epitopes` - Highlight Epitope Regions

Color only the predicted epitope regions (above threshold), not the full gradient.

**Syntax:**
```
isis epitopes <structures> [method <name>] [color <color>]
```

**Examples:**
```
# First run a prediction
isis predict #1 method emini

# Highlight epitopes in red
isis epitopes #1 color red

# Different method and color
isis epitopes #1 method parker color orange

# Hex color
isis epitopes #1 color #ff6600
```

---

### `isis list` - List Available Methods

Show all available prediction methods.

**Example:**
```
isis list
```

**Output:**
```
Available ISIS prediction methods:
  emini: Emini Surface Accessibility
  parker: Parker Hydrophilicity
  chou-fasman: Chou-Fasman Beta-Turn
  kolaskar-tongaonkar: Kolaskar-Tongaonkar Antigenicity
  karplus-schulz: Karplus-Schulz Flexibility
```

---

### `isis clear` - Remove Predictions

Remove all ISIS prediction attributes from structures.

**Syntax:**
```
isis clear <structures>
```

**Examples:**
```
# Clear all ISIS attributes
isis clear #1

# Clear multiple structures
isis clear #1-3

# Start fresh
isis clear #1
color #1 white
```

---

## Command Line Reference

### `isis predict`

```bash
isis predict <input> [options]

Arguments:
  input                 FASTA file or - for stdin

Options:
  -m, --methods TEXT    Comma-separated methods (default: all)
  -w, --window INT      Window size (default: method-specific)
  -f, --format TEXT     Output format: table, csv, json, epitopes
  -o, --output FILE     Output file (default: stdout)
```

**Examples:**

```bash
# Basic prediction with table output
isis predict protein.fasta

# Specific methods
isis predict protein.fasta -m emini,parker

# CSV output for spreadsheet
isis predict protein.fasta -f csv -o results.csv

# JSON for programming
isis predict protein.fasta -f json -o results.json

# Compact epitope-only output
isis predict protein.fasta -f epitopes

# From stdin
cat protein.fasta | isis predict -

# Custom window size
isis predict protein.fasta -w 9
```

### `isis list-methods`

```bash
isis list-methods
```

---

## Python API Reference

### `predict()`

```python
from isis import predict

result = predict(
    sequence,           # Amino acid sequence (str)
    method="emini",     # Prediction method
    window_size=None,   # Window size (default: method-specific)
    threshold=None,     # Score threshold (default: method-specific)
    sequence_name="Seq" # Identifier
)

# Access results
result.scores          # numpy array of per-position scores
result.positions       # numpy array of positions (1-indexed)
result.threshold       # threshold used
result.epitopes        # list of Epitope objects
result.score_at(10)    # score at position 10
result.to_dict()       # serialize to dictionary
```

### `predict_all()`

```python
from isis import predict_all

results = predict_all(
    sequence,
    methods=["emini", "parker"],  # or None for all
    window_size=None
)

for method, result in results.items():
    print(f"{method}: {len(result.epitopes)} epitopes")
```

### `Epitope` Object

```python
epitope.start      # Start position (1-indexed)
epitope.end        # End position (inclusive)
epitope.sequence   # Amino acid sequence
epitope.score      # Average score
epitope.length     # Length
```

---

## Worked Examples

### Example 1: Quick Consensus Analysis

```
open 1ubq
isis consensus #1
```

Output:
```
Running consensus analysis on 1ubq...
  Chain A: Found consensus with threshold multiplier 0.9
  Chain A: 2 consensus epitope(s)
    12-18: TITLEVP (avg 3.2 methods)
    44-52: LEDGRTLSD (avg 2.8 methods)
```

### Example 2: Compare Individual Methods

```
open 3sgb

# Run each method
isis predict #1 method emini
isis predict #1 method kolaskar-tongaonkar

# View Emini scores
isis color #1 method emini palette white:blue

# Switch to Kolaskar
isis color #1 method kolaskar-tongaonkar palette white:red
```

### Example 3: Surface Visualization

```
open 4hhb

# Predict and color
isis consensus #1 min_methods 3

# Add transparent surface
surface #1
transparency #1 60

# Save image
view all
save ~/epitopes.png supersample 3
```

### Example 4: Batch Processing (CLI)

```bash
# Process all FASTA files
for f in *.fasta; do
    isis predict "$f" -f json -o "${f%.fasta}.json"
done

# Or with parallel
ls *.fasta | parallel 'isis predict {} -f csv -o {.}.csv'
```

### Example 5: Python Consensus Analysis

```python
from isis import predict_all, available_methods
from collections import defaultdict

sequence = "MKTAYIAKQRQISFVKSHFSRQLEEALCLSLHR"
results = predict_all(sequence)

# Count method agreement per position
votes = defaultdict(int)
for method, result in results.items():
    for ep in result.epitopes:
        for pos in range(ep.start, ep.end + 1):
            votes[pos] += 1

# Find consensus positions (3+ methods)
consensus = [pos for pos, count in votes.items() if count >= 3]
print(f"Consensus positions: {consensus}")
```

---

## Troubleshooting

### ChimeraX: "Unknown command: isis"

1. Ensure bundle is installed: `toolshed list installed`
2. Restart ChimeraX after installation
3. Check for errors: `log show`

### ChimeraX: "ISIS library not installed"

Install ISIS into ChimeraX's Python:
```bash
/Applications/ChimeraX-1.10.app/Contents/bin/python3.11 -m pip install /path/to/ISIS
```

### ChimeraX: Commands still not working

Manually register in Python shell (Tools → General → Shell):
```python
from chimerax.isis.cmd import register_all_commands
register_all_commands(session)
```

### "Sequence too short"

Sequences must be at least 6 amino acids (the minimum epitope length).

### No epitopes predicted

- The sequence may lack strong epitope signals
- Try `isis consensus` which auto-loosens thresholds
- Or manually lower threshold: `isis predict #1 threshold 0.8`

---

## Legacy Code

The original Python 2 / UCSF Chimera code is preserved on the `legacy/v1-python2` branch:

```bash
git checkout legacy/v1-python2
```

---

## Citation

If you use ISIS in your research, please cite the original method papers:

- **Emini:** Emini EA et al. J Virol. 1985;55(3):836-839
- **Parker:** Parker JM et al. Biochemistry. 1986;25(19):5425-5432
- **Chou-Fasman:** Chou PY, Fasman GD. Adv Enzymol. 1978;47:45-148
- **Kolaskar-Tongaonkar:** Kolaskar AS, Tongaonkar PC. FEBS Lett. 1990;276(1-2):172-174
- **Karplus-Schulz:** Karplus PA, Schulz GE. Naturwissenschaften. 1985;72:212-213

---

## License

MIT License

## Contributing

Issues and pull requests welcome at https://github.com/jrjhealey/ISIS
