# ISIS - In Silico Immunogenicity Studies

B-cell epitope prediction using established amino acid property scales. Predict which regions of a protein are likely to be recognized by antibodies.

## Features

- **5 prediction methods** based on peer-reviewed amino acid scales
- **Python API** for scripting and pipelines
- **Command-line interface** for batch processing
- **ChimeraX plugin** for 3D visualization on protein structures

---

## Installation

### Core Library

```bash
# From PyPI (when published)
pip install isis-epitope

# From source
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
4. Restart ChimeraX to complete

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

# Single method
result = predict("MKTAYIAKQRQISFVKSHFSRQLE", method="emini")
print(f"Found {len(result.epitopes)} epitopes")
for ep in result.epitopes:
    print(f"  {ep.start}-{ep.end}: {ep.sequence} (score={ep.score:.2f})")

# All methods
results = predict_all("MKTAYIAKQRQISFVKSHFSRQLE")
for method, result in results.items():
    print(f"{method}: {len(result.epitopes)} epitopes")
```

### Command Line

```bash
# Basic prediction
isis predict sequence.fasta

# Multiple methods, CSV output
isis predict sequence.fasta --methods emini,parker,kolaskar-tongaonkar --format csv

# JSON output for downstream processing  
isis predict sequence.fasta --format json --output results.json

# List available methods
isis list-methods
```

### ChimeraX

```
open 1ubq
isis predict #1
isis color #1
```

---

## Prediction Methods

| Method | Property | Default Window | Threshold | Reference |
|--------|----------|----------------|-----------|-----------|
| `emini` | Surface accessibility | 6 | 1.0 | Emini et al. 1985 |
| `parker` | Hydrophilicity | 7 | average | Parker et al. 1986 |
| `chou-fasman` | Beta-turn propensity | 7 | average | Chou & Fasman 1978 |
| `kolaskar-tongaonkar` | Antigenicity | 7 | 1.0 | Kolaskar & Tongaonkar 1990 |
| `karplus-schulz` | Flexibility | 7 | 1.0 | Karplus & Schulz 1985 |

### Method Selection Guide

- **`emini`** - Best for identifying surface-exposed regions. Good starting point.
- **`parker`** - Hydrophilic regions tend to be antigenic. Complements Emini.
- **`kolaskar-tongaonkar`** - Semi-empirical, ~75% accuracy on known epitopes.
- **`chou-fasman`** - Beta-turns are often in loop regions accessible to antibodies.
- **`karplus-schulz`** - Flexible regions can adapt to antibody binding.

**Recommendation:** Run multiple methods and look for consensus regions.

---

## Python API Reference

### `predict(sequence, method, window_size, threshold, sequence_name)`

Predict epitopes for a single sequence.

**Parameters:**
- `sequence` (str): Amino acid sequence
- `method` (str): Prediction method (default: "emini")
- `window_size` (int): Sliding window size (default: method-specific)
- `threshold` (float): Score threshold for epitope calls (default: method-specific)
- `sequence_name` (str): Identifier for the sequence

**Returns:** `Prediction` object

```python
from isis import predict

result = predict(
    "MKTAYIAKQRQISFVKSHFSRQLE",
    method="emini",
    window_size=7,
    threshold=1.2
)

# Access results
print(result.scores)      # numpy array of scores
print(result.positions)   # numpy array of positions (1-indexed)
print(result.threshold)   # threshold used
print(result.epitopes)    # list of Epitope objects
```

### `predict_all(sequence, methods, window_size, sequence_name)`

Run multiple prediction methods.

```python
from isis import predict_all

results = predict_all(
    "MKTAYIAKQRQISFVKSHFSRQLE",
    methods=["emini", "parker", "kolaskar-tongaonkar"]
)

for method, result in results.items():
    print(f"{method}: {len(result.epitopes)} epitopes")
```

### `Prediction` Object

```python
result.method          # Method name
result.sequence        # Input sequence
result.sequence_name   # Sequence identifier
result.window_size     # Window size used
result.threshold       # Threshold used
result.positions       # numpy array, 1-indexed center positions
result.scores          # numpy array, per-position scores
result.epitopes        # List of Epitope objects above threshold

result.score_at(10)    # Get score at position 10 (or None)
result.to_dict()       # Serialize to dictionary
```

### `Epitope` Object

```python
epitope.start      # Start position (1-indexed)
epitope.end        # End position (1-indexed, inclusive)
epitope.sequence   # Amino acid sequence
epitope.score      # Average score
epitope.length     # Length of epitope
```

---

## Command Line Reference

### `isis predict`

```bash
isis predict <input> [options]

Arguments:
  input                 FASTA file or - for stdin

Options:
  -m, --methods TEXT    Comma-separated methods (default: emini,parker,chou-fasman,kolaskar-tongaonkar)
  -w, --window INT      Window size (default: method-specific)
  -f, --format TEXT     Output format: table, csv, json, epitopes (default: table)
  -o, --output FILE     Output file (default: stdout)
```

**Examples:**

```bash
# Default table output
isis predict protein.fasta

# Specific methods
isis predict protein.fasta -m emini,parker

# CSV for spreadsheet
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

Shows all available methods with descriptions.

---

## ChimeraX Plugin Reference

### Commands

#### `isis predict <structures> [method <name>] [window <int>] [threshold <float>]`

Run epitope prediction on structure sequences. Scores are stored as residue attributes.

```
isis predict #1
isis predict #1 method emini
isis predict #1 method parker window 9
isis predict #1 method kolaskar-tongaonkar threshold 1.1
```

#### `isis color <structures> [method <name>] [palette <colors>]`

Color structure by prediction scores using a gradient.

```
isis color #1
isis color #1 method emini
isis color #1 palette white:yellow:red
isis color #1 palette blue:white:red
isis color #1 method parker palette cyan:magenta
```

**Palette format:** `low:high` or `low:mid:high`

#### `isis epitopes <structures> [method <name>] [color <color>]`

Highlight only the predicted epitope regions.

```
isis epitopes #1 color red
isis epitopes #1 method parker color orange
isis epitopes #1 color #ff6600
```

#### `isis list`

Show available prediction methods.

```
isis list
```

#### `isis consensus <structures> [min_methods <int>] [threshold <float>]`

Run ALL methods and create a consensus score (0-5) based on how many methods agree each position is epitopic. Automatically colors the structure.

```
isis consensus #1
isis consensus #1 min_methods 3
isis consensus #1 threshold 1.2
```

**Scoring:**
- 0 = no methods predict epitope
- 1-2 = weak consensus
- 3-4 = moderate consensus
- 5 = all methods agree (highest confidence)

**Colors:** white (0) → yellow (1-2) → orange (3-4) → red (5)

#### `isis clear <structures>`

Remove all ISIS prediction attributes from structures.

```
isis clear #1
```

### Structure Specifiers

```
#1              First open model
#2              Second open model
#1-3            Models 1, 2, and 3
#1/A            Chain A of model 1
#1/A,B          Chains A and B of model 1
```

---

## Worked Examples

### Example 1: Analyze a Single Protein

**Goal:** Find potential B-cell epitopes in ubiquitin.

**Python:**
```python
from isis import predict_all

# Ubiquitin sequence
ubiquitin = """
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ
KESTLHLVLRLRGG
"""

results = predict_all(ubiquitin.replace("\n", ""))

print("Predicted epitopes by method:\n")
for method, result in results.items():
    print(f"{method}:")
    if result.epitopes:
        for ep in result.epitopes:
            print(f"  {ep.start:3d}-{ep.end:3d}: {ep.sequence}")
    else:
        print("  No epitopes above threshold")
    print()
```

**CLI:**
```bash
echo ">ubiquitin
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ
KESTLHLVLRLRGG" > ubiquitin.fasta

isis predict ubiquitin.fasta
```

**ChimeraX:**
```
open 1ubq
isis predict #1
isis color #1 palette white:yellow:red
```

---

### Example 2: Compare Multiple Methods

**Goal:** Find consensus epitope regions across methods.

**Python:**
```python
from isis import predict_all
from collections import defaultdict

sequence = "MKTAYIAKQRQISFVKSHFSRQLEEALCLSLHRALQFGPRVLVVS"

results = predict_all(sequence)

# Count how many methods predict each position as epitopic
position_votes = defaultdict(int)

for method, result in results.items():
    for ep in result.epitopes:
        for pos in range(ep.start, ep.end + 1):
            position_votes[pos] += 1

# Find consensus regions (predicted by 3+ methods)
print("Consensus epitope positions (3+ methods):")
consensus = [pos for pos, votes in position_votes.items() if votes >= 3]
if consensus:
    # Group consecutive positions
    start = consensus[0]
    for i in range(1, len(consensus)):
        if consensus[i] != consensus[i-1] + 1:
            print(f"  {start}-{consensus[i-1]}: {sequence[start-1:consensus[i-1]]}")
            start = consensus[i]
    print(f"  {start}-{consensus[-1]}: {sequence[start-1:consensus[-1]]}")
else:
    print("  No consensus regions found")
```

**ChimeraX:**
```
open 3sgb

# Run multiple methods
isis predict #1 method emini
isis predict #1 method parker  
isis predict #1 method kolaskar-tongaonkar

# View each one
isis color #1 method emini palette white:blue
# screenshot ~/emini.png

isis color #1 method parker palette white:green
# screenshot ~/parker.png

isis color #1 method kolaskar-tongaonkar palette white:red
# screenshot ~/kolaskar.png
```

---

### Example 3: Batch Processing Multiple Sequences

**Python:**
```python
from isis import predict
from Bio import SeqIO  # requires biopython

results = []

for record in SeqIO.parse("proteins.fasta", "fasta"):
    result = predict(str(record.seq), method="emini", sequence_name=record.id)
    results.append({
        "id": record.id,
        "length": len(record.seq),
        "n_epitopes": len(result.epitopes),
        "epitopes": [(ep.start, ep.end, ep.sequence) for ep in result.epitopes]
    })

# Summary
for r in results:
    print(f"{r['id']}: {r['n_epitopes']} epitopes in {r['length']} aa")
```

**CLI:**
```bash
# Process all FASTA files in a directory
for f in *.fasta; do
    echo "=== $f ==="
    isis predict "$f" -f epitopes
done > all_epitopes.txt

# Or with parallel
ls *.fasta | parallel 'isis predict {} -f json -o {.}.json'
```

---

### Example 4: Visualize on Protein Structure

**ChimeraX workflow:**

```
# 1. Open structure
open 4hhb  # Hemoglobin

# 2. Predict with stringent threshold
isis predict #1 method emini threshold 1.3

# 3. Style the structure
preset cartoons/nucleotides
color #1 white

# 4. Color epitopes
isis epitopes #1 color red

# 5. Add surface
surface #1
transparency #1 70

# 6. Save image
view all
windowsize 1920 1080
save ~/hemoglobin_epitopes.png supersample 3
```

---

### Example 5: Multi-chain Analysis

**ChimeraX:**
```
# Open antibody-antigen complex
open 1igc

# Predict on all chains
isis predict #1

# Color each chain differently
isis color #1/A method emini palette white:red
isis color #1/B method emini palette white:blue

# Or just highlight epitopes on antigen chain
isis epitopes #1/A color yellow
```

**Python:**
```python
from isis import predict
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure("complex", "1igc.pdb")

for model in structure:
    for chain in model:
        # Build sequence from residues
        residues = [r for r in chain.get_residues() if r.id[0] == " "]
        seq = "".join(r.resname[0] if len(r.resname) == 3 else "X" for r in residues)
        
        if len(seq) < 10:
            continue
            
        result = predict(seq, method="emini")
        print(f"\nChain {chain.id} ({len(seq)} aa):")
        for ep in result.epitopes:
            print(f"  {ep.start}-{ep.end}: {ep.sequence}")
```

---

### Example 6: Export for External Tools

**JSON export for web visualization:**
```bash
isis predict protein.fasta -f json -o epitopes.json
```

```python
import json

with open("epitopes.json") as f:
    data = json.load(f)

# Access structured data
for entry in data:
    seq_name = entry["sequence_name"]
    for method, pred in entry["predictions"].items():
        scores = pred["scores"]
        epitopes = pred["epitopes"]
        # ... process for visualization
```

**CSV for spreadsheet analysis:**
```bash
isis predict protein.fasta -f csv -o epitopes.csv
```

---

## Output Formats

### Table (default)
Human-readable format with scores and epitope summary.

### CSV
Spreadsheet-compatible with position, residue, and scores per method.

### JSON
Structured data including all scores, positions, and epitope details.

### Epitopes
Compact TSV with only epitope calls: `method start end sequence score`

---

## Interpreting Results

### Scores
- **Emini, Kolaskar-Tongaonkar, Karplus-Schulz:** Values > 1.0 suggest epitope potential
- **Parker, Chou-Fasman:** Above-average values suggest epitope potential

### Epitopes
Contiguous regions of 6+ residues above the threshold are reported as epitopes.

### Confidence
- Single method prediction: Low confidence
- Multiple methods agree: Higher confidence
- Consensus across 3+ methods: Highest confidence

### Limitations
- These are **linear epitope** predictions only
- Conformational epitopes require 3D structure analysis
- Predictions are probabilistic, not definitive
- Experimental validation is always recommended

---

## Troubleshooting

### ChimeraX: "Unknown command: isis"
1. Ensure the bundle is installed: `toolshed list installed`
2. Restart ChimeraX after installation
3. Check for errors: `log show`

### ChimeraX: "ISIS library not installed"
Install ISIS into ChimeraX's Python:
```bash
/Applications/ChimeraX-1.10.app/Contents/bin/python3.11 -m pip install /path/to/ISIS
```

### "Sequence too short"
Sequences must be at least as long as the window size (typically 6-7 amino acids).

### No epitopes predicted
- Try a lower threshold: `isis predict #1 threshold 0.8`
- Try a different method
- The sequence may genuinely lack strong epitope signals

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
