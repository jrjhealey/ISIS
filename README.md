# ISIS - In Silico Immunogenicity Studies

Immunogenicity prediction across B-cell, T-cell and innate recognition, with
3D visualisation in ChimeraX. Every prediction class has been benchmarked
against experimental data, and the measured accuracy is reported below --
including where methods do not work.

## Features

- **B-cell linear** -- 5 amino-acid property scales, with consensus voting
- **B-cell conformational** -- DiscoTope / ElliPro / SEPPA styles, from 3D structure
- **T-cell** -- MHC class I and II binding from models trained on IEDB, plus
  proteasomal cleavage and TAP transport
- **Innate** -- N-/O-glycosylation sequons, signal peptides, TLR ligand motifs
- **Structural features** -- SASA, protrusion index, contact number, B-factors
- **Publication-style plots** -- one consistent visual language, from CLI or ChimeraX
- **Benchmarks included** -- reproducible, with the datasets committed
- **Python API**, **command line**, and **ChimeraX plugin**

## Installation

ISIS is two pieces, and this catches people out:

| Piece | What it is | Where it must live |
|---|---|---|
| `isis-epitope` | the prediction library (all the science) | **ChimeraX's own Python**, for the plugin |
| `ChimeraX-ISIS` | the plugin (commands, colouring, plots) | ChimeraX |

Installing only the plugin leaves it running against whatever `isis-epitope` is
already present. If that copy is old, the linear scales still work while T-cell,
conformational and innate prediction all report themselves `UNAVAILABLE` -- which
looks like missing features rather than a broken install. The plugin now declares
`isis-epitope >= 2.1` as a dependency and checks the version at start-up, but if
you hit it, **`isis doctor`** prints the state and the exact repair command.

### For the ChimeraX plugin

```bash
git clone https://github.com/jrjhealey/ISIS.git
cd ISIS
./install_chimerax.sh
```

The script installs the library into ChimeraX's Python **with the `[ml,plot]`
extras**, force-reinstalls the bundle, then verifies by running `isis doctor` and
printing what ChimeraX actually sees. It exits non-zero if verification fails, so
a failed install cannot look like a successful one.

If ChimeraX is somewhere unusual, or you have more than one installed:

```bash
CHIMERAX_PYTHON_OVERRIDE=/Applications/ChimeraX.app/Contents/bin/python3.11 ./install_chimerax.sh
CHIMERAX_PYTHON_OVERRIDE=/usr/lib/ucsf-chimerax/bin/python3.11 ./install_chimerax.sh
```

**Restart ChimeraX afterwards**, then confirm:

```
isis doctor
```

Every component should read `OK` and the core library `2.1.0` or newer.

#### Why not just `devel install`?

`devel install` (and `toolshed install ... reinstall true`) **refuse to replace a
bundle already installed at the same version** -- they stop at "already installed
with the same version" and leave the previous code running while reporting
success. If you are iterating on the plugin, or reinstalling after a change that
did not bump the version, use the script, or force it yourself:

```bash
/path/to/chimerax/bin/python3.11 -m pip install --force-reinstall --no-deps \
    src/isis_chimerax/dist/chimerax_isis-*.whl
```

### For CLI / Python use only

```bash
git clone https://github.com/jrjhealey/ISIS.git
cd ISIS
pip install -e ".[ml,plot]"
```

Extras: `ml` pulls scikit-learn and joblib (needed for the MHC models),
`plot` pulls matplotlib (needed for `isis plot`). ChimeraX 1.12 already ships
matplotlib, so `plot` is only required for CLI plotting.

### Upgrading

Upgrade **both** pieces, or the version check will tell you off:

```
pip install isis-epitope upgrade true     # in ChimeraX
devel install /path/to/ISIS/src/isis_chimerax
```

### Troubleshooting installation

| Symptom | Cause | Fix |
|---|---|---|
| `isis list` shows only the 5 linear scales | old `isis-epitope` without the methods module | `pip install isis-epitope upgrade true` |
| Everything `UNAVAILABLE` | `isis-epitope` not in ChimeraX's Python | same as above |
| Toolshed shows an old version | `devel install` refuses same-version replacement | re-run `./install_chimerax.sh` |
| Installer said OK but nothing changed | same-version bundle not replaced (fixed in 2.1.0) | re-run `./install_chimerax.sh` |
| `isis plot` says matplotlib missing | no matplotlib in ChimeraX's Python | `pip install matplotlib` |

Run `isis doctor` first in every case -- it names the problem and the command.

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

### Measured accuracy

Benchmarks are in `benchmark/`, with datasets committed so results reproduce
without re-downloading. Read this table before trusting any single method.

| Class | Best method | Benchmark | Result |
|---|---|---|---|
| **T-cell MHC-I** | per-allele RF models | held-out IEDB peptides | **AUC 0.82-0.95** (mean 0.87) |
| **T-cell MHC-II** | HLA-DRB1*01:01 | held-out IEDB peptides | AUC 0.76 |
| **B-cell conformational** | `ellipro` | 47 antibody-antigen complexes | AUC 0.68 |
| **B-cell linear** | `kolaskar-tongaonkar` | 49 antigens, IEDB positions | MCC +0.15 |
| **N-glycosylation** | motif scan | spike glycan sites | recovers all 5 known sites |
| Signal peptide | -- | 6 pos/neg controls | 4/6 -- weak |
| O-glycosylation | -- | control proteins | flags 58-83% of S/T -- not usable |
| proteasome, TAP, TLR | -- | -- | not yet benchmarked |

Three findings worth knowing before you pick a method:

1. **MHC-I is the most trustworthy component.** Note that `HLA-A*02:01`, the
   most requested allele, is among the *weakest* installed (AUC 0.82) despite
   having the most training data.
2. **No conformational method beats plain solvent accessibility.** `sasa_only`
   scores AUC 0.669 against ElliPro's 0.684 (p=0.057, not significant), and
   nothing built from these features exceeds ~0.69 under leave-one-antigen-out
   CV. Published SEMA reaches 0.76 using inverse-folding embeddings, so closing
   that gap needs different features rather than reweighting these.
3. **"Best B-cell method" depends on which epitope you mean.**
   Kolaskar-Tongaonkar is the only linear scale beating chance on linear IEDB
   epitopes, yet it is *below* chance on structural antibody-contact patches
   (AUC 0.418) because it favours buried hydrophobic residues. Four of the five
   linear scales sit at or below chance (MCC -0.10 to -0.17) on linear epitopes.

### B-cell linear (sequence only)

| Method | Property | Window | Threshold |
|--------|----------|--------|-----------|
| `emini` | Surface accessibility | 6 | 1.0 |
| `parker` | Hydrophilicity | 7 | average |
| `chou-fasman` | Beta-turn propensity | 7 | average |
| `kolaskar-tongaonkar` | Antigenicity | 7 | 1.0 |
| `karplus-schulz` | Flexibility | 7 | 1.0 |

### B-cell conformational (needs 3D structure)

| Method | Basis |
|--------|-------|
| `discotope` | propensity + contact number + log SASA |
| `ellipro` | protrusion index from the centroid |
| `seppa` | scoring over surface patches |

These need SASA and coordinates, so they run inside ChimeraX (which computes
both from the open structure) rather than from a bare sequence.

### T-cell

`mhc1` and `mhc2` use Random Forest models trained on IEDB binding data.
Run `isis list` (ChimeraX) or `isis list-methods` (CLI) for the alleles actually
installed -- both read the model files rather than a hard-coded list.

| Method | Predicts |
|--------|----------|
| `mhc1` | MHC class I binding, 8-11mers, per allele |
| `mhc2` | MHC class II binding, 15mers with 9mer core |
| `proteasome` | proteasomal cleavage sites |
| `tap` | TAP transport efficiency |
| `consensus` | weighted MHC + cleavage + TAP pipeline |

### Innate

| Method | Predicts |
|--------|----------|
| `glyco` | N-glycosylation sequons (N-X-S/T) and O-glycosylation sites |
| `signal` | signal peptide and cleavage position |
| `tlr` | TLR ligand motifs, LPS-binding and basic/hydrophobic patches |

### Chain and residue specifications

Every ChimeraX command takes a residue spec, so a chain selector is honoured:

```
isis bcell linear #1        # every chain
isis bcell linear #1/A      # chain A only
isis bcell linear #1/A,C    # chains A and C
```

Colouring, export and clearing are confined to the chains named. A spec that
clips a chain part-way (`#1/A:100-200`) selects that chain but does **not**
truncate it -- scoring always runs over the full chain sequence, because a
sliding window over a fragment would not reproduce the scores it has in the
whole protein.
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

### `isis bcell conformational` - Structure-Based Epitopes

Requires 3D structure; SASA and coordinates are computed from the open model.

```
isis bcell conformational #1/A method discotope
isis bcell conformational #1/A method ellipro
isis bcell conformational #1/A method seppa
```

`discotope` derives its cut-off from the protein's own score distribution
(top 15% of residues). It does **not** use the published DiscoTope 2.0 value of
-7.7: that belongs to a different score scale, and against this implementation's
range it labelled every residue an epitope.

### `isis tcell` - T-cell Epitopes

```
isis tcell mhc1 #1/A allele HLA-A*02:01
isis tcell mhc2 #1/A allele HLA-DRB1*01:01
isis tcell proteasome #1/A
isis tcell tap #1/A
isis tcell consensus #1/A allele HLA-A*02:01
```

Run `isis list` for installed alleles. Requesting one that is not installed
raises an error naming the available options.

### `isis innate` - Innate Recognition

```
isis innate glyco #1/A
isis innate glyco #1/A glycoType o
isis innate signal #1/A
isis innate tlr #1/A
isis innate consensus #1/A
```

### `isis structure` - Structural Features

```
isis structure sasa #1/A
isis structure protrusion #1/A
isis structure contacts #1/A cutoff 8.0
isis structure bfactor #1/A
```

### `isis bcell consensus` - Combine All B-cell Methods

```
isis bcell consensus #1/A
isis bcell consensus #1/A minMethods 3 minLength 8
```

Votes across the linear scales and, when a structure is present, the
conformational methods too, storing the count per residue. The benchmark shows
consensus voting does not rescue the weak linear scales, so treat a high vote
count as agreement rather than as evidence of accuracy.

### `isis export` - Scores to File

```
isis export #1/A format csv
isis export #1/A format json output scores.json
```

Exports every ISIS score stored on the selected chains, one column per method,
aligned to sequence position. Positions with no modelled residue export as 0.

### `isis plot` - Figures From the Structure

No FASTA needed -- the sequence comes from the structure's own chain sequence.

```
isis plot #1                       # all chains, all methods
isis plot #1/A                     # chain A only
isis plot #1/A method emini,parker window 9 outdir figures/ prefix myprotein
```

Writes a per-residue profile, a method call-matrix and a consensus track per
chain, logging clickable links. Needs matplotlib in ChimeraX's Python (it ships
with ChimeraX 1.12).

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

---

## Benchmarks

Reproducible, with datasets committed so nothing needs re-downloading.

### Linear B-cell epitopes -- 49 antigens

```bash
python3 benchmark/run_benchmark.py
```

49 antigens across 32 organisms (24,183 residues, 10,657 IEDB-confirmed epitope
positions). Sequences fetched from NCBI, positions verified against the actual
sequence.

| Method | Sens | Spec | F1 | MCC |
|---|---|---|---|---|
| kolaskar-tongaonkar | 0.63 | 0.53 | 0.52 | **+0.148** |
| karplus-schulz | 0.35 | 0.51 | 0.32 | -0.135 |
| chou-fasman | 0.28 | 0.56 | 0.30 | -0.133 |
| parker | 0.28 | 0.53 | 0.29 | -0.169 |
| emini | 0.16 | 0.75 | 0.20 | -0.104 |

MCC is 0 for random guessing. Consensus voting does not rescue the weak scales
at any threshold from 1 to 5.

To rebuild the dataset from scratch (needs the IEDB dump and network):

```bash
python3 benchmark/scan_iedb.py
python3 benchmark/fetch_sequences.py
python3 benchmark/select_final_set.py
```

### Conformational epitopes -- 47 antibody-antigen complexes

```bash
python3 benchmark/conformational/build_dataset.py     # fetches from RCSB
python3 benchmark/conformational/run_benchmark.py
python3 benchmark/conformational/fit_consensus.py
```

Complexes from Table 2 of Ponomarenko & Bourne 2007, the canonical benchmark for
this task. Whole-protein antigens (61-613 aa: lysozyme, hemagglutinin,
neuraminidase, EGFR, VEGF-A, prion protein, OspA, CD3 and others), 8,706
residues, 832 epitope residues.

Ground truth is an observed physical contact -- an antigen residue with any atom
within 4 Å of any antibody atom -- not a curated annotation.

Structural features are computed on the antigen **alone**, with the antibody
deleted. Computing SASA on the complex would bury exactly the epitope residues,
making the feature encode the label.

| Method | AUC | vs plain SASA |
|---|---|---|
| ellipro | 0.684 | +0.015 (p=0.057, ns) |
| *sasa_only (control)* | *0.669* | -- |
| discotope | 0.659 | -0.010 (ns) |
| seppa | 0.612 | -0.056 (worse) |
| best fitted combination | 0.689 | +0.021 (p=0.117, ns) |

Reference: SEMA 0.76, random 0.50.

### Rendering

```bash
python3 benchmark/render_ghost.py            # inside ChimeraX
python3 benchmark/make_contact_sheets.py
```

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
