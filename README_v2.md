# ISIS v2.0 - In Silico Immunogenicity Studies

B-cell epitope prediction using established amino acid property scales.

## Installation

```bash
pip install isis-epitope
```

For ChimeraX integration:
```bash
# From ChimeraX command line:
toolshed install isis
```

## Quick Start

### Python API

```python
from isis import predict, predict_all

# Single method
result = predict("MKTAYIAKQRQISFVKSHFSRQLE", method="emini")
print(f"Threshold: {result.threshold}")
for epitope in result.epitopes:
    print(f"  {epitope.start}-{epitope.end}: {epitope.sequence}")

# All methods
results = predict_all("MKTAYIAKQRQISFVKSHFSRQLE")
for method, result in results.items():
    print(f"{method}: {len(result.epitopes)} epitopes")
```

### Command Line

```bash
# Predict epitopes from FASTA file
isis predict sequence.fasta

# Specific methods and output format
isis predict sequence.fasta --methods emini,parker --format csv --output results.csv

# List available methods
isis list-methods
```

### ChimeraX Commands

```
# Open a structure, then:
isis predict #1                          # Predict with default method
isis predict #1 method emini window 7    # Custom parameters
isis color #1 method emini               # Color by scores
isis epitopes #1 method emini color red  # Highlight epitopes only
isis clear #1                            # Remove predictions
isis list                                # Show available methods
```

## Prediction Methods

| Method | Property | Default Window | Reference |
|--------|----------|----------------|-----------|
| `emini` | Surface accessibility | 6 | Emini et al. 1985 |
| `parker` | Hydrophilicity | 7 | Parker et al. 1986 |
| `chou-fasman` | Beta-turn propensity | 7 | Chou & Fasman 1978 |
| `kolaskar-tongaonkar` | Antigenicity | 7 | Kolaskar & Tongaonkar 1990 |
| `karplus-schulz` | Flexibility | 7 | Karplus & Schulz 1985 |

## API Reference

### `predict(sequence, method="emini", window_size=None, threshold=None)`

Predict epitopes for a single sequence.

**Returns:** `Prediction` object with:
- `.scores` - numpy array of per-position scores
- `.positions` - numpy array of 1-indexed positions
- `.threshold` - threshold used for epitope calling
- `.epitopes` - list of `Epitope` objects above threshold
- `.to_dict()` - serialize to dictionary

### `predict_all(sequence, methods=None, window_size=None)`

Run multiple prediction methods.

**Returns:** dict mapping method name to `Prediction`

### `Epitope`

Dataclass representing a predicted epitope:
- `.start` - 1-indexed start position
- `.end` - 1-indexed end position (inclusive)
- `.sequence` - amino acid sequence
- `.score` - average score
- `.length` - epitope length

## Development

```bash
git clone https://github.com/nanosyrinx/isis
cd isis
pip install -e ".[dev]"
pytest
```

## License

MIT License
