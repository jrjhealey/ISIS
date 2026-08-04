#!/usr/bin/env python3
"""Fetch full sequences for candidate antigens and validate epitope positions."""
import json
import subprocess
import time
import re

with open('benchmark/dataset/candidate_antigens.json') as f:
    candidates = json.load(f)

EXCLUDE_TERMS = ['polyprotein', 'hypothetical', 'unknown', 'unnamed', 'conserved hypothetical']
filtered = []
for acc, data in candidates.items():
    name = data['source_mol'].lower()
    if any(term in name for term in EXCLUDE_TERMS):
        continue
    filtered.append((acc, data))
filtered.sort(key=lambda x: -x[1]['positive_count'])

# Take more than 50 as buffer since some fetches/validations will fail
pool = filtered[:84]

def fetch_fasta(accession):
    """Fetch a FASTA sequence from NCBI via curl."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={accession}&rettype=fasta&retmode=text"
    try:
        result = subprocess.run(['curl', '-s', '--max-time', '20', url],
                                 capture_output=True, text=True, timeout=25)
        lines = result.stdout.strip().split('\n')
        if not lines or not lines[0].startswith('>'):
            return None
        seq = ''.join(lines[1:]).replace(' ', '').upper()
        seq = re.sub(r'[^A-Z]', '', seq)
        return seq
    except Exception:
        return None

validated = []
print(f"Fetching {len(pool)} candidate sequences from NCBI...\n")

for i, (acc, data) in enumerate(pool):
    seq = fetch_fasta(acc)
    time.sleep(0.34)  # NCBI rate limit ~3/sec

    if not seq or len(seq) < 60 or len(seq) > 1800:
        print(f"  [{i+1}/{len(pool)}] {acc}: SKIP (len={len(seq) if seq else 0})")
        continue

    # Validate peptide positions against fetched sequence
    matches = 0
    valid_positions = set()
    for pep, start, end in data['peptides']:
        if start < 1 or end > len(seq) or start > end:
            continue
        candidate_slice = seq[start-1:end]
        if candidate_slice == pep:
            matches += 1
            for p in range(start-1, end):
                valid_positions.add(p)
        elif pep in seq:
            # position off but peptide exists somewhere - use actual position
            pos = seq.find(pep)
            matches += 1
            for p in range(pos, pos+len(pep)):
                valid_positions.add(p)

    match_rate = matches / len(data['peptides']) if data['peptides'] else 0

    if matches >= 4 and len(valid_positions) >= 8:
        validated.append({
            'accession': acc,
            'organism': data['organism'],
            'source_mol': data['source_mol'],
            'sequence': seq,
            'seq_len': len(seq),
            'epitope_positions': sorted(valid_positions),
            'n_peptides_matched': matches,
            'n_peptides_total': len(data['peptides']),
            'match_rate': match_rate,
        })
        print(f"  [{i+1}/{len(pool)}] {acc}: OK len={len(seq)} matched={matches}/{len(data['peptides'])} gt_pos={len(valid_positions)}")
    else:
        print(f"  [{i+1}/{len(pool)}] {acc}: LOW MATCH ({matches}/{len(data['peptides'])})")

print(f"\n\nValidated {len(validated)} antigens with confirmed ground truth positions")

with open('benchmark/dataset/validated_antigens.json', 'w') as f:
    json.dump(validated, f)

print("Saved to benchmark/dataset/validated_antigens.json")
