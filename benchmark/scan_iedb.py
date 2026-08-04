#!/usr/bin/env python3
"""
Scan the local IEDB mhc_ligand_full.csv dump for antigens that have an NCBI
protein accession and enough positively-tested peptides to serve as a
benchmark ground-truth source. Excludes generic host/model organisms and
unnamed/polyprotein entries so the final list is diverse, well-characterized
proteins with a fetchable full-length sequence.

Output: benchmark/dataset/candidate_antigens.json
"""
import csv
import json
import re
from collections import defaultdict

IEDB_CSV = 'data/iedb/mhc_ligand_full.csv'
OUTPUT = 'benchmark/dataset/candidate_antigens.json'
MAX_ROWS = 3_000_000  # full file is ~5.7M rows; this sample is sufficient for antigen discovery

COL_PEPTIDE = 11
COL_START = 15
COL_END = 16
COL_SOURCE_MOL = 19
COL_SOURCE_IRI = 20
COL_ORGANISM = 23
COL_QUAL = 94

EXCLUDE_ORGANISMS = {
    'homo sapiens', 'mus musculus', 'rattus norvegicus', 'unidentified',
    'canis lupus familiaris', 'bos taurus', 'gallus gallus', 'mus musculus c57bl/6',
}
EXCLUDE_NAME_TERMS = ['polyprotein', 'hypothetical', 'unknown', 'unnamed', 'conserved hypothetical']


def main():
    antigen_data = defaultdict(lambda: {
        'accession': None, 'organism': '', 'source_mol': '',
        'peptides': [], 'positive_count': 0,
    })

    print(f"Scanning {IEDB_CSV} for antigens with NCBI accessions...")

    with open(IEDB_CSV, 'r', encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f)
        next(reader)
        next(reader)

        count = 0
        for row in reader:
            count += 1
            if count % 500000 == 0:
                print(f"  {count} rows processed, {len(antigen_data)} antigens so far...")
            if count > MAX_ROWS:
                break

            try:
                if len(row) < 108:
                    continue

                peptide = row[COL_PEPTIDE].strip()
                start = row[COL_START].strip()
                end = row[COL_END].strip()
                source_mol = row[COL_SOURCE_MOL].strip()
                source_iri = row[COL_SOURCE_IRI].strip()
                organism = row[COL_ORGANISM].strip()
                qual = row[COL_QUAL].strip().lower()

                if not peptide or not peptide.isalpha() or not peptide.isupper():
                    continue
                if not (8 <= len(peptide) <= 25):
                    continue
                if not start or not end:
                    continue
                if organism.lower() in EXCLUDE_ORGANISMS:
                    continue
                if not source_mol:
                    continue

                m = re.search(r'ncbi\.nlm\.nih\.gov/protein/([A-Za-z0-9_.]+)', source_iri)
                if not m:
                    continue
                accession = m.group(1)

                key = source_mol[:60]
                entry = antigen_data[key]
                entry['accession'] = accession
                entry['organism'] = organism[:40]
                entry['source_mol'] = source_mol[:60]

                is_positive = 'positive' in qual or qual == ''
                if is_positive:
                    entry['positive_count'] += 1
                    entry['peptides'].append((peptide, int(start), int(end)))

            except Exception:
                continue

    print(f"\nTotal rows scanned: {count}")
    print(f"Unique antigens with NCBI accessions: {len(antigen_data)}")

    # Dedupe by accession, keep richest entry, require >=5 positive epitopes
    by_accession = {}
    for key, data in antigen_data.items():
        if data['positive_count'] < 5:
            continue
        acc = data['accession']
        if acc not in by_accession or data['positive_count'] > by_accession[acc]['positive_count']:
            by_accession[acc] = data

    # Drop generic/polyprotein/hypothetical names
    filtered = {
        acc: data for acc, data in by_accession.items()
        if not any(term in data['source_mol'].lower() for term in EXCLUDE_NAME_TERMS)
    }

    sorted_antigens = sorted(filtered.items(), key=lambda x: -x[1]['positive_count'])
    top = dict(sorted_antigens[:120])

    print(f"\nCandidates after filtering: {len(filtered)}, keeping top {len(top)}")

    with open(OUTPUT, 'w') as f:
        json.dump(top, f)
    print(f"Saved to {OUTPUT}")


if __name__ == '__main__':
    main()
