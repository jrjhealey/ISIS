#!/usr/bin/env python3
"""
Select the final benchmark antigen set from validated candidates.

Filters out antigens whose IEDB-confirmed epitope coverage is too extreme
(near 0% or near 100% of the sequence) since those make specificity/
sensitivity undefined or meaningless, then caps per-organism representation
so the final set is diverse rather than dominated by e.g. Vaccinia proteins.

Input:  benchmark/dataset/validated_antigens.json
Output: benchmark/dataset/final_benchmark_set.json
"""
import json
from collections import defaultdict

INPUT = 'benchmark/dataset/validated_antigens.json'
OUTPUT = 'benchmark/dataset/final_benchmark_set.json'

MIN_COVERAGE = 0.03
MAX_COVERAGE = 0.90
MAX_PER_ORGANISM = 6
TARGET_N = 50


def main():
    with open(INPUT) as f:
        validated = json.load(f)

    for v in validated:
        v['coverage'] = len(v['epitope_positions']) / v['seq_len']

    usable = [v for v in validated if MIN_COVERAGE <= v['coverage'] <= MAX_COVERAGE]
    print(f"Usable ({MIN_COVERAGE:.0%}-{MAX_COVERAGE:.0%} coverage): {len(usable)} of {len(validated)}")

    usable.sort(key=lambda x: -x['n_peptides_matched'])

    def select(max_per_organism):
        organism_count = defaultdict(int)
        selection = []
        for v in usable:
            org = v['organism']
            if organism_count[org] < max_per_organism:
                selection.append(v)
                organism_count[org] += 1
        return selection

    final_selection = select(MAX_PER_ORGANISM)
    if len(final_selection) < TARGET_N:
        final_selection = select(10)  # relax cap if still short

    final_selection = final_selection[:TARGET_N]

    print(f"\nFinal selection: {len(final_selection)} antigens")
    organisms = set(v['organism'] for v in final_selection)
    print(f"Unique organisms: {len(organisms)}")
    print(f"Total sequence length: {sum(v['seq_len'] for v in final_selection):,} aa")
    print(f"Total ground-truth epitope positions: {sum(len(v['epitope_positions']) for v in final_selection):,}")

    with open(OUTPUT, 'w') as f:
        json.dump(final_selection, f)
    print(f"\nSaved to {OUTPUT}")


if __name__ == '__main__':
    main()
