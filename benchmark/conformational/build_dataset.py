#!/usr/bin/env python3
"""
Build a conformational B-cell epitope benchmark from antibody-antigen complexes.

Ground truth here is not a curated annotation but a directly observed physical
contact: an epitope residue is one the antibody is touching in a solved
structure. That is the definition used by the canonical benchmark for this task
(Ponomarenko & Bourne 2007, BMC Struct Biol) - a protein antigen residue with
any atom within 4 A of any antibody atom - and it is what DiscoTope, ElliPro
and SEPPA were themselves evaluated against.

The complex list is Table 2 of that paper, i.e. literature-selected
non-redundant antibody-protein complexes, so the antigens are whole proteins
with experimentally validated epitopes rather than sequence fragments.

One methodological point drives the whole design: structural features are
computed on the antigen chain ALONE, with the antibody deleted. Computing SASA
on the complex would bury exactly the residues that form the epitope, so the
"most buried" residues would be the answer - the feature would encode the label
and any benchmark built on it would be meaningless.

Output: benchmark/conformational/dataset/complexes.json
"""
import json
import os
import subprocess
import time
import warnings

import numpy as np

from Bio.PDB import PDBParser, NeighborSearch, Selection
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1

warnings.filterwarnings("ignore")

OUTDIR = "benchmark/conformational/dataset"
PDBDIR = "/tmp/ab_ag_pdb"
CONTACT_CUTOFF = 4.0      # Angstrom, per Ponomarenko & Bourne
MIN_ANTIGEN_LEN = 50      # a whole protein antigen, not a bound peptide
MIN_EPITOPE_RES = 4

# Table 2 of Ponomarenko & Bourne 2007 - antibody/protein-antigen complexes.
COMPLEXES = """
2ADF 1AFV 1BGX 1E6J 1EGJ 1FSK 1H0D 1I9R 1IQD 1JRH 1LK3 1MHP 1NL0 1NSN 1OAZ
1ORQ 1ORS 1PKQ 1RJL 1SY6 1TZI 1WEJ 1YJD 1YNT 1YY9 1ZA3 1ZTX 2JEL 1A14 1RJC
1JPS 1AR1 1EZV 1OSP 1FNS 1G9M 1R3J 1N8Z 1NFD 1TQB 1TXV 1V7M 1XIW 1Z3G 2AEP
1R0A 1QFU 1EO8 1KEN 2VIR 2VIS 1MLC 1OB1 1P2C 1QFW
""".split()

AB_TERMS = (
    "immunoglobulin", "antibody", "fab ", "fab", "fv ", "scfv", "igg", "ige",
    "iga", "igm", "heavy chain", "light chain", "kappa", "lambda", "nanobody",
    "vhh", "single-domain", "hc fragment", "fragment b", "myeloma",
)
# Descriptions that contain an antibody word but are the antigen itself
AB_EXCLUDE = ("receptor", "binding protein",)

THREE_TO_ONE = {k.upper(): v for k, v in protein_letters_3to1.items()}


def fetch(url, dest=None, retries=3):
    for attempt in range(retries):
        try:
            cmd = ["curl", "-sL", "--max-time", "40", url]
            if dest:
                cmd += ["-o", dest]
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=50)
            if dest:
                if os.path.exists(dest) and os.path.getsize(dest) > 0:
                    return True
            elif r.stdout.strip():
                return r.stdout
        except Exception:
            pass
        time.sleep(1 + attempt)
    return None


def entity_descriptions(pdb_id):
    """Map auth chain id -> (description, length) using the RCSB Data API."""
    entry = fetch(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
    if not entry:
        return {}
    try:
        ids = json.loads(entry)["rcsb_entry_container_identifiers"]["polymer_entity_ids"]
    except Exception:
        return {}

    chains = {}
    for eid in ids:
        raw = fetch(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{eid}")
        time.sleep(0.12)
        if not raw:
            continue
        try:
            d = json.loads(raw)
            desc = (d.get("rcsb_polymer_entity", {}) or {}).get("pdbx_description", "") or ""
            ident = d.get("rcsb_polymer_entity_container_identifiers", {}) or {}
            length = (d.get("entity_poly", {}) or {}).get(
                "rcsb_sample_sequence_length") or 0
            for ch in (ident.get("auth_asym_ids") or []):
                chains[ch] = (desc, int(length or 0))
        except Exception:
            continue
    return chains


def is_antibody(desc):
    d = desc.lower()
    if any(x in d for x in AB_EXCLUDE):
        return False
    return any(t in d for t in AB_TERMS)


def chain_sequence_and_residues(chain):
    """Modelled standard amino acids, in order, with their 1-letter codes."""
    residues, seq = [], []
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        if "CA" not in res:
            continue
        one = THREE_TO_ONE.get(res.get_resname().upper())
        if one is None:
            continue
        residues.append(res)
        seq.append(one)
    return "".join(seq), residues


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(PDBDIR, exist_ok=True)
    parser = PDBParser(QUIET=True)
    sr = ShrakeRupley()

    dataset, skipped = [], []
    unique_ids = sorted(set(COMPLEXES))
    print(f"Processing {len(unique_ids)} antibody-antigen complexes\n")

    for i, pdb_id in enumerate(unique_ids, 1):
        path = os.path.join(PDBDIR, f"{pdb_id}.pdb")
        if not os.path.exists(path):
            if not fetch(f"https://files.rcsb.org/download/{pdb_id}.pdb", dest=path):
                skipped.append((pdb_id, "download failed"))
                print(f"[{i}/{len(unique_ids)}] {pdb_id}: download failed")
                continue

        chains_meta = entity_descriptions(pdb_id)
        if not chains_meta:
            skipped.append((pdb_id, "no entity metadata"))
            print(f"[{i}/{len(unique_ids)}] {pdb_id}: no entity metadata")
            continue

        try:
            model = parser.get_structure(pdb_id, path)[0]
        except Exception as e:
            skipped.append((pdb_id, f"parse error {e}"))
            continue

        ab_chains = [c for c, (d, _l) in chains_meta.items() if is_antibody(d)]
        ag_candidates = [
            (l, c, d) for c, (d, l) in chains_meta.items()
            if not is_antibody(d) and l >= MIN_ANTIGEN_LEN and c in model
        ]
        if not ab_chains or not ag_candidates:
            skipped.append((pdb_id, "no Ab/Ag split"))
            print(f"[{i}/{len(unique_ids)}] {pdb_id}: no Ab/Ag split "
                  f"({len(ab_chains)} Ab, {len(ag_candidates)} Ag)")
            continue

        # Largest non-antibody protein chain is the antigen
        ag_candidates.sort(reverse=True)
        _len, ag_id, ag_desc = ag_candidates[0]

        ab_atoms = [a for c in ab_chains if c in model
                    for a in Selection.unfold_entities(model[c], "A")]
        if not ab_atoms:
            skipped.append((pdb_id, "no antibody atoms"))
            continue

        seq, ag_residues = chain_sequence_and_residues(model[ag_id])
        if len(seq) < MIN_ANTIGEN_LEN:
            skipped.append((pdb_id, f"antigen too short ({len(seq)})"))
            print(f"[{i}/{len(unique_ids)}] {pdb_id}: antigen too short")
            continue

        # --- ground truth: residues touching the antibody (<= 4 A, any atom) ---
        ns = NeighborSearch(ab_atoms)
        epitope_idx = []
        for idx, res in enumerate(ag_residues):
            if any(ns.search(atom.coord, CONTACT_CUTOFF, level="A")
                   for atom in res):
                epitope_idx.append(idx)

        if len(epitope_idx) < MIN_EPITOPE_RES:
            skipped.append((pdb_id, f"epitope too small ({len(epitope_idx)})"))
            print(f"[{i}/{len(unique_ids)}] {pdb_id}: epitope too small")
            continue

        # --- features on the UNBOUND antigen: delete every other chain first ---
        structure = parser.get_structure(pdb_id, path)
        unbound = structure[0]
        for ch in [c.id for c in unbound]:
            if ch != ag_id:
                unbound.detach_child(ch)
        sr.compute(unbound, level="R")

        sasa, coords = [], []
        for res in ag_residues:
            key = res.get_id()
            ures = unbound[ag_id][key] if key in unbound[ag_id] else None
            sasa.append(float(getattr(ures, "sasa", np.nan)) if ures is not None else np.nan)
            coords.append([float(x) for x in res["CA"].coord])

        dataset.append({
            "pdb_id": pdb_id,
            "antigen_chain": ag_id,
            "antigen_desc": ag_desc[:80],
            "antibody_chains": sorted(ab_chains),
            "sequence": seq,
            "n_residues": len(seq),
            "resnums": [r.get_id()[1] for r in ag_residues],
            "epitope_idx": epitope_idx,
            "sasa": sasa,
            "ca_coords": coords,
        })
        print(f"[{i}/{len(unique_ids)}] {pdb_id}: antigen {ag_id} "
              f"({len(seq)} aa) epitope {len(epitope_idx)} res "
              f"({len(epitope_idx)/len(seq):.1%}) | Ab {','.join(sorted(ab_chains))} "
              f"| {ag_desc[:40]}")

    with open(os.path.join(OUTDIR, "complexes.json"), "w") as f:
        json.dump(dataset, f)

    print(f"\n{'='*70}")
    print(f"Built {len(dataset)} complexes, skipped {len(skipped)}")
    if dataset:
        cov = [len(d['epitope_idx']) / d['n_residues'] for d in dataset]
        print(f"Antigen size: {min(d['n_residues'] for d in dataset)}-"
              f"{max(d['n_residues'] for d in dataset)} aa")
        print(f"Epitope coverage: {np.mean(cov):.1%} mean "
              f"({min(cov):.1%}-{max(cov):.1%})")
        print(f"Total residues: {sum(d['n_residues'] for d in dataset):,} | "
              f"epitope residues: {sum(len(d['epitope_idx']) for d in dataset):,}")
    for pid, why in skipped:
        print(f"  skipped {pid}: {why}")


if __name__ == "__main__":
    main()
