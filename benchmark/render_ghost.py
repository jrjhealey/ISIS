#!/usr/bin/env python3
"""
Render ISIS predictions in a "ghost + sketch" style, matching the look of
cryo-EM figure panels: the protein appears as a translucent white silhouette
with pencil-like black outlines, and the predicted epitope spans are picked
out in colour as cartoon, with 3D text labels naming each span.

Run inside ChimeraX:
    chimerax --nogui --offscreen --exit --script benchmark/render_ghost.py

Outputs PNGs to benchmark/output/renders/, four camera aspects per
(molecule, method).
"""
import os
import numpy as np
from chimerax.core.commands import run

OUTDIR = 'benchmark/output/renders'

# Structures to render. All three are real antigens; NP and LF are members of
# the 49-antigen IEDB benchmark set.
MOLECULES = [
    {'path': '/tmp/isis_chimerax_test/6vxx.pdb', 'tag': 'spike',
     'name': 'SARS-CoV-2 spike (6VXX)', 'label_chain': 'A'},
    {'path': '/tmp/isis_chimerax_test/2IQH.pdb', 'tag': 'influenza_np',
     'name': 'Influenza A nucleoprotein (2IQH)', 'label_chain': 'A'},
    {'path': '/tmp/isis_chimerax_test/1J7N.pdb', 'tag': 'anthrax_lf',
     'name': 'Anthrax lethal factor (1J7N)', 'label_chain': 'A'},
]

# One entry per prediction method. `colors` cycles across the spans that
# method finds, echoing the multi-colour subunit palette of the reference
# figure rather than a single hue.
#
# `mode` matters: most methods produce a continuous per-residue score where
# contiguous high-scoring stretches are the feature of interest ('span').
# N-glycosylation is different - it is a point feature marking individual
# sequon asparagines, so a percentile cut is meaningless there (most residues
# score 0, making the 80th percentile 0 and flagging the whole protein).
# Those are found by absolute threshold and drawn as discrete sites.
METHODS = [
    {'tag': 'bcell_consensus', 'cmd': 'isis bcell consensus',
     'attr': 'isis_bcell_consensus', 'mode': 'span', 'min_span': 5,
     'colors': ['#8B3A62', '#5B4B8A', '#1B6B6B', '#3E7B3E', '#6B8E23']},
    {'tag': 'bcell_discotope', 'cmd': 'isis bcell conformational',
     'attr': 'isis_bcell_conf_discotope', 'mode': 'span', 'min_span': 3,
     'colors': ['#1B6B8A', '#2E8B8B', '#3E6B9E', '#4A7BA7', '#5B8FA8']},
    {'tag': 'tcell_mhc1', 'cmd': 'isis tcell mhc1',
     'attr': 'isis_tcell_mhc1', 'mode': 'span', 'min_span': 5,
     'colors': ['#2E7D32', '#388E3C', '#43A047', '#558B2F', '#33691E']},
    {'tag': 'innate_glyco', 'cmd': 'isis innate glyco',
     'attr': 'isis_innate_nglyco', 'mode': 'site', 'threshold': 0.5,
     'colors': ['#8E24AA', '#AB47BC', '#7B1FA2', '#9C27B0', '#6A1B9A']},
    {'tag': 'structure_protrusion', 'cmd': 'isis structure protrusion',
     'attr': 'isis_structure_protrusion', 'mode': 'span', 'min_span': 5,
     'colors': ['#00838F', '#0097A7', '#00ACC1', '#26A69A', '#00695C']},
]

# Camera aspects. Each is a list of commands applied from the reset view.
ASPECTS = [
    ('front', []),
    ('side', ['turn y 90']),
    ('back', ['turn y 180']),
    ('top', ['turn x 90']),
]

PERCENTILE = 80      # residues above this percentile of the score are "of interest"
MAX_SPANS = 5        # label/colour at most this many spans/sites
SITE_PAD = 1         # residues either side of a point site, so it is visible


def apply_ghost_style(session):
    """
    The 'ghost + pencil sketch' look.

    A translucent white molecular surface supplies the cryo-EM-like ghost
    volume; ChimeraX silhouettes draw the dark outline that gives the
    sketched effect; flat lighting keeps it graphic rather than glossy.
    """
    run(session, 'set bgColor white')
    run(session, 'lighting flat')
    run(session, 'lighting shadows false')
    run(session, 'graphics silhouettes true width 2 color black')
    run(session, 'hide atoms')
    run(session, 'show cartoon')
    # Missing-density gaps are drawn as dashed pseudobonds carrying "N residues"
    # labels by default; they clutter the sketch, so drop them.
    run(session, 'hide #1 pseudobonds')
    # Molecular surface (not an EM density map) supplies the ghost volume.
    run(session, 'surface #1')
    run(session, 'color #1 white target s')
    run(session, 'transparency #1 72 s')
    # Base cartoon is a pale grey so the fold reads faintly through the ghost;
    # coloured spans then stand out as the only saturated thing in frame.
    run(session, 'color #1 #B4B4B4 target c')


def find_spans(session, structure, method, chain_id):
    """
    Find the residues of interest for one chain, given a method's config.

    Works off the per-residue attribute the isis command just stored, so it
    reflects exactly what the plugin computed.

    'span' mode uses a percentile cut, because the continuous methods are on
    different scales (vote counts, log-odds, 0-1 probabilities) and no single
    absolute threshold suits them all.

    'site' mode uses an absolute threshold and emits discrete padded sites -
    correct for point features like glycosylation sequons, where a percentile
    cut would flag the entire protein.
    """
    attr = method['attr']
    mode = method.get('mode', 'span')

    for chain in structure.chains:
        if chain.chain_id != chain_id:
            continue

        # chain.residues is index-parallel with chain.characters; None marks
        # positions with no modelled residue (unresolved density).
        vals = []
        for res in chain.residues:
            v = getattr(res, attr, None) if res is not None else None
            vals.append(v)

        present = [v for v in vals if v is not None]
        if not present:
            return []

        if mode == 'site':
            cut = method.get('threshold', 0.5)
            hits = [(v, i) for i, v in enumerate(vals)
                    if v is not None and v >= cut]
            hits.sort(reverse=True)

            out = []
            for score, i in hits[:MAX_SPANS]:
                lo = max(0, i - SITE_PAD)
                hi = min(len(vals) - 1, i + SITE_PAD)
                resnums = [chain.residues[j].number for j in range(lo, hi + 1)
                           if chain.residues[j] is not None]
                site_res = chain.residues[i]
                if resnums and site_res is not None:
                    out.append({'start': min(resnums), 'end': max(resnums),
                                'score': score, 'site': site_res.number})
            return out

        min_span = method.get('min_span', 5)
        if len(present) < min_span:
            return []

        cut = float(np.percentile(present, PERCENTILE))

        spans = []
        start_idx = None
        for i, v in enumerate(vals):
            hot = v is not None and v >= cut
            if hot and start_idx is None:
                start_idx = i
            elif not hot and start_idx is not None:
                if i - start_idx >= min_span:
                    spans.append((start_idx, i - 1))
                start_idx = None
        if start_idx is not None and len(vals) - start_idx >= min_span:
            spans.append((start_idx, len(vals) - 1))

        # Rank spans by mean score, keep the strongest few
        scored = []
        for a, b in spans:
            window = [vals[i] for i in range(a, b + 1) if vals[i] is not None]
            if window:
                scored.append((float(np.mean(window)), a, b))
        scored.sort(reverse=True)

        # Convert sequence indices back to author residue numbers for the
        # colour/label specs, skipping any gap positions at the edges.
        out = []
        for mean_score, a, b in scored[:MAX_SPANS]:
            resnums = [chain.residues[i].number for i in range(a, b + 1)
                       if chain.residues[i] is not None]
            if resnums:
                out.append({'start': min(resnums), 'end': max(resnums),
                            'score': mean_score})
        return out
    return []


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    manifest = []

    for mol in MOLECULES:
        for method in METHODS:
            run(session, 'close all')
            run(session, f"open {mol['path']}")
            structure = session.models[0]

            # Run the prediction, then style. Styling after the command means
            # the command's own auto-colouring is overwritten by the ghost
            # preset, leaving us in control of what gets colour.
            try:
                run(session, f"{method['cmd']} #1")
            except Exception as e:
                print(f"SKIP {mol['tag']}/{method['tag']}: {e}")
                continue

            apply_ghost_style(session)

            spans = find_spans(session, structure, method, mol['label_chain'])
            if not spans:
                print(f"NO SPANS {mol['tag']}/{method['tag']}")

            # Colour the spans on every chain (subunits are equivalent), but
            # label only one chain so the figure stays readable.
            for i, span in enumerate(spans):
                colour = method['colors'][i % len(method['colors'])]
                run(session, f"color #1:{span['start']}-{span['end']} {colour} target c")

            run(session, 'view')
            for i, span in enumerate(spans):
                colour = method['colors'][i % len(method['colors'])]
                if 'site' in span:
                    # Point feature: name the sequon residue itself, e.g. "N165"
                    anchor = span['site']
                    text = f"N{span['site']}"
                else:
                    anchor = (span['start'] + span['end']) // 2
                    text = f"{span['start']}-{span['end']}"
                spec = f"#1/{mol['label_chain']}:{anchor}"
                try:
                    run(session, f"label {spec} residues text \"{text}\" "
                                 f"height 4 color {colour} bgColor white")
                except Exception:
                    pass

            for aspect_name, cmds in ASPECTS:
                run(session, 'view')
                for c in cmds:
                    run(session, c)
                out = f"{OUTDIR}/{mol['tag']}_{method['tag']}_{aspect_name}.png"
                run(session, f"save {out} width 1400 height 1100 supersample 2")
                manifest.append(out)

            print(f"DONE {mol['tag']}/{method['tag']}: {len(spans)} spans -> "
                  + ", ".join(f"{s['start']}-{s['end']}" for s in spans))

    print(f"\nWrote {len(manifest)} images to {OUTDIR}")


main()
