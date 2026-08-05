"""
ISIS - In Silico Immunogenicity Studies: ChimeraX commands.

Predicts and visualises immunogenic features directly on an open structure.

B-cell linear (sequence-based)
    isis bcell linear <spec> [method <name>] [window <n>] [threshold <x>]
    isis bcell consensus <spec> [minMethods <n>] [minLength <n>]

B-cell conformational (uses SASA/contacts/protrusion from the structure)
    isis bcell conformational <spec> [method discotope|ellipro|seppa]

T-cell
    isis tcell mhc1 <spec> [allele <name>] [length <n>]
    isis tcell mhc2 <spec> [allele <name>]
    isis tcell proteasome <spec>
    isis tcell tap <spec>
    isis tcell consensus <spec> [allele <name>]

Innate immunity
    isis innate glyco <spec> [glycoType n|o]
    isis innate signal <spec>
    isis innate tlr <spec>
    isis innate consensus <spec>

Structural features
    isis structure sasa <spec>
    isis structure protrusion <spec>
    isis structure contacts <spec> [cutoff <A>]
    isis structure bfactor <spec>

Output and utilities
    isis plot <spec> [method <names>] [window <n>] [minMethods <n>]
              [outdir <path>] [prefix <name>]
    isis export <spec> [format csv|json] [output <path>]
    isis color <spec> [method <name>] [palette <colors>]
    isis clear <spec>
    isis list                 - methods, installed alleles, benchmarked accuracy
    isis doctor               - installation report and how to repair it

Legacy aliases (still supported)
    isis predict <spec> [method <name>]
    isis consensus <spec>
    isis epitopes <spec> [method <name>] [color <name>]

Specifications
--------------
Commands take a residue spec, so a chain selector is honoured: `#1/A` scores
chain A alone, `#1` covers every chain, `#1/A,C` two chains. Colouring, export
and clearing are confined to the chains named.

A spec that clips a chain part-way (`#1/A:100-200`) selects that chain but does
not truncate it - scoring always runs over the chain's full sequence, because a
sliding window over a fragment does not reproduce the scores it would have in
the whole protein.

Installation
------------
The prediction code lives in the separate `isis-epitope` package, which must be
importable from ChimeraX's own Python. If methods report themselves UNAVAILABLE,
run `isis doctor`: it prints the state and the exact command to repair it.

Accuracy
--------
Every class has been benchmarked against experimental data; `isis list` reports
the numbers, including for the components that perform poorly. Briefly: MHC-I
binding is the most reliable (AUC 0.82-0.95 per allele), no conformational
method beats plain solvent accessibility, and four of the five linear B-cell
scales do not exceed chance.
"""

from chimerax.core.commands import CmdDesc, register, StringArg, IntArg, FloatArg, BoolArg
from chimerax.atomic import ResiduesArg

import numpy as np

# The prediction code lives in the separate `isis-epitope` package, which must
# be importable from ChimeraX's OWN bundled Python - not the system Python.
#
# This is the failure mode worth being loud about: the bundle installs fine
# against a stale or absent core, every non-linear method then reports itself
# "not installed", and the listing reads like the tool only supports the five
# original linear scales. It looks like a missing feature rather than a broken
# install. Hence the version check and the explicit remedy below.
MIN_CORE_VERSION = (2, 1)

CORE_VERSION = None
CORE_PROBLEM = None  # human-readable reason, or None when everything is fine

try:
    from isis import predict, predict_all, available_methods, METHOD_INFO
    from isis import __version__ as _core_version_str
    CORE_VERSION = _core_version_str
    ISIS_AVAILABLE = True
except ImportError as _e:
    ISIS_AVAILABLE = False
    CORE_PROBLEM = f"the isis-epitope package is not installed ({_e})"


def _version_tuple(text):
    parts = []
    for chunk in str(text).split("."):
        digits = "".join(c for c in chunk if c.isdigit())
        parts.append(int(digits) if digits else 0)
    return tuple(parts)


try:
    from isis.methods.bcell_conformational import (
        DiscoTopePredictor, ElliProPredictor, SEPPAPredictor,
        calculate_protrusion_index
    )
    from isis.methods.tcell import TcellPredictor, AVAILABLE_ALLELES
    from isis.methods.innate import InnatePredictor
    from isis.methods.structural import (
        calculate_contact_number, calculate_residue_depth,
        normalize_bfactors, calculate_structural_epitope_score,
        detect_surface_patches
    )
    METHODS_AVAILABLE = True
except ImportError as _e:
    METHODS_AVAILABLE = False
    if ISIS_AVAILABLE and CORE_PROBLEM is None:
        # Core imports but the method modules do not: the classic symptom of an
        # old B-cell-linear-only install left behind by an earlier version.
        CORE_PROBLEM = (
            f"isis-epitope {CORE_VERSION} is installed but too old - it has no "
            f"methods module, so only the linear B-cell scales are available "
            f"(need >= {'.'.join(str(v) for v in MIN_CORE_VERSION)})")

if (ISIS_AVAILABLE and CORE_PROBLEM is None
        and _version_tuple(CORE_VERSION) < MIN_CORE_VERSION):
    CORE_PROBLEM = (
        f"isis-epitope {CORE_VERSION} is older than the "
        f"{'.'.join(str(v) for v in MIN_CORE_VERSION)} this bundle expects")

# Plotting is optional - the bundle stays usable for colouring structures even
# if matplotlib is missing from ChimeraX's Python.
try:
    from isis import plotting as isis_plotting
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False


def _chimerax_python():
    """
    Path to ChimeraX's own interpreter, which is where the core must live.

    Derived from sys.prefix rather than sys.executable: inside ChimeraX the
    latter is the launcher binary (/usr/bin/chimerax), so using it would print
    a repair command that cannot work.
    """
    import os as _os
    import sys as _sys
    candidate = _os.path.join(
        _sys.prefix, "bin",
        f"python{_sys.version_info.major}.{_sys.version_info.minor}")
    return candidate if _os.path.exists(candidate) else _sys.executable


def _fix_instructions():
    """The exact commands that repair a missing or stale core install."""
    return [
        "Fix, typed straight into the ChimeraX command line:",
        "    pip install isis-epitope upgrade true",
        "From a clone of the repository instead:",
        "    pip install /path/to/ISIS upgrade true",
        "Or from a shell, using ChimeraX's own interpreter:",
        f'    "{_chimerax_python()}" -m pip install --upgrade isis-epitope',
        "Restart ChimeraX afterwards, then run 'isis doctor' to confirm.",
    ]


def _warn_if_broken(session):
    """Report a missing/stale core loudly, once, with the remedy."""
    if CORE_PROBLEM is None:
        return False
    session.logger.warning(f"ISIS: {CORE_PROBLEM}")
    for line in _fix_instructions():
        session.logger.warning(f"ISIS: {line}")
    return True

ATTR_PREFIX = "isis_"


# =============================================================================
# Helper Functions
# =============================================================================

def _get_ca_coords(residues):
    """Extract CA coordinates from residues (None entries -> NaN row for gaps)."""
    coords = []
    for res in residues:
        ca = res.find_atom("CA") if res is not None else None
        if ca is not None:
            coords.append(ca.coord)
        else:
            coords.append([np.nan, np.nan, np.nan])
    return np.array(coords)


def _get_sasa_values(session, structure, residues):
    """Get SASA values for residues using ChimeraX (None entries -> NaN for gaps)."""
    from chimerax.core.commands import run

    # Run SASA calculation
    run(session, f"measure sasa {structure.atomspec}", log=False)

    # Extract per-residue SASA from atoms (ChimeraX stores per-atom SASA as `.area`)
    sasa_values = []
    for res in residues:
        if res is None:
            sasa_values.append(np.nan)
        else:
            sasa_values.append(sum(getattr(a, 'area', 0) or 0 for a in res.atoms))
    return np.array(sasa_values)


def _get_bfactors(residues):
    """Extract average B-factors per residue (None entries -> NaN for gaps)."""
    bfactors = []
    for res in residues:
        if res is None:
            bfactors.append(np.nan)
            continue
        ca = res.find_atom("CA")
        if ca is not None:
            bfactors.append(ca.bfactor)
        else:
            atoms_bf = [a.bfactor for a in res.atoms if hasattr(a, 'bfactor')]
            bfactors.append(np.mean(atoms_bf) if atoms_bf else np.nan)
    return np.array(bfactors)


def _register_attr(session, attr_name):
    """Register a residue attribute with ChimeraX."""
    from chimerax.atomic import Residue
    if not hasattr(Residue, attr_name):
        Residue.register_attr(session, attr_name, "ISIS", attr_type=float)


def _store_scores(residues, scores, attr_name):
    """Store scores as residue attributes (skips gap positions, res is None)."""
    for i, res in enumerate(residues):
        if res is not None and i < len(scores) and not np.isnan(scores[i]):
            setattr(res, attr_name, float(scores[i]))


def _dense_from_tcell_result(result, seq_len):
    """
    Map a TcellPredictionResult onto a dense per-residue array.

    TcellPredictionResult.peptides is a list of plain peptide strings (not
    dicts), parallel to .scores and .positions (each peptide's 1-indexed
    start position). Each peptide's score is spread (via max) across the
    residue span it covers.
    """
    scores = np.zeros(seq_len)
    positions = result.positions if result.positions else list(range(1, len(result.peptides) + 1))
    for pep, start, score in zip(result.peptides, positions, result.scores):
        end = start + len(pep) - 1
        for pos in range(start, end + 1):
            if 0 <= pos - 1 < seq_len:
                scores[pos - 1] = max(scores[pos - 1], score)
    return scores


def _chain_groups(sel):
    """
    Group a residue selection into per-structure lists of whole chains.

    Yields (structure, [(chain_id, sequence, chain_residues), ...]) for every
    chain the specification touches, in first-seen order.

    Each chain contributes its FULL sequence and full index-parallel residue
    list, not merely the selected slice: a spec that clips a chain part-way
    selects that chain, it does not truncate it. A sliding window run over a
    fragment would not give the same scores as over the whole protein, so
    honouring the clip literally would silently change the numbers.

    This is what makes `#1/A` mean chain A. Passing a whole-structure spec
    (`#1`) still yields every chain, so existing usage is unaffected.
    """
    order = []
    chain_ids = {}
    for res in sel:
        st = res.structure
        if st not in chain_ids:
            chain_ids[st] = []
            order.append(st)
        if res.chain_id not in chain_ids[st]:
            chain_ids[st].append(res.chain_id)

    for st in order:
        by_id = {c.chain_id: c for c in st.chains}
        chains = []
        for cid in chain_ids[st]:
            chain = by_id.get(cid)
            if chain is None:
                continue
            seq = chain.characters
            if seq and len(seq) >= 6:
                chains.append((cid, seq, chain.residues))
        if chains:
            yield st, chains


def _chain_residues(chains):
    """Flat list of modelled residues across the given chains (skips gaps)."""
    out = []
    for _cid, _seq, residues in chains:
        out.extend(r for r in residues if r is not None)
    return out


def _chain_spec(structure, chains):
    """Atom spec covering just the given chains of a structure."""
    if not chains:
        return structure.atomspec
    ids = ",".join(cid for cid, _seq, _res in chains)
    return f"{structure.atomspec}/{ids}"


def _auto_color(session, structure, attr_name, palette="white:yellow:red",
                chains=None):
    """
    Colour by attribute, restricted to the chains that were selected.

    Scoping matters: colouring the whole structure when only one chain was
    asked for would repaint the other chains using their (absent) attribute
    values, visibly changing subunits the user did not ask about.
    """
    from chimerax.core.commands import run
    spec = _chain_spec(structure, chains)
    run(session, f"color byattribute {attr_name} {spec} palette {palette}")


def _find_epitope_regions(scores, threshold, min_length=6):
    """Find contiguous regions above threshold."""
    regions = []
    in_region = False
    start = 0
    region_scores = []

    for i, score in enumerate(scores):
        pos = i + 1
        if score >= threshold:
            if not in_region:
                in_region = True
                start = pos
                region_scores = [score]
            else:
                region_scores.append(score)
        else:
            if in_region:
                in_region = False
                length = (pos - 1) - start + 1
                if length >= min_length:
                    avg_score = sum(region_scores) / len(region_scores)
                    regions.append((start, pos - 1, avg_score))
                region_scores = []

    if in_region:
        length = len(scores) - start + 1
        if length >= min_length:
            avg_score = sum(region_scores) / len(region_scores)
            regions.append((start, len(scores), avg_score))

    return regions


# =============================================================================
# B-cell Linear Epitopes (existing methods)
# =============================================================================

def isis_bcell_linear(session, sel, method="emini", window=None, threshold=None):
    """Run B-cell linear epitope prediction."""
    if not ISIS_AVAILABLE:
        session.logger.error("ISIS library not installed")
        return

    from chimerax.atomic import Residue

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"B-cell linear prediction ({method}) on {structure.name}...")
        attr_name = f"{ATTR_PREFIX}bcell_{method.replace('-', '_')}"
        _register_attr(session, attr_name)

        for chain_id, seq, residues in _chains:
            try:
                pred = predict(seq, method=method, window_size=window, threshold=threshold)
                scores = np.array([pred.score_at(i+1) or 0 for i in range(len(seq))])
                _store_scores(residues, scores, attr_name)

                session.logger.info(f"  Chain {chain_id}: {len(pred.epitopes)} epitopes")
                for ep in pred.epitopes:
                    session.logger.info(f"    {ep.start}-{ep.end}: {ep.sequence}")
            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, chains=_chains)
        session.logger.info(f"Scores stored as: {attr_name}")


def isis_bcell_conformational(session, sel, method="discotope"):
    """Run B-cell conformational epitope prediction (requires 3D structure)."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"B-cell conformational prediction ({method}) on {structure.name}...")
        attr_name = f"{ATTR_PREFIX}bcell_conf_{method}"
        _register_attr(session, attr_name)

        for chain_id, seq, residues in _chains:
            try:
                coords = _get_ca_coords(residues)
                sasa = _get_sasa_values(session, structure, residues)
                contacts = calculate_contact_number(coords)

                if method == "discotope":
                    predictor = DiscoTopePredictor()
                    result = predictor.predict_with_structure(seq, sasa, contacts, coords)
                elif method == "ellipro":
                    predictor = ElliProPredictor()
                    result = predictor.predict_with_structure(seq, sasa, contacts, coords)
                elif method == "seppa":
                    predictor = SEPPAPredictor()
                    result = predictor.predict_with_structure(seq, sasa, contacts, coords)
                else:
                    session.logger.error(f"Unknown method: {method}")
                    continue

                _store_scores(residues, result['scores'], attr_name)

                epitopes = result.get('epitopes', [])
                session.logger.info(f"  Chain {chain_id}: {len(epitopes)} epitopes")
                for ep in epitopes:
                    # ConformationalEpitope.residue_indices are 0-indexed and
                    # may be non-contiguous (spatial, not sequence, clusters)
                    span = f"{min(ep.residue_indices)+1}-{max(ep.residue_indices)+1}"
                    session.logger.info(f"    {span} ({ep.size} residues): {ep.residues} score={ep.score:.2f}")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")
                import traceback
                session.logger.error(traceback.format_exc())

        _auto_color(session, structure, attr_name, palette="white:cyan:blue", chains=_chains)
        session.logger.info(f"Scores stored as: {attr_name}")


def isis_bcell_consensus(session, sel, min_methods=2, min_length=6):
    """Consensus across all B-cell methods (linear + conformational)."""
    if not ISIS_AVAILABLE:
        session.logger.error("ISIS library not installed")
        return

    linear_methods = available_methods()
    conf_methods = ["discotope", "ellipro"] if METHODS_AVAILABLE else []

    attr_name = f"{ATTR_PREFIX}bcell_consensus"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"B-cell consensus on {structure.name}...")

        for chain_id, seq, residues in _chains:
            votes = np.zeros(len(seq))
            methods_run = 0

            # Linear methods
            for method in linear_methods:
                try:
                    pred = predict(seq, method=method)
                    for ep in pred.epitopes:
                        if ep.length >= min_length:
                            for pos in range(ep.start, ep.end + 1):
                                if 1 <= pos <= len(seq):
                                    votes[pos - 1] += 1
                    methods_run += 1
                except:
                    pass

            # Conformational methods
            if METHODS_AVAILABLE:
                try:
                    coords = _get_ca_coords(residues)
                    sasa = _get_sasa_values(session, structure, residues)
                    contacts = calculate_contact_number(coords)

                    for method in conf_methods:
                        try:
                            if method == "discotope":
                                predictor = DiscoTopePredictor()
                            else:
                                predictor = ElliProPredictor()
                            result = predictor.predict_with_structure(seq, sasa, contacts, coords)
                            for ep in result.get('epitopes', []):
                                # residue_indices are 0-indexed and may be
                                # non-contiguous - vote at exactly those
                                # positions, not a start-end range
                                for idx in ep.residue_indices:
                                    if 0 <= idx < len(seq):
                                        votes[idx] += 1
                            methods_run += 1
                        except Exception as e:
                            session.logger.warning(f"  Chain {chain_id}: {method} failed: {e}")
                except Exception as e:
                    session.logger.warning(f"  Chain {chain_id}: conformational setup failed: {e}")

            _store_scores(residues, votes, attr_name)

            # Report consensus regions
            regions = _find_epitope_regions(votes, min_methods, min_length)
            session.logger.info(f"  Chain {chain_id}: {len(regions)} consensus epitopes ({methods_run} methods)")
            for start, end, avg in regions:
                peptide = seq[start-1:end]
                session.logger.info(f"    {start}-{end}: {peptide} (avg {avg:.1f} methods)")

        _auto_color(session, structure, attr_name, palette="white:yellow:orange:red", chains=_chains)


# =============================================================================
# T-cell Epitopes
# =============================================================================

def isis_tcell_mhc1(session, sel, allele="HLA-A*02:01", length=9):
    """Predict MHC Class I binding."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = TcellPredictor()
    attr_name = f"{ATTR_PREFIX}tcell_mhc1"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"MHC-I prediction ({allele}) on {structure.name}...")

        for chain_id, seq, residues in _chains:
            try:
                result = predictor.predict_mhc1(seq, allele=allele, peptide_length=length)

                scores = _dense_from_tcell_result(result, len(seq))
                _store_scores(residues, scores, attr_name)

                # Report top binders
                binders = result.get_binders()
                session.logger.info(f"  Chain {chain_id}: {len(binders)} top binders (>= threshold)")
                for pep in binders[:5]:
                    end = pep['position'] + len(pep['peptide']) - 1
                    ic50 = pep.get('ic50', 'N/A')
                    session.logger.info(f"    {pep['position']}-{end}: {pep['peptide']} (IC50={ic50})")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="white:green:darkgreen", chains=_chains)


def isis_tcell_mhc2(session, sel, allele="HLA-DRB1*01:01"):
    """Predict MHC Class II binding."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = TcellPredictor()
    attr_name = f"{ATTR_PREFIX}tcell_mhc2"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"MHC-II prediction ({allele}) on {structure.name}...")

        for chain_id, seq, residues in _chains:
            try:
                result = predictor.predict_mhc2(seq, allele=allele)

                scores = _dense_from_tcell_result(result, len(seq))
                _store_scores(residues, scores, attr_name)

                binders = result.get_binders()
                session.logger.info(f"  Chain {chain_id}: {len(binders)} top binders")
                for pep in binders[:5]:
                    end = pep['position'] + len(pep['peptide']) - 1
                    session.logger.info(f"    {pep['position']}-{end}: {pep['peptide']}")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="white:purple:darkviolet", chains=_chains)


def isis_tcell_proteasome(session, sel):
    """Predict proteasomal cleavage sites."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = TcellPredictor()
    attr_name = f"{ATTR_PREFIX}tcell_cleavage"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"Proteasomal cleavage prediction on {structure.name}...")

        for chain_id, seq, residues in _chains:
            try:
                result = predictor.predict_cleavage(seq)
                _store_scores(residues, result.scores, attr_name)

                # Report top cleavage sites
                top_sites = sorted(enumerate(result.scores), key=lambda x: -x[1])[:10]
                session.logger.info(f"  Chain {chain_id}: Top cleavage sites:")
                for pos, score in top_sites[:5]:
                    session.logger.info(f"    Position {pos+1}: {seq[pos]} (score={score:.2f})")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="white:orange:red", chains=_chains)


def isis_tcell_tap(session, sel):
    """Predict TAP transport efficiency."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = TcellPredictor()
    attr_name = f"{ATTR_PREFIX}tcell_tap"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"TAP transport prediction on {structure.name}...")

        for chain_id, seq, residues in _chains:
            try:
                result = predictor.predict_tap(seq)

                scores = _dense_from_tcell_result(result, len(seq))
                _store_scores(residues, scores, attr_name)
                session.logger.info(f"  Chain {chain_id}: TAP scores computed")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="white:lightblue:blue", chains=_chains)


def isis_tcell_consensus(session, sel, allele="HLA-A*02:01"):
    """Combined T-cell epitope prediction pipeline."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = TcellPredictor()
    attr_name = f"{ATTR_PREFIX}tcell_consensus"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"T-cell consensus prediction on {structure.name}...")

        for chain_id, seq, residues in _chains:
            try:
                result = predictor.predict_pipeline(seq, allele=allele)

                scores = _dense_from_tcell_result(result, len(seq))
                _store_scores(residues, scores, attr_name)

                # Report top candidates (result.get_binders() sorts by score, descending)
                top = result.get_binders()[:10]
                session.logger.info(f"  Chain {chain_id}: Top T-cell epitope candidates:")
                for pep in top[:5]:
                    end = pep['position'] + len(pep['peptide']) - 1
                    session.logger.info(f"    {pep['position']}-{end}: {pep['peptide']} (score={pep['score']:.2f})")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="white:lime:green", chains=_chains)


# =============================================================================
# Innate Immunity
# =============================================================================

def isis_innate_glyco(session, sel, glyco_type="n"):
    """Predict glycosylation sites."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = InnatePredictor()
    attr_name = f"{ATTR_PREFIX}innate_{glyco_type}glyco"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"{glyco_type.upper()}-glycosylation prediction on {structure.name}...")

        for chain_id, seq, residues in _chains:
            try:
                if glyco_type == "n":
                    result = predictor.predict_n_glyco(seq)
                else:
                    result = predictor.predict_o_glyco(seq)

                # predict_n_glyco/predict_o_glyco return sparse per-site lists
                # (one entry per detected site), not one score per residue -
                # expand into a dense per-residue array before storing.
                dense_scores = np.zeros(len(seq))
                for pos, score in zip(result['positions'], result['scores']):
                    if 1 <= pos <= len(seq):
                        dense_scores[pos - 1] = score
                _store_scores(residues, dense_scores, attr_name)

                # Report sites (result['sites'] holds the full dict per site)
                sites = result.get('sites', [])
                high_conf = [s for s in sites if s['score'] > 0.5]
                session.logger.info(f"  Chain {chain_id}: {len(high_conf)} potential sites")
                for site in sites[:10]:
                    session.logger.info(f"    Position {site['position']}: {site['motif']} ({site.get('site_type', 'canonical')})")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="white:pink:magenta", chains=_chains)


def isis_innate_signal(session, sel):
    """Predict signal peptide."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = InnatePredictor()
    attr_name = f"{ATTR_PREFIX}innate_signal"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"Signal peptide prediction on {structure.name}...")

        for chain_id, seq, residues in _chains:
            try:
                result = predictor.predict_signal_peptide(seq)

                # Create score array marking signal peptide region
                scores = np.zeros(len(seq))
                if result['is_signal_peptide']:
                    cleavage = result['cleavage_position']
                    for i in range(min(cleavage, len(scores))):
                        scores[i] = result['confidence']

                _store_scores(residues, scores, attr_name)

                if result['is_signal_peptide']:
                    session.logger.info(f"  Chain {chain_id}: Signal peptide detected")
                    session.logger.info(f"    Cleavage at position {result['cleavage_position']}")
                    session.logger.info(f"    Confidence: {result['confidence']:.2f}")
                else:
                    session.logger.info(f"  Chain {chain_id}: No signal peptide detected")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="white:yellow:orange", chains=_chains)


def isis_innate_tlr(session, sel):
    """Predict TLR ligand motifs."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = InnatePredictor()
    attr_name = f"{ATTR_PREFIX}innate_tlr"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"TLR motif prediction on {structure.name}...")

        for chain_id, seq, residues in _chains:
            try:
                result = predictor.predict_tlr_motifs(seq)

                # Score positions by TLR relevance (result['matches'] holds the
                # full per-match dicts; result['motifs'] is just matched sequence text)
                matches = result.get('matches', [])
                scores = np.zeros(len(seq))
                for match in matches:
                    start, end = match['position'], match['end']
                    for pos in range(start, end + 1):
                        if 0 <= pos - 1 < len(scores):
                            scores[pos - 1] = max(scores[pos - 1], match.get('score', 1.0))

                _store_scores(residues, scores, attr_name)

                session.logger.info(f"  Chain {chain_id}: {len(matches)} TLR-relevant motifs")
                for match in matches[:5]:
                    session.logger.info(f"    {match['motif_name']}: positions {match['position']}-{match['end']}")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="white:coral:red", chains=_chains)


def isis_innate_consensus(session, sel):
    """Consensus across all innate immunity predictions."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = InnatePredictor()
    attr_name = f"{ATTR_PREFIX}innate_consensus"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"Innate immunity consensus on {structure.name}...")

        for chain_id, seq, residues in _chains:
            try:
                # Combine all innate predictions
                scores = np.zeros(len(seq))

                # N-glycosylation (sparse per-site results -> dense per-residue)
                n_glyco = predictor.predict_n_glyco(seq)
                for pos, score in zip(n_glyco['positions'], n_glyco['scores']):
                    if 1 <= pos <= len(seq):
                        scores[pos - 1] += 0.3 * score

                # O-glycosylation (sparse per-site results -> dense per-residue)
                o_glyco = predictor.predict_o_glyco(seq)
                for pos, score in zip(o_glyco['positions'], o_glyco['scores']):
                    if 1 <= pos <= len(seq):
                        scores[pos - 1] += 0.3 * score

                # TLR motifs (result['matches'] holds the full per-match dicts)
                tlr = predictor.predict_tlr_motifs(seq)
                for match in tlr.get('matches', []):
                    for pos in range(match['position'], match['end'] + 1):
                        if 0 <= pos - 1 < len(scores):
                            scores[pos - 1] += 0.4 * match.get('score', 1.0)

                _store_scores(residues, scores, attr_name)
                session.logger.info(f"  Chain {chain_id}: Innate consensus computed")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="white:salmon:darkred", chains=_chains)


# =============================================================================
# Structural Analysis
# =============================================================================

def isis_structure_sasa(session, sel):
    """Calculate and display solvent accessible surface area."""
    attr_name = f"{ATTR_PREFIX}structure_sasa"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"SASA calculation on {structure.name}...")

        for chain_id, seq, residues in _chains:
            try:
                sasa = _get_sasa_values(session, structure, residues)
                _store_scores(residues, sasa, attr_name)

                avg_sasa = np.mean(sasa)
                session.logger.info(f"  Chain {chain_id}: avg SASA = {avg_sasa:.1f} Å²")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="blue:white:red", chains=_chains)


def isis_structure_protrusion(session, sel):
    """Calculate protrusion index."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    attr_name = f"{ATTR_PREFIX}structure_protrusion"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"Protrusion index on {structure.name}...")

        for chain_id, seq, residues in _chains:
            try:
                coords = _get_ca_coords(residues)
                protrusion = calculate_protrusion_index(coords)
                _store_scores(residues, protrusion, attr_name)

                session.logger.info(f"  Chain {chain_id}: protrusion index computed")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="white:cyan:blue", chains=_chains)


def isis_structure_contacts(session, sel, cutoff=8.0):
    """Calculate contact numbers."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    attr_name = f"{ATTR_PREFIX}structure_contacts"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"Contact number calculation on {structure.name}...")

        for chain_id, seq, residues in _chains:
            try:
                coords = _get_ca_coords(residues)
                contacts = calculate_contact_number(coords, cutoff=cutoff)
                _store_scores(residues, contacts, attr_name)

                avg_contacts = np.mean(contacts)
                session.logger.info(f"  Chain {chain_id}: avg contacts = {avg_contacts:.1f}")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="red:white:blue", chains=_chains)


def isis_structure_bfactor(session, sel):
    """Display normalized B-factors."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    attr_name = f"{ATTR_PREFIX}structure_bfactor"
    _register_attr(session, attr_name)

    for structure, _chains in _chain_groups(sel):
        session.logger.info(f"B-factor analysis on {structure.name}...")

        for chain_id, seq, residues in _chains:
            try:
                bfactors = _get_bfactors(residues)
                normalized = normalize_bfactors(bfactors)
                _store_scores(residues, normalized, attr_name)

                session.logger.info(f"  Chain {chain_id}: B-factors normalized")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="blue:white:red", chains=_chains)


# =============================================================================
# Export and Utility Commands
# =============================================================================

def isis_export(session, sel, format="csv", output=None):
    """Export all ISIS predictions to file."""
    import json
    import csv
    from io import StringIO

    for structure, _chains in _chain_groups(sel):
        data = {"structure": structure.name, "chains": []}

        for chain_id, seq, residues in _chains:
            chain_data = {
                "chain_id": chain_id,
                "sequence": seq,
                "predictions": {}
            }

            # First pass: discover which ISIS attributes are present anywhere
            # in this chain (gap/None positions never carry attributes).
            attr_names = set()
            for res in residues:
                if res is not None:
                    attr_names.update(a for a in dir(res) if a.startswith(ATTR_PREFIX))

            # Second pass: build one dense, position-aligned list per attribute
            # so index i always corresponds to seq[i], even across gaps.
            for attr in attr_names:
                values = []
                for res in residues:
                    val = getattr(res, attr, None) if res is not None else None
                    values.append(val if val else 0)
                chain_data["predictions"][attr] = values

            data["chains"].append(chain_data)

        if format == "json":
            result = json.dumps(data, indent=2)
        else:  # csv
            output_lines = []
            for chain in data["chains"]:
                output_lines.append(f"# Chain {chain['chain_id']}")
                headers = ["position", "residue"] + list(chain["predictions"].keys())
                output_lines.append(",".join(headers))
                for i, aa in enumerate(chain["sequence"]):
                    row = [str(i+1), aa]
                    for attr in chain["predictions"]:
                        scores = chain["predictions"][attr]
                        row.append(f"{scores[i]:.4f}" if i < len(scores) else "")
                    output_lines.append(",".join(row))
            result = "\n".join(output_lines)

        if output:
            with open(output, 'w') as f:
                f.write(result)
            session.logger.info(f"Exported to {output}")
        else:
            session.logger.info(f"Export ({format}):\n{result[:500]}...")


def isis_color(session, sel, method=None, palette="white:yellow:red"):
    """Color structure by ISIS prediction scores."""
    from chimerax.core.commands import run

    # If no method specified, try to find any ISIS attribute
    if method is None:
        method = "consensus"

    attr_name = f"{ATTR_PREFIX}{method.replace('-', '_')}"

    # Also try category-specific names
    possible_attrs = [
        attr_name,
        f"{ATTR_PREFIX}bcell_{method.replace('-', '_')}",
        f"{ATTR_PREFIX}tcell_{method}",
        f"{ATTR_PREFIX}innate_{method}",
        f"{ATTR_PREFIX}structure_{method}",
    ]

    for structure, _chains in _chain_groups(sel):
        # Look for the attribute only on the chains that were selected: a
        # prediction stored on chain B is not evidence that chain A can be
        # coloured by it.
        candidates = _chain_residues(_chains)
        found_attr = None
        for attr in possible_attrs:
            if any(hasattr(res, attr) for res in candidates):
                found_attr = attr
                break

        if not found_attr:
            session.logger.warning(f"{structure.name}: No predictions found for '{method}'")
            continue

        _auto_color(session, structure, found_attr, palette=palette, chains=_chains)
        session.logger.info(f"Colored {structure.name} by {found_attr}")


def isis_plot(session, sel, method=None, window=None, minMethods=3,
              outdir=None, prefix=None):
    """
    Render prediction figures straight from an open structure.

    No FASTA needed: the sequence comes from the structure's own chain
    sequence, which ChimeraX has already reconciled against the model.
    """
    if not ISIS_AVAILABLE:
        session.logger.error("ISIS library not installed")
        return
    if not PLOTTING_AVAILABLE:
        session.logger.error(
            "Plotting needs matplotlib in ChimeraX's Python. Install with:\n"
            "    pip install matplotlib      (in the ChimeraX command line)")
        return

    import os

    methods = ([m.strip() for m in method.split(",")] if method
               else list(available_methods()))

    outdir = outdir or os.getcwd()
    os.makedirs(outdir, exist_ok=True)

    plotted = 0
    for structure, chains in _chain_groups(sel):
        for chain_id, seq, _chain_res in chains:
            try:
                results = {m: predict(seq, method=m, window_size=window)
                           for m in methods}
            except Exception as e:
                session.logger.error(f"{structure.name}/{chain_id}: {e}")
                continue

            base = prefix or os.path.splitext(structure.name)[0]
            slug = "".join(c if c.isalnum() or c in "-_" else "_"
                           for c in f"{base}_{chain_id}")
            label = f"{structure.name} chain {chain_id}"

            written = [
                isis_plotting.save_figure(
                    isis_plotting.plot_profile(
                        seq, results,
                        subtitle=f"{label} · {len(seq)} residues"),
                    os.path.join(outdir, f"{slug}_profile.png")),
                isis_plotting.save_figure(
                    isis_plotting.plot_call_matrix(seq, results),
                    os.path.join(outdir, f"{slug}_calls.png")),
            ]

            votes = np.zeros(len(seq))
            for pred_result in results.values():
                for ep in pred_result.epitopes:
                    votes[ep.start - 1:ep.end] += 1
            written.append(isis_plotting.save_figure(
                isis_plotting.plot_consensus(seq, votes,
                                             min_methods=minMethods,
                                             n_methods=len(results)),
                os.path.join(outdir, f"{slug}_consensus.png")))

            plotted += 1
            session.logger.info(
                f"{structure.name}/{chain_id}: {len(seq)} residues, "
                f"{len(methods)} methods")
            for path in written:
                # Clickable in the ChimeraX log, so the figure is one click away.
                session.logger.info(
                    f'&nbsp;&nbsp;<a href="file://{path}">{os.path.basename(path)}</a>',
                    is_html=True)

    if plotted == 0:
        session.logger.warning("No polymer chains in the given specification")


def isis_clear(session, sel):
    """Remove all ISIS prediction attributes."""
    for structure, _chains in _chain_groups(sel):
        # Clear only the selected chains, so `isis clear #1/A` leaves the
        # predictions on B and C intact.
        cleared = set()
        for res in _chain_residues(_chains):
            for attr in list(vars(res).keys()):
                if attr.startswith(ATTR_PREFIX):
                    delattr(res, attr)
                    cleared.add(attr)
        chain_list = ",".join(cid for cid, _s, _r in _chains)
        if cleared:
            session.logger.info(
                f"{structure.name}/{chain_list}: Cleared {len(cleared)} attributes")
        else:
            session.logger.info(
                f"{structure.name}/{chain_list}: No ISIS attributes found")


def isis_doctor(session):
    """
    Report the installation state and how to repair it.

    Exists because the common failure is silent: the bundle loads, the linear
    scales work, and everything else reports itself missing, which reads as
    "this tool only does B-cell linear" rather than "the library is stale".
    """
    import sys as _sys
    log = session.logger.info

    log("ISIS installation report")
    log("")
    log(f"  ChimeraX Python : {_chimerax_python()}")
    log(f"  Launcher        : {_sys.executable}")
    log(f"  Python version  : {_sys.version.split()[0]}")
    # Where this module was actually imported from, plus every other copy on
    # sys.path. ChimeraX searches its user data directory BEFORE the application
    # directory, so a bundle installed there by an earlier `devel install`
    # shadows a newer one installed into the app - the new code is on disk but
    # never loaded, so new commands are "Unknown command" and the listing stays
    # old. Reporting the load path makes that visible instead of baffling.
    import os as _os
    log(f"  Bundle loaded   : {_os.path.dirname(__file__)}")
    _copies = []
    for _p in _sys.path:
        _cand = _os.path.join(_p, "chimerax", "isis")
        if _os.path.isdir(_cand):
            _real = _os.path.realpath(_cand)
            if _real not in _copies:
                _copies.append(_real)
    if len(_copies) > 1:
        session.logger.warning(
            f"  {len(_copies)} copies of the bundle are installed. The first on "
            f"sys.path wins and shadows the others:")
        for _c in _copies:
            session.logger.warning(f"    {_c}")
        session.logger.warning(
            "  Re-run install_chimerax.sh, which removes every copy before "
            "installing.")
    log(f"  Core library    : isis-epitope {CORE_VERSION or 'NOT FOUND'}")
    log(f"  Expected core   : >= {'.'.join(str(v) for v in MIN_CORE_VERSION)}")
    log("")

    log("  Component status")
    rows = [
        ("B-cell linear scales", ISIS_AVAILABLE),
        ("B-cell conformational", METHODS_AVAILABLE),
        ("T-cell / MHC", METHODS_AVAILABLE),
        ("Innate immunity", METHODS_AVAILABLE),
        ("Structural features", True),
        ("Plotting (matplotlib)", PLOTTING_AVAILABLE),
    ]
    for label, ok in rows:
        log(f"    {'OK     ' if ok else 'MISSING'}  {label}")

    if METHODS_AVAILABLE:
        try:
            from isis.models.mhc_predictor import MHCPredictor
            mp = MHCPredictor()
            log(f"    OK       MHC models: {len(mp.available_alleles(1))} class I, "
                f"{len(mp.available_alleles(2))} class II alleles")
        except Exception as e:
            log(f"    MISSING  MHC models ({e})")

    log("")
    if CORE_PROBLEM is None and METHODS_AVAILABLE:
        log("  Everything required is present. Run 'isis list' for the methods.")
    else:
        if CORE_PROBLEM:
            session.logger.warning(f"  Problem: {CORE_PROBLEM}")
        for line in _fix_instructions():
            log(f"  {line}")


def isis_list(session):
    """
    List available prediction methods and the alleles actually installed.

    Allele lists are queried from the loaded model files rather than written
    out here. A hard-coded list drifts the moment models are retrained - this
    one previously advertised HLA-DRB1*04:01, which was never shipped, while
    omitting six MHC-I alleles that were.
    """
    log = session.logger.info

    log(f"ISIS Prediction Methods  (core library isis-epitope "
        f"{CORE_VERSION or 'NOT INSTALLED'})")
    log("")
    if _warn_if_broken(session):
        log("")
        log("Methods shown UNAVAILABLE below are a broken install, not")
        log("missing features. Run 'isis doctor' for the full report.")
        log("")

    log("B-cell Linear (sequence-based) - isis bcell linear:")
    if ISIS_AVAILABLE:
        for method in available_methods():
            info = METHOD_INFO.get(method, {})
            log(f"  {method}: {info.get('name', method)}")
    else:
        log("  UNAVAILABLE - isis-epitope missing; run 'isis doctor'")

    log("")
    log("B-cell Conformational (needs 3D structure) - isis bcell conformational:")
    if METHODS_AVAILABLE:
        for key, predictor in (("discotope", DiscoTopePredictor),
                               ("ellipro", ElliProPredictor),
                               ("seppa", SEPPAPredictor)):
            log(f"  {key}: {getattr(predictor, 'description', key)}")
    else:
        log("  UNAVAILABLE - core library missing or too old; run 'isis doctor'")

    log("")
    log("T-cell (MHC binding) - isis tcell mhc1 / mhc2:")
    if METHODS_AVAILABLE:
        try:
            from isis.models.mhc_predictor import MHCPredictor
            mp = MHCPredictor()
            mhc1 = mp.available_alleles(1)
            mhc2 = mp.available_alleles(2)
            log(f"  mhc1: {len(mhc1)} alleles installed")
            for a in mhc1:
                log(f"    {a}")
            log(f"  mhc2: {len(mhc2)} alleles installed")
            for a in mhc2:
                log(f"    {a}")
        except Exception as e:
            log(f"  (MHC models unavailable: {e})")
        log("  proteasome: Proteasomal cleavage sites")
        log("  tap: TAP transport efficiency")
        log("  consensus: combined MHC + cleavage + TAP pipeline")
    else:
        log("  UNAVAILABLE - core library missing or too old; run 'isis doctor'")

    log("")
    log("Innate Immunity - isis innate:")
    if METHODS_AVAILABLE:
        log("  glyco: N- and O-glycosylation sequons")
        log("  signal: Signal peptide detection")
        log("  tlr: TLR ligand motifs")
        log("  consensus: weighted combination of the above")
    else:
        log("  UNAVAILABLE - core library missing or too old; run 'isis doctor'")

    log("")
    log("Structural Analysis - isis structure:")
    log("  sasa: Solvent accessible surface area")
    log("  protrusion: Protrusion index")
    log("  contacts: Residue contact numbers")
    log("  bfactor: Normalized B-factors")

    log("")
    log("Output and utilities:")
    log("  isis plot: per-residue figures from the structure's own sequence"
        + ("" if PLOTTING_AVAILABLE else "  (matplotlib missing)"))
    log("  isis export: scores to CSV/JSON")
    log("  isis color / isis clear: recolour or remove stored scores")

    log("")
    log("Specs: any command accepts a residue spec, so #1/A means chain A only.")
    log("A spec that clips a chain selects it without truncating its sequence.")

    log("")
    log("Run 'isis doctor' to check the installation.")
    log("")
    log("Benchmarked accuracy (see benchmark/ in the repository):")
    log("  MHC-I binding      AUC 0.82-0.95 per allele (mean 0.87) - most reliable")
    log("  Conformational     AUC 0.61-0.68; none beats plain SASA (0.67)")
    log("  B-cell linear      only kolaskar-tongaonkar exceeds chance on linear")
    log("                     IEDB epitopes (MCC +0.15); it is BELOW chance for")
    log("                     structural contact patches - the targets differ")
    log("  Signal peptide     4/6 on positive/negative controls - weak")
    log("  O-glycosylation    flags 58-83% of all S/T - not discriminative")
    log("  proteasome, tap, tlr: not yet benchmarked")


# =============================================================================
# Legacy Commands (backward compatibility)
# =============================================================================

def isis_predict(session, sel, method="emini", window=None, threshold=None):
    """Legacy predict command - redirects to isis_bcell_linear."""
    isis_bcell_linear(session, sel, method=method, window=window, threshold=threshold)


def isis_consensus(session, sel, min_methods=2, min_length=6):
    """Legacy consensus command - redirects to isis_bcell_consensus."""
    isis_bcell_consensus(session, sel, min_methods=min_methods, min_length=min_length)


def isis_epitopes(session, sel, method="emini", color="red"):
    """Color only predicted epitope regions."""
    if not ISIS_AVAILABLE:
        session.logger.error("ISIS library not installed")
        return

    from chimerax.core.commands import run

    for structure, _chains in _chain_groups(sel):
        for chain_id, seq, residues in _chains:
            try:
                pred = predict(seq, method=method)
                for epitope in pred.epitopes:
                    start, end = epitope.start, epitope.end
                    spec = f"{structure.atomspec}/{chain_id}:{start}-{end}"
                    run(session, f"color {spec} {color}")
            except Exception as e:
                session.logger.error(f"Chain {chain_id}: {e}")


# =============================================================================
# Command Registration
# =============================================================================

def register_all_commands(session):
    """Register all ISIS commands."""
    from chimerax.core.commands import register as reg

    try:
        session.logger.info("ISIS: Registering commands...")

        # B-cell commands
        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("method", StringArg), ("window", IntArg), ("threshold", FloatArg)],
            synopsis="B-cell linear epitope prediction"
        )
        reg("isis bcell linear", desc, isis_bcell_linear)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("method", StringArg)],
            synopsis="B-cell conformational epitope prediction"
        )
        reg("isis bcell conformational", desc, isis_bcell_conformational)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("min_methods", IntArg), ("min_length", IntArg)],
            synopsis="B-cell consensus prediction"
        )
        reg("isis bcell consensus", desc, isis_bcell_consensus)

        # T-cell commands
        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("allele", StringArg), ("length", IntArg)],
            synopsis="MHC Class I binding prediction"
        )
        reg("isis tcell mhc1", desc, isis_tcell_mhc1)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("allele", StringArg)],
            synopsis="MHC Class II binding prediction"
        )
        reg("isis tcell mhc2", desc, isis_tcell_mhc2)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            synopsis="Proteasomal cleavage prediction"
        )
        reg("isis tcell proteasome", desc, isis_tcell_proteasome)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            synopsis="TAP transport prediction"
        )
        reg("isis tcell tap", desc, isis_tcell_tap)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("allele", StringArg)],
            synopsis="T-cell consensus prediction"
        )
        reg("isis tcell consensus", desc, isis_tcell_consensus)

        # Innate immunity commands
        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("glyco_type", StringArg)],
            synopsis="Glycosylation site prediction"
        )
        reg("isis innate glyco", desc, isis_innate_glyco)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            synopsis="Signal peptide prediction"
        )
        reg("isis innate signal", desc, isis_innate_signal)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            synopsis="TLR ligand motif prediction"
        )
        reg("isis innate tlr", desc, isis_innate_tlr)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            synopsis="Innate immunity consensus"
        )
        reg("isis innate consensus", desc, isis_innate_consensus)

        # Structural commands
        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            synopsis="Calculate SASA"
        )
        reg("isis structure sasa", desc, isis_structure_sasa)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            synopsis="Calculate protrusion index"
        )
        reg("isis structure protrusion", desc, isis_structure_protrusion)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("cutoff", FloatArg)],
            synopsis="Calculate contact numbers"
        )
        reg("isis structure contacts", desc, isis_structure_contacts)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            synopsis="Analyze B-factors"
        )
        reg("isis structure bfactor", desc, isis_structure_bfactor)

        # Utility commands
        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("format", StringArg), ("output", StringArg)],
            synopsis="Export predictions"
        )
        reg("isis export", desc, isis_export)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("method", StringArg), ("palette", StringArg)],
            synopsis="Color by prediction scores"
        )
        reg("isis color", desc, isis_color)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            synopsis="Clear ISIS attributes"
        )
        reg("isis clear", desc, isis_clear)

        desc = CmdDesc(synopsis="Report ISIS installation state and how to fix it")
        reg("isis doctor", desc, isis_doctor)

        # Residue spec (not structure spec) so `isis plot #1/A` honours the chain.
        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("method", StringArg), ("window", IntArg),
                     ("minMethods", IntArg), ("outdir", StringArg),
                     ("prefix", StringArg)],
            synopsis="Plot prediction profiles for a structure's own sequence"
        )
        reg("isis plot", desc, isis_plot)

        desc = CmdDesc(synopsis="List prediction methods")
        reg("isis list", desc, isis_list)

        # Legacy commands (backward compatibility)
        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("method", StringArg), ("window", IntArg), ("threshold", FloatArg)],
            synopsis="Predict B-cell epitopes (legacy)"
        )
        reg("isis predict", desc, isis_predict)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("min_methods", IntArg), ("min_length", IntArg)],
            synopsis="Consensus prediction (legacy)"
        )
        reg("isis consensus", desc, isis_consensus)

        desc = CmdDesc(
            required=[("sel", ResiduesArg)],
            keyword=[("method", StringArg), ("color", StringArg)],
            synopsis="Highlight epitope regions"
        )
        reg("isis epitopes", desc, isis_epitopes)

        session.logger.info("ISIS: All commands registered successfully")

    except Exception as e:
        session.logger.error(f"ISIS: Failed to register commands: {e}")
        import traceback
        session.logger.error(traceback.format_exc())
