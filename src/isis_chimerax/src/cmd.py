"""
ChimeraX commands for ISIS epitope prediction.

Command structure:
    isis bcell linear #1 [method <name>]
    isis bcell conformational #1 [method <name>]
    isis bcell consensus #1

    isis tcell mhc1 #1 [allele <name>]
    isis tcell mhc2 #1 [allele <name>]
    isis tcell proteasome #1
    isis tcell tap #1
    isis tcell consensus #1

    isis innate glyco #1
    isis innate signal #1
    isis innate tlr #1
    isis innate consensus #1

    isis structure sasa #1
    isis structure protrusion #1
    isis structure contacts #1
    isis structure bfactor #1

    isis color #1 [method <name>] [palette <colors>]
    isis export #1 [format csv|json]
    isis clear #1
    isis list

Legacy commands (still work):
    isis predict #1 [method <name>]
    isis consensus #1
"""

from chimerax.core.commands import CmdDesc, register, StringArg, IntArg, FloatArg, BoolArg
from chimerax.atomic import AtomicStructuresArg

import numpy as np

# Import ISIS core modules
try:
    from isis import predict, predict_all, available_methods, METHOD_INFO
    ISIS_AVAILABLE = True
except ImportError:
    ISIS_AVAILABLE = False

# Import new method modules
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
except ImportError:
    METHODS_AVAILABLE = False

ATTR_PREFIX = "isis_"


# =============================================================================
# Helper Functions
# =============================================================================

def _get_chains(structure):
    """
    Extract chains with sequences from structure.

    Uses chain.residues (not chain.existing_residues): the former is
    index-parallel with chain.characters, with None standing in for any
    position not present in the model (common in cryo-EM/crystal structures
    with disordered/unresolved loops). chain.existing_residues is compacted
    and drops those gaps, which silently misaligns sequence position i with
    the wrong residue for any chain that isn't fully modeled.
    """
    chains = []
    for chain in structure.chains:
        seq = chain.characters
        if seq and len(seq) >= 6:
            chains.append((chain.chain_id, seq, chain.residues))
    return chains


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


def _auto_color(session, structure, attr_name, palette="white:yellow:red"):
    """Auto-color structure by attribute."""
    from chimerax.core.commands import run
    spec = structure.atomspec
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

def isis_bcell_linear(session, structures, method="emini", window=None, threshold=None):
    """Run B-cell linear epitope prediction."""
    if not ISIS_AVAILABLE:
        session.logger.error("ISIS library not installed")
        return

    from chimerax.atomic import Residue

    for structure in structures:
        session.logger.info(f"B-cell linear prediction ({method}) on {structure.name}...")
        attr_name = f"{ATTR_PREFIX}bcell_{method.replace('-', '_')}"
        _register_attr(session, attr_name)

        for chain_id, seq, residues in _get_chains(structure):
            try:
                pred = predict(seq, method=method, window_size=window, threshold=threshold)
                scores = np.array([pred.score_at(i+1) or 0 for i in range(len(seq))])
                _store_scores(residues, scores, attr_name)

                session.logger.info(f"  Chain {chain_id}: {len(pred.epitopes)} epitopes")
                for ep in pred.epitopes:
                    session.logger.info(f"    {ep.start}-{ep.end}: {ep.sequence}")
            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name)
        session.logger.info(f"Scores stored as: {attr_name}")


def isis_bcell_conformational(session, structures, method="discotope"):
    """Run B-cell conformational epitope prediction (requires 3D structure)."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    for structure in structures:
        session.logger.info(f"B-cell conformational prediction ({method}) on {structure.name}...")
        attr_name = f"{ATTR_PREFIX}bcell_conf_{method}"
        _register_attr(session, attr_name)

        for chain_id, seq, residues in _get_chains(structure):
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

        _auto_color(session, structure, attr_name, palette="white:cyan:blue")
        session.logger.info(f"Scores stored as: {attr_name}")


def isis_bcell_consensus(session, structures, min_methods=2, min_length=6):
    """Consensus across all B-cell methods (linear + conformational)."""
    if not ISIS_AVAILABLE:
        session.logger.error("ISIS library not installed")
        return

    linear_methods = available_methods()
    conf_methods = ["discotope", "ellipro"] if METHODS_AVAILABLE else []

    attr_name = f"{ATTR_PREFIX}bcell_consensus"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"B-cell consensus on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
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

        _auto_color(session, structure, attr_name, palette="white:yellow:orange:red")


# =============================================================================
# T-cell Epitopes
# =============================================================================

def isis_tcell_mhc1(session, structures, allele="HLA-A*02:01", length=9):
    """Predict MHC Class I binding."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = TcellPredictor()
    attr_name = f"{ATTR_PREFIX}tcell_mhc1"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"MHC-I prediction ({allele}) on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
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

        _auto_color(session, structure, attr_name, palette="white:green:darkgreen")


def isis_tcell_mhc2(session, structures, allele="HLA-DRB1*01:01"):
    """Predict MHC Class II binding."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = TcellPredictor()
    attr_name = f"{ATTR_PREFIX}tcell_mhc2"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"MHC-II prediction ({allele}) on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
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

        _auto_color(session, structure, attr_name, palette="white:purple:darkviolet")


def isis_tcell_proteasome(session, structures):
    """Predict proteasomal cleavage sites."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = TcellPredictor()
    attr_name = f"{ATTR_PREFIX}tcell_cleavage"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"Proteasomal cleavage prediction on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
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

        _auto_color(session, structure, attr_name, palette="white:orange:red")


def isis_tcell_tap(session, structures):
    """Predict TAP transport efficiency."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = TcellPredictor()
    attr_name = f"{ATTR_PREFIX}tcell_tap"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"TAP transport prediction on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
            try:
                result = predictor.predict_tap(seq)

                scores = _dense_from_tcell_result(result, len(seq))
                _store_scores(residues, scores, attr_name)
                session.logger.info(f"  Chain {chain_id}: TAP scores computed")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="white:lightblue:blue")


def isis_tcell_consensus(session, structures, allele="HLA-A*02:01"):
    """Combined T-cell epitope prediction pipeline."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = TcellPredictor()
    attr_name = f"{ATTR_PREFIX}tcell_consensus"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"T-cell consensus prediction on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
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

        _auto_color(session, structure, attr_name, palette="white:lime:green")


# =============================================================================
# Innate Immunity
# =============================================================================

def isis_innate_glyco(session, structures, glyco_type="n"):
    """Predict glycosylation sites."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = InnatePredictor()
    attr_name = f"{ATTR_PREFIX}innate_{glyco_type}glyco"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"{glyco_type.upper()}-glycosylation prediction on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
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

        _auto_color(session, structure, attr_name, palette="white:pink:magenta")


def isis_innate_signal(session, structures):
    """Predict signal peptide."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = InnatePredictor()
    attr_name = f"{ATTR_PREFIX}innate_signal"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"Signal peptide prediction on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
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

        _auto_color(session, structure, attr_name, palette="white:yellow:orange")


def isis_innate_tlr(session, structures):
    """Predict TLR ligand motifs."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = InnatePredictor()
    attr_name = f"{ATTR_PREFIX}innate_tlr"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"TLR motif prediction on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
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

        _auto_color(session, structure, attr_name, palette="white:coral:red")


def isis_innate_consensus(session, structures):
    """Consensus across all innate immunity predictions."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    predictor = InnatePredictor()
    attr_name = f"{ATTR_PREFIX}innate_consensus"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"Innate immunity consensus on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
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

        _auto_color(session, structure, attr_name, palette="white:salmon:darkred")


# =============================================================================
# Structural Analysis
# =============================================================================

def isis_structure_sasa(session, structures):
    """Calculate and display solvent accessible surface area."""
    attr_name = f"{ATTR_PREFIX}structure_sasa"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"SASA calculation on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
            try:
                sasa = _get_sasa_values(session, structure, residues)
                _store_scores(residues, sasa, attr_name)

                avg_sasa = np.mean(sasa)
                session.logger.info(f"  Chain {chain_id}: avg SASA = {avg_sasa:.1f} Å²")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="blue:white:red")


def isis_structure_protrusion(session, structures):
    """Calculate protrusion index."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    attr_name = f"{ATTR_PREFIX}structure_protrusion"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"Protrusion index on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
            try:
                coords = _get_ca_coords(residues)
                protrusion = calculate_protrusion_index(coords)
                _store_scores(residues, protrusion, attr_name)

                session.logger.info(f"  Chain {chain_id}: protrusion index computed")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="white:cyan:blue")


def isis_structure_contacts(session, structures, cutoff=8.0):
    """Calculate contact numbers."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    attr_name = f"{ATTR_PREFIX}structure_contacts"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"Contact number calculation on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
            try:
                coords = _get_ca_coords(residues)
                contacts = calculate_contact_number(coords, cutoff=cutoff)
                _store_scores(residues, contacts, attr_name)

                avg_contacts = np.mean(contacts)
                session.logger.info(f"  Chain {chain_id}: avg contacts = {avg_contacts:.1f}")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="red:white:blue")


def isis_structure_bfactor(session, structures):
    """Display normalized B-factors."""
    if not METHODS_AVAILABLE:
        session.logger.error("ISIS methods module not available")
        return

    attr_name = f"{ATTR_PREFIX}structure_bfactor"
    _register_attr(session, attr_name)

    for structure in structures:
        session.logger.info(f"B-factor analysis on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
            try:
                bfactors = _get_bfactors(residues)
                normalized = normalize_bfactors(bfactors)
                _store_scores(residues, normalized, attr_name)

                session.logger.info(f"  Chain {chain_id}: B-factors normalized")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        _auto_color(session, structure, attr_name, palette="blue:white:red")


# =============================================================================
# Export and Utility Commands
# =============================================================================

def isis_export(session, structures, format="csv", output=None):
    """Export all ISIS predictions to file."""
    import json
    import csv
    from io import StringIO

    for structure in structures:
        data = {"structure": structure.name, "chains": []}

        for chain_id, seq, residues in _get_chains(structure):
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


def isis_color(session, structures, method=None, palette="white:yellow:red"):
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

    for structure in structures:
        found_attr = None
        for attr in possible_attrs:
            if any(hasattr(res, attr) for res in structure.residues):
                found_attr = attr
                break

        if not found_attr:
            session.logger.warning(f"{structure.name}: No predictions found for '{method}'")
            continue

        _auto_color(session, structure, found_attr, palette)
        session.logger.info(f"Colored {structure.name} by {found_attr}")


def isis_clear(session, structures):
    """Remove all ISIS prediction attributes."""
    for structure in structures:
        cleared = set()
        for res in structure.residues:
            for attr in list(vars(res).keys()):
                if attr.startswith(ATTR_PREFIX):
                    delattr(res, attr)
                    cleared.add(attr)
        if cleared:
            session.logger.info(f"{structure.name}: Cleared {len(cleared)} attributes")
        else:
            session.logger.info(f"{structure.name}: No ISIS attributes found")


def isis_list(session):
    """List all available prediction methods."""
    session.logger.info("ISIS Prediction Methods:")
    session.logger.info("")

    session.logger.info("B-cell Linear (sequence-based):")
    if ISIS_AVAILABLE:
        for method in available_methods():
            info = METHOD_INFO.get(method, {})
            session.logger.info(f"  {method}: {info.get('name', method)}")
    else:
        session.logger.info("  (ISIS core not installed)")

    session.logger.info("")
    session.logger.info("B-cell Conformational (structure-based):")
    if METHODS_AVAILABLE:
        session.logger.info("  discotope: DiscoTope-style (SASA + propensity + contacts)")
        session.logger.info("  ellipro: ElliPro-style (protrusion index)")
        session.logger.info("  seppa: SEPPA-style (surface patches)")
    else:
        session.logger.info("  (methods module not installed)")

    session.logger.info("")
    session.logger.info("T-cell (MHC binding):")
    if METHODS_AVAILABLE:
        session.logger.info("  mhc1: MHC Class I binding (HLA-A*02:01, HLA-A*01:01, HLA-B*07:02)")
        session.logger.info("  mhc2: MHC Class II binding (HLA-DRB1*01:01, HLA-DRB1*04:01)")
        session.logger.info("  proteasome: Proteasomal cleavage sites")
        session.logger.info("  tap: TAP transport efficiency")
    else:
        session.logger.info("  (methods module not installed)")

    session.logger.info("")
    session.logger.info("Innate Immunity:")
    if METHODS_AVAILABLE:
        session.logger.info("  glyco: N- and O-glycosylation sites")
        session.logger.info("  signal: Signal peptide detection")
        session.logger.info("  tlr: TLR ligand motifs")
    else:
        session.logger.info("  (methods module not installed)")

    session.logger.info("")
    session.logger.info("Structural Analysis:")
    session.logger.info("  sasa: Solvent accessible surface area")
    session.logger.info("  protrusion: Protrusion index")
    session.logger.info("  contacts: Residue contact numbers")
    session.logger.info("  bfactor: Normalized B-factors")


# =============================================================================
# Legacy Commands (backward compatibility)
# =============================================================================

def isis_predict(session, structures, method="emini", window=None, threshold=None):
    """Legacy predict command - redirects to isis_bcell_linear."""
    isis_bcell_linear(session, structures, method=method, window=window, threshold=threshold)


def isis_consensus(session, structures, min_methods=2, min_length=6):
    """Legacy consensus command - redirects to isis_bcell_consensus."""
    isis_bcell_consensus(session, structures, min_methods=min_methods, min_length=min_length)


def isis_epitopes(session, structures, method="emini", color="red"):
    """Color only predicted epitope regions."""
    if not ISIS_AVAILABLE:
        session.logger.error("ISIS library not installed")
        return

    from chimerax.core.commands import run

    for structure in structures:
        for chain_id, seq, residues in _get_chains(structure):
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
            required=[("structures", AtomicStructuresArg)],
            keyword=[("method", StringArg), ("window", IntArg), ("threshold", FloatArg)],
            synopsis="B-cell linear epitope prediction"
        )
        reg("isis bcell linear", desc, isis_bcell_linear)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[("method", StringArg)],
            synopsis="B-cell conformational epitope prediction"
        )
        reg("isis bcell conformational", desc, isis_bcell_conformational)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[("min_methods", IntArg), ("min_length", IntArg)],
            synopsis="B-cell consensus prediction"
        )
        reg("isis bcell consensus", desc, isis_bcell_consensus)

        # T-cell commands
        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[("allele", StringArg), ("length", IntArg)],
            synopsis="MHC Class I binding prediction"
        )
        reg("isis tcell mhc1", desc, isis_tcell_mhc1)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[("allele", StringArg)],
            synopsis="MHC Class II binding prediction"
        )
        reg("isis tcell mhc2", desc, isis_tcell_mhc2)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            synopsis="Proteasomal cleavage prediction"
        )
        reg("isis tcell proteasome", desc, isis_tcell_proteasome)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            synopsis="TAP transport prediction"
        )
        reg("isis tcell tap", desc, isis_tcell_tap)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[("allele", StringArg)],
            synopsis="T-cell consensus prediction"
        )
        reg("isis tcell consensus", desc, isis_tcell_consensus)

        # Innate immunity commands
        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[("glyco_type", StringArg)],
            synopsis="Glycosylation site prediction"
        )
        reg("isis innate glyco", desc, isis_innate_glyco)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            synopsis="Signal peptide prediction"
        )
        reg("isis innate signal", desc, isis_innate_signal)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            synopsis="TLR ligand motif prediction"
        )
        reg("isis innate tlr", desc, isis_innate_tlr)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            synopsis="Innate immunity consensus"
        )
        reg("isis innate consensus", desc, isis_innate_consensus)

        # Structural commands
        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            synopsis="Calculate SASA"
        )
        reg("isis structure sasa", desc, isis_structure_sasa)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            synopsis="Calculate protrusion index"
        )
        reg("isis structure protrusion", desc, isis_structure_protrusion)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[("cutoff", FloatArg)],
            synopsis="Calculate contact numbers"
        )
        reg("isis structure contacts", desc, isis_structure_contacts)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            synopsis="Analyze B-factors"
        )
        reg("isis structure bfactor", desc, isis_structure_bfactor)

        # Utility commands
        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[("format", StringArg), ("output", StringArg)],
            synopsis="Export predictions"
        )
        reg("isis export", desc, isis_export)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[("method", StringArg), ("palette", StringArg)],
            synopsis="Color by prediction scores"
        )
        reg("isis color", desc, isis_color)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            synopsis="Clear ISIS attributes"
        )
        reg("isis clear", desc, isis_clear)

        desc = CmdDesc(synopsis="List prediction methods")
        reg("isis list", desc, isis_list)

        # Legacy commands (backward compatibility)
        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[("method", StringArg), ("window", IntArg), ("threshold", FloatArg)],
            synopsis="Predict B-cell epitopes (legacy)"
        )
        reg("isis predict", desc, isis_predict)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[("min_methods", IntArg), ("min_length", IntArg)],
            synopsis="Consensus prediction (legacy)"
        )
        reg("isis consensus", desc, isis_consensus)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[("method", StringArg), ("color", StringArg)],
            synopsis="Highlight epitope regions"
        )
        reg("isis epitopes", desc, isis_epitopes)

        session.logger.info("ISIS: All commands registered successfully")

    except Exception as e:
        session.logger.error(f"ISIS: Failed to register commands: {e}")
        import traceback
        session.logger.error(traceback.format_exc())
