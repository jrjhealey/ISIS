"""
ChimeraX commands for ISIS epitope prediction.

If commands don't register automatically, run in ChimeraX Python shell:
    from chimerax.isis.cmd import register_all_commands
    register_all_commands(session)
"""

from chimerax.core.commands import CmdDesc, register, StringArg, IntArg, FloatArg
from chimerax.atomic import AtomicStructuresArg

# Import ISIS core
try:
    from isis import predict, predict_all, available_methods, METHOD_INFO
    ISIS_AVAILABLE = True
except ImportError:
    ISIS_AVAILABLE = False

ATTR_PREFIX = "isis_"


def _get_chains(structure):
    """Extract chains with sequences from structure."""
    chains = []
    for chain in structure.chains:
        seq = chain.characters
        if seq and len(seq) >= 6:
            chains.append((chain.chain_id, seq, chain.existing_residues))
    return chains


def isis_predict(session, structures, method="emini", window=None, threshold=None):
    """Run B-cell epitope prediction on structure sequences."""
    if not ISIS_AVAILABLE:
        session.logger.error("ISIS library not installed. Run: pip install isis-epitope")
        return

    if not structures:
        session.logger.error("No structures specified")
        return

    from chimerax.atomic import Residue

    for structure in structures:
        session.logger.info(f"Predicting epitopes for {structure.name}...")
        attr_name = f"{ATTR_PREFIX}{method.replace('-', '_')}"

        # Register the attribute with ChimeraX
        if not hasattr(Residue, attr_name):
            Residue.register_attr(session, attr_name, "ISIS", attr_type=float)

        for chain_id, seq, residues in _get_chains(structure):
            try:
                pred = predict(
                    seq,
                    method=method,
                    window_size=window,
                    threshold=threshold,
                    sequence_name=f"{structure.name}_{chain_id}"
                )

                # Store scores as residue attributes
                scores_set = 0
                for i, res in enumerate(residues):
                    pos = i + 1
                    score = pred.score_at(pos)
                    if score is not None:
                        setattr(res, attr_name, float(score))
                        scores_set += 1

                n_epitopes = len(pred.epitopes)
                session.logger.info(f"  Chain {chain_id}: {len(seq)} aa, {scores_set} scores, {n_epitopes} epitopes")
                for ep in pred.epitopes:
                    session.logger.info(f"    {ep.start}-{ep.end}: {ep.sequence}")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")
                import traceback
                session.logger.error(traceback.format_exc())

        session.logger.info(f"Scores stored as attribute: {attr_name}")


def isis_color(session, structures, method="emini", palette="white:yellow:red"):
    """Color structure by ISIS prediction scores."""
    from chimerax.core.commands import run

    attr_name = f"{ATTR_PREFIX}{method.replace('-', '_')}"

    for structure in structures:
        has_attr = any(hasattr(res, attr_name) for res in structure.residues)
        if not has_attr:
            session.logger.warning(
                f"{structure.name}: No {method} predictions. Run 'isis predict' first."
            )
            continue

        spec = structure.atomspec
        cmd = f"color byattribute {attr_name} {spec} palette {palette}"
        run(session, cmd)
        session.logger.info(f"Colored {structure.name} by {method}")


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


def isis_consensus(session, structures, min_methods=2, min_length=6):
    """
    Run all prediction methods and create consensus score.

    The consensus score is the number of methods (0-5) that predict
    each position as part of an epitope. Only considers epitopes of
    min_length or longer. If no consensus epitopes are found, iteratively
    loosens thresholds until consensus regions are identified.
    """
    if not ISIS_AVAILABLE:
        session.logger.error("ISIS library not installed")
        return

    from chimerax.atomic import Residue
    from chimerax.core.commands import run

    methods = available_methods()
    attr_name = f"{ATTR_PREFIX}consensus"

    # Register the attribute
    if not hasattr(Residue, attr_name):
        Residue.register_attr(session, attr_name, "ISIS", attr_type=float)

    # Threshold multipliers to try (progressively looser)
    threshold_multipliers = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]

    for structure in structures:
        session.logger.info(f"Running consensus analysis on {structure.name}...")

        for chain_id, seq, residues in _get_chains(structure):
            best_votes = None
            best_multiplier = 1.0
            found_consensus = False

            # Try progressively looser thresholds until we find consensus
            for multiplier in threshold_multipliers:
                votes = [0] * len(seq)

                for method in methods:
                    try:
                        # Get default threshold for method and apply multiplier
                        info = METHOD_INFO.get(method, {})
                        default_thresh = info.get('default_threshold')

                        if default_thresh is not None:
                            thresh = default_thresh * multiplier
                        else:
                            thresh = None  # Will use average

                        pred = predict(seq, method=method, threshold=thresh)

                        # Only count epitopes that meet minimum length
                        for epitope in pred.epitopes:
                            if epitope.length >= min_length:
                                for pos in range(epitope.start, epitope.end + 1):
                                    if 1 <= pos <= len(seq):
                                        votes[pos - 1] += 1
                    except Exception as e:
                        session.logger.warning(f"  {method} failed: {e}")

                # Check if we have any consensus regions
                consensus_regions = _find_consensus_regions(votes, min_methods, min_length)

                if consensus_regions:
                    best_votes = votes
                    best_multiplier = multiplier
                    found_consensus = True
                    break

                # Save first attempt even if no consensus
                if best_votes is None:
                    best_votes = votes

            if found_consensus:
                if best_multiplier < 1.0:
                    session.logger.info(f"  Chain {chain_id}: Found consensus with threshold multiplier {best_multiplier}")
            else:
                session.logger.info(f"  Chain {chain_id}: No consensus epitopes found (showing raw scores)")

            # Store consensus scores
            for i, res in enumerate(residues):
                if i < len(best_votes):
                    setattr(res, attr_name, float(best_votes[i]))

            # Report consensus regions
            consensus_regions = _find_consensus_regions(best_votes, min_methods, min_length)
            session.logger.info(f"  Chain {chain_id}: {len(consensus_regions)} consensus epitope(s)")
            for start, end, avg_score in consensus_regions:
                peptide = seq[start-1:end]
                session.logger.info(f"    {start}-{end}: {peptide} (avg {avg_score:.1f} methods)")

        # Auto-color
        spec = structure.atomspec
        run(session, f"color byattribute {attr_name} {spec} palette white:yellow:orange:red")

        session.logger.info(f"Consensus scores stored as: {attr_name}")
        session.logger.info(f"  Color: white(0) → yellow(1-2) → orange(3-4) → red(5)")


def _find_consensus_regions(votes, min_methods, min_length):
    """Find contiguous regions where min_methods or more agree."""
    regions = []
    in_region = False
    start = 0
    scores = []

    for i, vote in enumerate(votes):
        pos = i + 1
        if vote >= min_methods:
            if not in_region:
                in_region = True
                start = pos
                scores = [vote]
            else:
                scores.append(vote)
        else:
            if in_region:
                in_region = False
                length = (pos - 1) - start + 1
                if length >= min_length:
                    avg_score = sum(scores) / len(scores)
                    regions.append((start, pos - 1, avg_score))
                scores = []

    # Handle region at end
    if in_region:
        length = len(votes) - start + 1
        if length >= min_length:
            avg_score = sum(scores) / len(scores)
            regions.append((start, len(votes), avg_score))

    return regions


def isis_clear(session, structures):
    """Remove ISIS prediction attributes."""
    for structure in structures:
        cleared = set()
        for res in structure.residues:
            for attr in list(vars(res).keys()):
                if attr.startswith(ATTR_PREFIX):
                    delattr(res, attr)
                    cleared.add(attr)
        if cleared:
            session.logger.info(f"{structure.name}: Cleared {', '.join(cleared)}")
        else:
            session.logger.info(f"{structure.name}: No ISIS attributes found")


def isis_list(session):
    """List available prediction methods."""
    if not ISIS_AVAILABLE:
        session.logger.error("ISIS library not installed")
        return

    session.logger.info("Available ISIS prediction methods:")
    for method in available_methods():
        info = METHOD_INFO[method]
        session.logger.info(f"  {method}: {info['name']}")


def register_command(session, command_name):
    """Register a single ISIS command."""
    logger = session.logger

    if command_name == "isis predict":
        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[
                ("method", StringArg),
                ("window", IntArg),
                ("threshold", FloatArg),
            ],
            synopsis="Predict B-cell epitopes"
        )
        register("isis predict", desc, isis_predict, logger=logger)

    elif command_name == "isis color":
        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[
                ("method", StringArg),
                ("palette", StringArg),
            ],
            synopsis="Color by epitope scores"
        )
        register("isis color", desc, isis_color, logger=logger)

    elif command_name == "isis epitopes":
        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[
                ("method", StringArg),
                ("color", StringArg),
            ],
            synopsis="Highlight epitope regions"
        )
        register("isis epitopes", desc, isis_epitopes, logger=logger)

    elif command_name == "isis clear":
        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            synopsis="Clear ISIS attributes"
        )
        register("isis clear", desc, isis_clear, logger=logger)

    elif command_name == "isis list":
        desc = CmdDesc(synopsis="List prediction methods")
        register("isis list", desc, isis_list, logger=logger)


def register_all_commands(session):
    """Register all ISIS commands at bundle initialization."""
    from chimerax.core.commands import register as reg

    try:
        session.logger.info("ISIS: Registering commands...")

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[
                ("method", StringArg),
                ("window", IntArg),
                ("threshold", FloatArg),
            ],
            synopsis="Predict B-cell epitopes"
        )
        reg("isis predict", desc, isis_predict)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[
                ("method", StringArg),
                ("palette", StringArg),
            ],
            synopsis="Color by epitope scores"
        )
        reg("isis color", desc, isis_color)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[
                ("method", StringArg),
                ("color", StringArg),
            ],
            synopsis="Highlight epitope regions"
        )
        reg("isis epitopes", desc, isis_epitopes)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            synopsis="Clear ISIS attributes"
        )
        reg("isis clear", desc, isis_clear)

        desc = CmdDesc(synopsis="List prediction methods")
        reg("isis list", desc, isis_list)

        desc = CmdDesc(
            required=[("structures", AtomicStructuresArg)],
            keyword=[
                ("min_methods", IntArg),
                ("min_length", IntArg),
            ],
            synopsis="Consensus epitope prediction (all methods)"
        )
        reg("isis consensus", desc, isis_consensus)

        session.logger.info("ISIS: Commands registered successfully")

    except Exception as e:
        session.logger.error(f"ISIS: Failed to register commands: {e}")
        import traceback
        session.logger.error(traceback.format_exc())
