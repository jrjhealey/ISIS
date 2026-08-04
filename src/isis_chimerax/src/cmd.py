"""
ChimeraX commands for ISIS epitope prediction.
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
METHOD_CHOICES = ["emini", "parker", "chou-fasman", "kolaskar-tongaonkar", "karplus-schulz"]


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

    for structure in structures:
        session.logger.info(f"Predicting epitopes for {structure.name}...")
        attr_name = f"{ATTR_PREFIX}{method.replace('-', '_')}"

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
                for i, res in enumerate(residues):
                    pos = i + 1
                    score = pred.score_at(pos)
                    if score is not None:
                        setattr(res, attr_name, score)

                n_epitopes = len(pred.epitopes)
                session.logger.info(f"  Chain {chain_id}: {len(seq)} aa, {n_epitopes} epitopes")
                for ep in pred.epitopes:
                    session.logger.info(f"    {ep.start}-{ep.end}: {ep.sequence}")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

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

    session.logger.info("Available ISIS prediction methods:\n")
    for method in available_methods():
        info = METHOD_INFO[method]
        session.logger.info(f"  {method}: {info['name']}")


def register_command(command_name, logger):
    """Register a single ISIS command."""

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
