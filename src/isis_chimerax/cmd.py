"""
ChimeraX commands for ISIS epitope prediction.

Commands:
    isis predict #1 [method emini] [window 7]
    isis color #1 [method emini] [palette white:yellow:red]
    isis clear #1
    isis list
"""

from chimerax.core.commands import (
    CmdDesc, register, StringArg, IntArg, FloatArg, ModelsArg, EnumOf
)
from chimerax.atomic import AtomicStructure

# Import ISIS core (handle case where it's not installed)
try:
    from isis import predict, predict_all, available_methods, METHOD_INFO
    ISIS_AVAILABLE = True
except ImportError:
    ISIS_AVAILABLE = False


ATTR_PREFIX = "isis_"


def _get_sequence(structure):
    """Extract sequence from a ChimeraX atomic structure."""
    sequences = []
    for chain in structure.chains:
        residues = chain.existing_residues
        seq = "".join(r.one_letter_code or "X" for r in residues)
        sequences.append((chain.chain_id, seq, residues))
    return sequences


def _store_predictions(session, structure, predictions, method):
    """Store prediction results as residue attributes."""
    attr_name = f"{ATTR_PREFIX}{method.replace('-', '_')}"

    for chain_id, seq, residues in _get_sequence(structure):
        if method not in predictions.get(chain_id, {}):
            continue

        pred = predictions[chain_id][method]

        for i, res in enumerate(residues):
            pos = i + 1
            score = pred.score_at(pos)
            if score is not None:
                setattr(res, attr_name, score)

    return attr_name


def isis_predict(session, structures, method="emini", window=None, threshold=None):
    """
    Run B-cell epitope prediction on structure sequences.

    Results are stored as residue attributes for visualization.
    """
    if not ISIS_AVAILABLE:
        session.logger.error("ISIS library not installed. Run: pip install isis-epitope")
        return

    if not structures:
        session.logger.error("No structures specified")
        return

    for structure in structures:
        session.logger.info(f"Predicting epitopes for {structure.name}...")

        chain_predictions = {}
        for chain_id, seq, residues in _get_sequence(structure):
            if len(seq) < 6:
                session.logger.warning(f"Chain {chain_id}: sequence too short ({len(seq)} aa)")
                continue

            try:
                pred = predict(
                    seq,
                    method=method,
                    window_size=window,
                    threshold=threshold,
                    sequence_name=f"{structure.name}_{chain_id}"
                )
                chain_predictions[chain_id] = {method: pred}

                n_epitopes = len(pred.epitopes)
                session.logger.info(
                    f"  Chain {chain_id}: {len(seq)} aa, {n_epitopes} epitopes predicted"
                )
                for ep in pred.epitopes:
                    session.logger.info(f"    {ep.start}-{ep.end}: {ep.sequence}")

            except Exception as e:
                session.logger.error(f"  Chain {chain_id}: {e}")

        attr_name = _store_predictions(session, structure, chain_predictions, method)
        session.logger.info(f"Scores stored as attribute: {attr_name}")


def isis_color(session, structures, method="emini", palette="white:yellow:red"):
    """
    Color structure by ISIS prediction scores.

    Uses rangecol to map scores to colors.
    """
    from chimerax.core.commands import run

    attr_name = f"{ATTR_PREFIX}{method.replace('-', '_')}"
    colors = palette.split(":")

    for structure in structures:
        spec = structure.atomspec

        # Check if attribute exists
        has_attr = False
        for res in structure.residues:
            if hasattr(res, attr_name):
                has_attr = True
                break

        if not has_attr:
            session.logger.warning(
                f"{structure.name}: No {method} predictions found. Run 'isis predict' first."
            )
            continue

        # Build rangecol command
        if len(colors) == 2:
            cmd = f"color byattribute {attr_name} {spec} palette {colors[0]}:{colors[1]}"
        elif len(colors) == 3:
            cmd = f"color byattribute {attr_name} {spec} palette {colors[0]}:{colors[1]}:{colors[2]}"
        else:
            cmd = f"color byattribute {attr_name} {spec} palette white:yellow:red"

        run(session, cmd)
        session.logger.info(f"Colored {structure.name} by {method}")


def isis_color_epitopes(session, structures, method="emini", color="red"):
    """
    Color only the predicted epitope regions.
    """
    if not ISIS_AVAILABLE:
        session.logger.error("ISIS library not installed")
        return

    from chimerax.core.commands import run

    for structure in structures:
        for chain_id, seq, residues in _get_sequence(structure):
            try:
                pred = predict(seq, method=method, sequence_name=chain_id)

                for epitope in pred.epitopes:
                    # Build residue selection
                    res_range = f":{epitope.start}-{epitope.end}"
                    if chain_id:
                        res_range = f"/{chain_id}{res_range}"

                    spec = f"{structure.atomspec}{res_range}"
                    run(session, f"color {spec} {color}")

            except Exception as e:
                session.logger.error(f"Chain {chain_id}: {e}")


def isis_clear(session, structures):
    """
    Remove ISIS prediction attributes from structures.
    """
    for structure in structures:
        cleared = []
        for res in structure.residues:
            for attr in list(vars(res).keys()):
                if attr.startswith(ATTR_PREFIX):
                    delattr(res, attr)
                    if attr not in cleared:
                        cleared.append(attr)

        if cleared:
            session.logger.info(f"{structure.name}: Cleared {', '.join(cleared)}")
        else:
            session.logger.info(f"{structure.name}: No ISIS attributes found")


def isis_list(session):
    """
    List available ISIS prediction methods.
    """
    if not ISIS_AVAILABLE:
        session.logger.error("ISIS library not installed")
        return

    session.logger.info("Available ISIS prediction methods:\n")
    for method in available_methods():
        info = METHOD_INFO[method]
        session.logger.info(f"  {method}")
        session.logger.info(f"    {info['name']}")
        session.logger.info(f"    Window: {info['default_window']}")
        session.logger.info(f"    {info['description']}\n")


def register_commands(command_info, logger):
    """Register ISIS commands with ChimeraX."""

    method_choices = ["emini", "parker", "chou-fasman", "kolaskar-tongaonkar", "karplus-schulz"]

    # isis predict
    predict_desc = CmdDesc(
        required=[("structures", ModelsArg)],
        keyword=[
            ("method", EnumOf(method_choices)),
            ("window", IntArg),
            ("threshold", FloatArg),
        ],
        synopsis="Predict B-cell epitopes on structure sequences"
    )
    register("isis predict", predict_desc, isis_predict, logger=logger)

    # isis color
    color_desc = CmdDesc(
        required=[("structures", ModelsArg)],
        keyword=[
            ("method", EnumOf(method_choices)),
            ("palette", StringArg),
        ],
        synopsis="Color structure by epitope prediction scores"
    )
    register("isis color", color_desc, isis_color, logger=logger)

    # isis epitopes (color only epitopes)
    epitopes_desc = CmdDesc(
        required=[("structures", ModelsArg)],
        keyword=[
            ("method", EnumOf(method_choices)),
            ("color", StringArg),
        ],
        synopsis="Highlight predicted epitope regions"
    )
    register("isis epitopes", epitopes_desc, isis_color_epitopes, logger=logger)

    # isis clear
    clear_desc = CmdDesc(
        required=[("structures", ModelsArg)],
        synopsis="Remove ISIS prediction attributes"
    )
    register("isis clear", clear_desc, isis_clear, logger=logger)

    # isis list
    list_desc = CmdDesc(synopsis="List available prediction methods")
    register("isis list", list_desc, isis_list, logger=logger)
