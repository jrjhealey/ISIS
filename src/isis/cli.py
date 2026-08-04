#!/usr/bin/env python3
"""
Command-line interface for ISIS epitope prediction.

Usage:
    isis predict sequence.fasta
    isis predict sequence.fasta --methods emini,parker --window 7
    isis predict sequence.fasta --output results.csv
    isis list-methods
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import List, Tuple

from . import predict, predict_all, available_methods, __version__


def parse_fasta(text: str) -> List[Tuple[str, str]]:
    """Parse FASTA format, return list of (name, sequence) tuples."""
    sequences = []
    current_name = None
    current_seq = []

    for line in text.strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name is not None:
                sequences.append((current_name, "".join(current_seq)))
            current_name = line[1:].split()[0] or f"Sequence_{len(sequences) + 1}"
            current_seq = []
        else:
            current_seq.append(line.replace(" ", ""))

    if current_name is not None:
        sequences.append((current_name, "".join(current_seq)))
    elif current_seq:
        sequences.append(("Sequence_1", "".join(current_seq)))

    return sequences


def cmd_predict(args):
    """Run epitope prediction."""
    # Read input
    if args.input == "-":
        text = sys.stdin.read()
    else:
        text = Path(args.input).read_text()

    sequences = parse_fasta(text)
    if not sequences:
        print("Error: No sequences found in input", file=sys.stderr)
        sys.exit(1)

    # Parse methods
    if args.methods:
        methods = [m.strip() for m in args.methods.split(",")]
    else:
        methods = ["emini", "parker", "chou-fasman", "kolaskar-tongaonkar"]

    # Run predictions
    all_results = []
    for name, seq in sequences:
        results = predict_all(seq, methods, args.window, sequence_name=name)
        all_results.append((name, seq, results))

    # Output
    if args.format == "table":
        output_table(all_results, args.output)
    elif args.format == "csv":
        output_csv(all_results, args.output)
    elif args.format == "json":
        output_json(all_results, args.output)
    elif args.format == "epitopes":
        output_epitopes(all_results, args.output)


def output_table(results, output_path):
    """Output as human-readable table."""
    lines = []

    for name, seq, predictions in results:
        lines.append(f"\n{'=' * 60}")
        lines.append(f"Sequence: {name}")
        lines.append(f"Length: {len(seq)}")
        lines.append(f"{'=' * 60}")

        # Build score table
        first_pred = next(iter(predictions.values()))
        positions = first_pred.positions

        header = ["Pos", "AA"] + list(predictions.keys())
        lines.append("  ".join(f"{h:>12}" for h in header))
        lines.append("-" * (14 * len(header)))

        for i, pos in enumerate(positions):
            aa = seq[int(pos) - 1]
            row = [f"{int(pos):>12}", f"{aa:>12}"]
            for method, pred in predictions.items():
                score = pred.scores[i] if i < len(pred.scores) else float("nan")
                row.append(f"{score:>12.3f}")
            lines.append("  ".join(row))

        # Epitope summary
        lines.append(f"\n{'Predicted Epitopes':^60}")
        lines.append("-" * 60)
        for method, pred in predictions.items():
            lines.append(f"\n{method} (threshold={pred.threshold:.3f}):")
            if pred.epitopes:
                for ep in pred.epitopes:
                    lines.append(f"  {ep}")
            else:
                lines.append("  No epitopes above threshold")

    output = "\n".join(lines)
    if output_path:
        Path(output_path).write_text(output)
    else:
        print(output)


def output_csv(results, output_path):
    """Output as CSV with all scores."""
    out = sys.stdout if not output_path else open(output_path, "w", newline="")
    writer = csv.writer(out)

    for name, seq, predictions in results:
        methods = list(predictions.keys())
        first_pred = next(iter(predictions.values()))

        writer.writerow(["Sequence", name])
        writer.writerow(["Position", "Residue"] + methods)

        for i, pos in enumerate(first_pred.positions):
            aa = seq[int(pos) - 1]
            row = [int(pos), aa]
            for method in methods:
                pred = predictions[method]
                score = pred.scores[i] if i < len(pred.scores) else ""
                row.append(f"{score:.4f}" if isinstance(score, float) else score)
            writer.writerow(row)

        writer.writerow([])

    if output_path:
        out.close()


def output_json(results, output_path):
    """Output as JSON."""
    data = []
    for name, seq, predictions in results:
        entry = {
            "sequence_name": name,
            "sequence": seq,
            "predictions": {m: p.to_dict() for m, p in predictions.items()},
        }
        data.append(entry)

    output = json.dumps(data, indent=2)
    if output_path:
        Path(output_path).write_text(output)
    else:
        print(output)


def output_epitopes(results, output_path):
    """Output only epitope calls (compact format)."""
    lines = []
    for name, seq, predictions in results:
        lines.append(f">{name}")
        for method, pred in predictions.items():
            for ep in pred.epitopes:
                lines.append(f"{method}\t{ep.start}\t{ep.end}\t{ep.sequence}\t{ep.score:.3f}")

    output = "\n".join(lines)
    if output_path:
        Path(output_path).write_text(output)
    else:
        print(output)


def cmd_list_methods(args):
    """List available prediction methods."""
    from .scales import METHOD_INFO

    print("Available prediction methods:\n")
    for method in available_methods():
        info = METHOD_INFO[method]
        print(f"  {method}")
        print(f"    {info['name']}")
        print(f"    Window: {info['default_window']}, Threshold: {info['default_threshold'] or 'average'}")
        print(f"    {info['description']}")
        print()


def main():
    parser = argparse.ArgumentParser(
        prog="isis",
        description="ISIS - In Silico Immunogenicity Studies: B-cell epitope prediction",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # predict command
    pred_parser = subparsers.add_parser("predict", help="Run epitope prediction")
    pred_parser.add_argument("input", help="Input FASTA file (or - for stdin)")
    pred_parser.add_argument(
        "-m", "--methods",
        help="Comma-separated methods (default: emini,parker,chou-fasman,kolaskar-tongaonkar)"
    )
    pred_parser.add_argument(
        "-w", "--window",
        type=int,
        help="Window size (default: method-specific)"
    )
    pred_parser.add_argument(
        "-f", "--format",
        choices=["table", "csv", "json", "epitopes"],
        default="table",
        help="Output format (default: table)"
    )
    pred_parser.add_argument(
        "-o", "--output",
        help="Output file (default: stdout)"
    )
    pred_parser.set_defaults(func=cmd_predict)

    # list-methods command
    list_parser = subparsers.add_parser("list-methods", help="List available methods")
    list_parser.set_defaults(func=cmd_list_methods)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
