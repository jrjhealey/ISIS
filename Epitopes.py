#!/usr/bin/env python2.7
"""
A wrapper script to run the various tools within the
in silico immunogenicity IEDB suite
"""
from __future__ import print_function
import argparse
import logging
import sys

from bcell_standalone.BCellRunner import BCellRunner

logging.basicConfig(format="[%(asctime)s] %(levelname)-8s->  %(message)s",
                    level=logging.NOTSET, datefmt='%d/%m/%Y %H:%M:%S %p')
logger = logging.getLogger(__name__)

def get_args():
    """Parse command line arguments"""
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-v",
            "--verbose",
            default="1",
            choices=["0", "1", "2"],
            help="Change the verbosity of output (info/debug/errors)."
        )
        parser.add_argument(
            "-i",
            "--infile",
            action="store",
            help="The sequence to analyse"
        )
        parser.add_argument(
            "-w",
            "--window_size",
            action="store",
            type=int,
            default=6,
            help="Window size for peptides"
        )
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")
        sys.exit(1)

    return parser.parse_args()


def main():
    args = get_args()
    assert int(args.verbose) < 3, ("verbose supports maximum 3 "
                                   "levels at present [0, 1, 2], "
                                   " got: {}".format(args.verbose))
    levels_dict = {"0": logging.WARNING, "1": logging.INFO, "2": logging.DEBUG}
    logging.getLogger().setLevel(levels_dict[args.verbose])
    logger.info("Launching {}...".format(__file__))
    logger.info("Calling BCellRunner...")

    BCR = BCellRunner(args.infile, args.window_size)
    epitopes = {method:BCR.results[method].values()[0]['epitopes'] for method in BCR.results}
    print(epitopes)


if __name__ == "__main__":
    main()
