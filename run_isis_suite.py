#!/usr/bin/env python3
"""
A wrapper script to run the various tools within the
in silico immunogenicity IEDB suite
"""
from __future__ import print_function
import argparse
import logging
import sys


NO_COLOR = "\33[m"
RED, GREEN, ORANGE, BLUE, PURPLE, LBLUE, GREY = map("\33[%dm".__mod__, range(31, 38))
logging.basicConfig(format="[%(asctime)s] %(levelname)-8s->  %(message)s",
                    level=logging.NOTSET, datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger(__name__)

def add_color(logger_method, color):
    def wrapper(message, *args, **kwargs):
        return logger_method(color+message+NO_COLOR, *args, **kwargs)
    return wrapper


def get_args():
    """Parse command line arguments"""
    desc = "A wrapper for the IEDBtools suite."
    epi = ("This script pulls in several objects/other scripts from the"
           "IEDB suite and collates their outputs to provide information")

    try:
        parser = argparse.ArgumentParser(description=desc,
            epilog=epi,
            prog=__file__)
        parser.add_argument(
            "-v",
            "--verbose",
            default="1",
            choices=["0", "1", "2"],
            metavar="Set verbosity level for logging.",
            help="Change the verbosity of output (info/debug/errors)."
        )
        parser.add_argument(
            "-i",
            "--infile",
            action="store",
            help="The sequence to analyse"
        )
        #parser.add_argument('outfile', action='store', default=None,
        #                    help='Output file of primers in the chosen format.')
        parser.add_argument(
            "-w",
            "--window_size",
            action="store",
            type=int,
            default=6,
            help="Window size must be defined for aggregating data, but "
                 "is actually defined on a per-method basis if set as False"
        )

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")
        sys.exit(1)

    return parser.parse_args()


def main():
    """Main to dispatch runners and begin the data collation"""
    args = get_args()

    assert int(args.verbose) < 3, ("verbose supports maximum 3 "
                                   "levels at present [0, 1, 2], "
                                   " got: {}".format(args.verbose))
    levels_dict = {"0": logging.WARNING, "1": logging.INFO, "2": logging.DEBUG}
    logging.getLogger().setLevel(levels_dict[args.verbose])
    logger.info("Launching {}...".format(__file__))

    # Run the BCellStandalone tools
    from bcell_standalone.BCellRunner import BCellRunner
    import pandas as pd
    pd.set_option('expand_frame_repr', False)

    BCR = BCellRunner(args.infile, args.window_size)

    # Since the first 5 columns will be the same, select these from the first result
    # The datastructure from this tool is also fucking insane so this code is gonna
    # be weird... sorry....
    columns = [col for col in BCR.results['Emini'].values()[0]['prediction_result'][0][0:5] if col != 'Residue']
    collated_df = pd.DataFrame(columns=columns)

    for i, value in enumerate(BCR.results['Emini'].values()[0]['prediction_result'][0][0:5]):
        if i == 1:
            # Skip the residue column since its meaningless
            continue
        else:
            collated_df[value] = [x[i] for x in BCR.results['Emini'].values()[0]['prediction_result'][1:]]

    # Now for each method, append new columns laterally with the scores

    # NOTE!!!
    # The 'residue' column is incorrect w.r.t. the Karplus method, as it considers
    # a slightly different position. For a given window size however, the window
    # is the same and correct.
    # NaNs introduced by different length arrays are made 0 (or maybe threshold?)
    for method in sorted(BCR.results.keys()):
        try:
            collated_df[method] = pd.DataFrame([x[5] for x in BCR.results[method].values()[0]['prediction_result'][1:]])
        except ValueError as e:
            logger.error(e + '\n' + "Different length arrays were returned "
                                    "from the methods, (usually Karplus). "
                                    "Try again with a fixed window size.")


    print(collated_df)









if __name__ == "__main__":
    main()
