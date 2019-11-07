"""
Created on: 03/18/2014
Updated on: 05/19/2017

@author: Dorjee Gyaltsen
@brief: antibody epitope prediction - standalone version
"""

import os
import sys

from argparse import ArgumentParser

from antibody_epitope_prediction import AntibodyEpitopePrediction
from util import *


def flatten_list(lis_o_lis):
    return [val for lis in lis_o_lis for val in lis]


def calculate(_args):
    if args.list:
        show_method_list()
        sys.exit(0)
    else:
        _prog = ArgumentParser().prog
        if not args.method:
            print(("arguments -m/--method is required\n{}\n"
                  "\nfor detail usage: $ python {} --help".format(msg(), _prog)))
            sys.exit(0)
        elif not args.filename and not args.swissprot:
            print(("either -f/--file or -s/--swissprot argument is expected\n{}\n"
                  "\nfor detail usage: $ python {} --help".format(msg(), _prog)))
            sys.exit(0)

        aep = AntibodyEpitopePrediction()
        output = aep.predict_antibody_epitope(
            method_name=args.method, swissprot=args.swissprot, filename=args.filename,
            window_size=args.window_size)

        print_chart_table(output)

        if args.path:
            generate_plot(args.method, output, args.path)


def is_valid_file(parser, arg):
    """
    Check if arg is a valid file that already exists on the file system.
    """
    arg = os.path.abspath(arg)
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


def msg():
    _prog = ArgumentParser().prog
    usage = """python {} [-h] -m [METHOD_NAME] [options] [SWISSPROT/FILE]""".format(_prog)

    return usage


def csv_list(s):
    return s.split(",")


def get_parser():
    """Get parser object for script predict_antibody_epitope.py."""

    parser = ArgumentParser(description=__doc__, usage=msg())

    req_argument = parser.add_argument_group('required arguments')

    req_argument.add_argument("-m", "--method",
                        dest="method",
                        help="select a method from available method options")

    req_argument.add_argument("-s", "--swissprot",
                        dest="swissprot",
                        help="use when the input is a SwissProt ID (default=file)")

    req_argument.add_argument("-f", "--file",
                        dest="filename",
                        type=lambda x: is_valid_file(parser, x),
                        metavar="FILE",
                        help="a file containing a list of sequence(s)",)

    parser.add_argument("-w", "--window",
                        dest="window_size",
                        default=False,
                        help="sets window size if specified (default=method specific. Eg: 6,7,...)")

    parser.add_argument("-l", "--list",
                        action="store_true",
                        help="show all available method options.")

    parser.add_argument("--plot",
                        dest="path",
                        help="generate a plot.")

    parser.add_argument('--version', action='version', version='%(prog)s v3.0')

    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    calculate(args)
