#!/usr/bin/env python2.7

import os
import sys
import argparse
import pandas as pd
from collections import OrderedDict

# def get_args():
#     """Parse command line arguments"""
#     desc = ("Part of the wrapper tools for the IEDBtools suite. "
#             "This script can be piped to in order to convert the dataframe "
#             "To a set of attribute files for use with UCSF Chimera")
#     epi = ("Convert wrapper output to UCSF Chimera attribute files "
#            "This enables 'Render by Attribute' functionality for immunogenicity.")
#
#     try:
#         parser = argparse.ArgumentParser(description=desc,
#             epilog=epi,
#             prog=__file__)
#         # parser.add_argument(
#         #     "-v",
#         #     "--verbose",
#         #     default="1",
#         #     choices=["0", "1", "2"],
#         #     help="Change the verbosity of output (info/debug/errors)."
#         # )
#         parser.add_argument(
#             "input",
#             type=argparse.FileType('r'),
#             default=sys.stdin,
#             action="store",
#             help="A positional argument which is a stream from STDIN."
#         )
#         parser.add_argument(
#             "outpath",
#             action="store",
#             help="File path to place created output files."
#         )
#         if len(sys.argv) == 1:
#             parser.print_help(sys.stderr)
#             sys.exit(1)
#
#     except NameError:
#         sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")
#         sys.exit(1)
#
#     return parser.parse_args()

def pad(lst, lamount, ramount, val):
    """Pad a list to the left (prepend) or right (extend) by a specified amount
       with a specific value/entry.
    """
    if lamount > 0:
        lst = ([val] * lamount) + lst
    if ramount > 0:
        lst = lst + ([val] * ramount)

    return lst



def attfile_writer(attribute, sequence, position, score):
    """Write a file containing the sequence-position-specific attributes in an
       easy to parse format for reading in to chimera.

       Note, chimera attribute names must start lower-case, and cannot contain 'unusual' characters
    """

    # Chimera is very strict on file content, so sanitise the attrib name:
    if "-" in attribute:
        attribute = attribute.replace("-", "")
    if any(char.isdigit() for char in attribute):
        raise ValueError("Digits are not allowed in attribute definitions.")
    # Ensure the first character is lowercase
    if attribute[0].isupper():
        attribute = attribute[:1].lower() + attribute[1:]

    container = OrderedDict()
    container['AttributeName'] = attribute
    container['Sequence'] = sequence
    container['Position'] = list(position)
    # Pad the score with leading and trailing zeros to keep scores in register with position
    tmp_list = pad(list(score), container['Position'][0]-1, 0, None)
    # Use 2 separate rounds of left and right padding so that the correct length for the subsequent
    # padding can be obtained.
    container['Score'] = pad(tmp_list, 0, len(sequence)-len(tmp_list), None)


    return container

def pep2seq(column):
    return "".join([f[0] for f in column] + [column[-1][1:]])


def main():

    df = pd.read_csv(
        sys.stdin,
        delim_whitespace=True, # sep='\s+' would be equivalent
        header=1,
        index_col=0
    )
    container = OrderedDict()
    container["Sequence"] = pep2seq(list(df['Peptide']))
    container["Position"] = list(df["Position"])

    for method in df.columns[4:]:
        attribute = method
        if "-" in attribute:
            attribute = attribute.replace("-", "")
        if any(char.isdigit() for char in attribute):
            raise ValueError("Digits are not allowed in attribute definitions.")
        # Ensure the first character is lowercase
        if attribute[0].isupper():
            attribute = attribute[:1].lower() + attribute[1:]
        attribute += "Score"

        tmp_list = pad(list(df[method]), container['Position'][0]-1, 0, None)
        container[attribute] = pad(tmp_list, 0, len(container["Sequence"])-len(tmp_list), None)

    print(container)



if __name__ == "__main__":
    main()
