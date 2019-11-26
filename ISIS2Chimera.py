#!/usr/bin/env python2.7

import os
import sys
import argparse
import pandas as pd

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

def attfile_writer(method, attribute, match_mode, recipient,
                   position, score):
    """Synthesise the text content for an attribute file as defined here:
    https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/defineattrib/defineattrib.html#attrfile

    Requires a list of positions and their corresponding scores as lists"""

    content = " ".join(["#",
                       "Immunogenicity profile for",
                       "dummy_string",
                       "with method:",
                       method]) + "\n"
    # Chimera is very strict on file content, so sanitise the attrib name:
    if "-" in attribute:
        attribute = attribute.replace("-", "_")
    elif any(char.isdigit() for char in attribute):
        raise ValueError("Digits are not allowed in attribute definitions.")

    content += "attribute: " + attribute.lower() + '\n'
    content += "match mode: " + match_mode + '\n'
    content += "recipient: " + recipient + '\n'
    for p, s in zip(position, score):
        content += "\t".join(["", ":"+str(p), str(s)]) + '\n'

    return content


def main():

    df = pd.read_csv(
        sys.stdin,
        delim_whitespace=True, # sep='\s+' would be equivalent
        header=1,
        index_col=0
    )

    for method in df.columns[4:]:
        #with open('./dummy_file' + method + ".att", 'w') as fh:
            #fh.write(attfile_writer(method,
        print(attfile_writer(method, "Immunogenicity_" + method,
                                    match_mode="1-to-1",
                                    recipient="residues",
                                    position=df['Position'],
                                    score=df[method]))






if __name__ == "__main__":
    main()
