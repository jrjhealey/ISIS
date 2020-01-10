from __future__ import print_function
import sys
from collections import OrderedDict


import chimera
from chimera import openModels
from chimera import Molecule
from chimera import runCommand as rc

from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

# May not need:
# def levenshtein_distance(s1, s2):
#     if len(s1) < len(s2):
#         return levenshtein_distance(s2, s1)
#
#     # len(s1) >= len(s2)
#     if len(s2) == 0:
#         return len(s1)
#
#     previous_row = range(len(s2) + 1)
#     for i, c1 in enumerate(s1):
#         current_row = [i + 1]
#         for j, c2 in enumerate(s2):
#             insertions = previous_row[j + 1] + 1
#             deletions = current_row[j] + 1
#             substitutions = previous_row[j] + (c1 != c2)
#             current_row.append(min(insertions, deletions, substitutions))
#         previous_row = current_row
#
#     return previous_row[-1]

# Read in the values for the attributes
# glob a folder of attribute files?

# Read data from file of OrderedDicts as literals
data = []
with open(sys.argv[1], 'r') as fh:
    for line in fh:
        print(line)
        data.append(eval(line))

# General datastructure:
# OrderedDict([('AttributeName', 'parker'),
#              ('Score', [2.867, 2.9219999999999997,...
#              ('Sequence', 'MSTSTSQIAVEYPIPVYRFIVSVGDEKIPFNSVSGLDISYDTIEYRDGVGNWFKMPGQSQSTNITLRKGVFPGKTELFDWINSIQLNQVEKKDITISLTNDAGTELLMTWNVSNAFPTSLTSPSFDATSNDIAVQEITLMADRVIMQAV'),
#              ('Position', [5, 6, 7, 8, 9, 10, 11, 12, 13...


assignments = {}

# Split the model
#rc("split")

# For each model
for model in openModels.list(modelTypes=Molecule):
    # For each sequence in each model
    for modelseq in model.sequences():
        # For each entry in the inputfile
        for entry in data:
            # First test for exact sequence matches
            if entry["Sequence"] == modelseq:
                # If exact, set residues directly
                # and add an entry to the dict of associations
                assignments[model.name] = entry["Score"]
                for score, res in zip(entry["Score"], model.residues):
                    setattr(res, entry["AttributeName"], score)
            else:
                # If model sequence is some contiguous subsequence of the query seq
                # Figure out where the subsequence starts (reciprocally)
                # and assign values from that point until the end of the vals
        		try:
        		    offest = str(modelseq).index(entry["Sequence"])
        		    print(m.id, m.subid, m.name, i)
        		except ValueError:
                    try:
            		    offset = str(entry["Sequence"].index(modelseq))
                    except ValueError():
                        continue

                for offset, val in enumerate(entry["Score"]):
                    r = seq.residues[i+offset]
                    setattr(r, entry["AttributeName"], val)
                # Else, align the sequence of the MODEL to the query ()
                # Treat the entry seq as reference, model seq as query
                alignment = pairwise2.align.globalds(entry["Sequence"], modelseq,
                                         blosum62, -10, -0.5)
                # Next, construct a CIGAR string for the alignment to apply the same set of 'moves'
                # to the scores arrays.
                # Based around the behaviour that a match (=) should set the attribute to the corresponding
                # residue. A mismatch should raise a warning (M) and not set an attribute  but probably not break.
                # A gap should cause the corresponding value to be dropped from the array.
                # Afterwards, collapse the sequence (remove all gaps) and then it should be the same length as the
                # CIGAR-altered score array.

                cigar = alignment2cigar(alignment[0][0], alignment[0][1])


def alignment2cigar(ref, qry):
    """Reconstruct a CIGAR string from a pair of pairwise-aligned sequence strings"""
    if len(ref) != len(qry):
        raise Exception('Unequal length alignment found.')
    cigar = []
    # Iterate both strings by character
    for i in range(len(ref)):
        r, q = ref[i], qry[i]
        # Ensure no columns of only gaps
        if r == '-' and q == '-':
            raise Exception('Found a column where both sequences have gaps. Check alignment.')
        # Set the correct character
        op = '=' if r == q else 'I' if r == '-' else 'D' if q == '-' else 'X'
        # Enumerate the CIGAR integers
        if len(cigar) > 0 and cigar[-1][1] is op:
            cigar[-1][0] += 1
        else:
            cigar.append([1, op])

    return "".join(map(lambda x: str(x[0]) + x[1], cigar))

def expand_cigar(cigar):
    """Given a CIGAR string, expand it to it's letter only representation
       e.g. 4=1M2D3=  ->  ====MDD===
    """
    return "".join([int(i) * j for i, j in zip(cigar[::2], cigar[1::2])])

# Alignment strategies:
#  After alignment, if the sequence is a perfect match, other than gaps, simply find the
#  indexes of all gap characters, and drop the corresponding indexes from the score array.
# something like
alignment = pairwise2.align.globalds(modelseq, entry["Sequence"], blosum62, -10, -0.5)
# alignment[0][0] = the first sequence of the first (best or joint best) alignment result (reference)
# alignment[0][1] = the second sequence of the first (best or joint best) alignment result (query)
for i, char in enumerate(alignment[0][0])):
# Iterate the aligned string, if there is a gap char, drop the score value
# this is then added as a new object to the entry, to preserve the original data
    temp_scores = entry["Score"]
    if char == "-":
        del temp_scores[i]
    entry["ReindexedScores"] = temp_scores

# This approach only works if gaps are the ONLY issue with the alignment, which should be fine
# for now, but isn't a robust/future proof approach.

#rc("rangecol parker min white mid white max red") #.format(entry['AttributeName'])
#
# def worms(residue, attribute):
#     residue.ribbonDrawMode = chimera.Residue.Ribbon_Round
#     residue.ribbonDisplay = True
#     if getattr(r, attribute, None) is None:
#         rad = 0.05
#     else:
#         rad = r.attribute/2.0 + 1.0
#         residue.ribbonStyle = chimera.RibbonStyleWorm([rad])
