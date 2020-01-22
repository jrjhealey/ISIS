from __future__ import print_function
import sys
import re
from collections import OrderedDict

import chimera
from chimera import openModels
from chimera import Molecule
from chimera import runCommand as rc

from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

#### Definitions of helper functions

def alignment2cigar(ref, qry):
    """Reconstruct a CIGAR string from a pair of pairwise-aligned sequence strings"""
    if len(ref) != len(qry):
        raise Exception('Unequal length alignment found.')
    cigar = []
    # Iterate both strings by character
    for i in range(len(ref)):
        r, q = ref[i], qry[i]
        # Set the correct character
        op = '-' if (r == "-" and q == "-") else '=' if r == q else \
             'I' if r == '-' else 'D' if q == '-' else 'X'
        # Enumerate the CIGAR integers
        if len(cigar) > 0 and cigar[-1][1] is op:
            cigar[-1][0] += 1
        else: cigar.append([1, op])
    return "".join(map(lambda x: str(x[0]) + x[1], cigar))


def expand_cigar(cigar):
    """Given a CIGAR string, expand it to it's letter only representation.
       This code uses a modified definition of CIGAR strings that permits dual
       gap columns as (-).

       M/=  ->  MATCH
       I    ->  Insertion to the reference
       D    ->  Deletion to the reference
       N    ->  Skipped region from the reference (currently unused)
       S    ->  Soft clipping (currently unused)
       H    ->  Hard clipping (currently unused)
       P    ->  Padding (silent deletion from padded reference) (currently unused)
       =    ->  Sequence match
       X    ->  Sequence mismatch
       -    ->  Dual gap column (missing data)

       e.g. 4=1M2D3=  ->  ====MDD===

       e.g. '1-2D13X3D1X1=11X1=7X1=5X1=2X1=1X5I6-' ->
            '-DDXXXXXXXXXXXXXDDDX=XXXXXXXXXXX=XXXXXXX=XXXXX=XX=XIIIII------'
    """
    regex = re.compile(r'(\d+)([-=IDMX])')
    return "".join([int(num) * op for num, op in
                   [submatch for submatch in re.findall(regex, cigar)]])


def cigar_guided_score_alignment(cigar, scores_list):
    """Taking a CIGAR representation of a pairwise alignment, expand and
       recapitulate the 'moveset', and apply corresponding list manipulations

       Behaviour map:
       D   ->   Sequence deletion: drop the score value
       M/= ->   Sequence match: propagate the score value
       I   ->   Sequence insertion: insert dummy data (0) (???)
                 (It should be unlikely the structure has more sequence)
       X   ->   Sequence mismatch: insert dummy data (0)
    """
    new_scores = []
    for i, (j, k) in enumerate(zip(expand_cigar(cigar), scores_list)):
        if j == "=":
            new_scores.append(scores_list[i])
        elif j == "D" or j == "-":
            new_scores.append("-")
        elif j == "I" or j == "X":
            new_scores.append(0)

            # Find out how chimera treats missing data for attributes. May need to set "" instead of 0
    return new_scores
# Read in the values for the attributes

# Read data from file of OrderedDicts as literals
data = []
with open(sys.argv[1], 'r') as fh:
    data = [eval(line) for line in fh]
# Consider changing to glob a folder of files once one is working

# General datastructure:
# OrderedDict([('AttributeName', 'parker'),
#              ('Score', [2.867, 2.9219999999999997,...
#              ('Sequence', 'MSTSTSQIAVEYPIPVYRFIVSVGDEK...
#              ('Position', [5, 6, 7, 8, 9, 10, 11, 12, 13...


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
                if offset
                    for offset, val in enumerate(entry["Score"]):
                        r = seq.residues[i+offset]
                        setattr(r, entry["AttributeName"], val)

                # Might be able to do away with the block above if the alignment strategy is solid

                # Else, align the sequence of the MODEL to the query ()
                # Treat the entry seq as reference, model seq as query
                alignment = pairwise2.align.globalds(entry["Sequence"], modelseq,
                                         blosum62, -10, -0.5) # may need to optimise alignment params

                # Next, construct a CIGAR string for the alignment to apply the same set of 'moves'
                # to the scores arrays.
                # Based around the behaviour that a match (=) should set the attribute to the corresponding
                # residue. A mismatch should raise a warning (M) and not set an attribute  but probably not break.
                # A gap should cause the corresponding value to be dropped from the array.
                # Afterwards, collapse the sequence (remove all gaps) and then it should be the same length as the
                # CIGAR-altered score array.

                cigar = alignment2cigar(alignment[0][0], alignment[0][1]) # Uses the best alignment
                new_scores = cigar_guided_score_alignment(cigar, entry["Score"])
                # Finally, drop all the deleted regions
                entry["StructureScore"] = list(filter(lambda a: a != "-", new_scores))
                # For debugging:
                # print("\n".join([alignment[0][0],expand_cigar(cigar),
                #                 alignment[0][1], "".join([str(x) for x in new_scores])]))
                # Doesn't play nice with floats

                for score, res in zip(entry["StructureScore"], model.residues):
                    setattr(r, entry["AttributeName"], val)

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
