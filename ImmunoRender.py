from __future__ import print_function
import sys
import re
import logging
from time import sleep
from collections import OrderedDict

import chimera
from chimera import openModels
from chimera import Molecule
from chimera import runCommand as rc

from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

logging.basicConfig(format="[%(asctime)s] %(levelname)-8s->  %(message)s",
                    level=logging.NOTSET, datefmt='%d/%m/%Y %I:%M:%S %p')
logger = logging.getLogger(__name__)

#### Definitions of helper functions

def worms(residue, attribute):
    """Alter chain thicknesses based on an attribute"""
    residue.ribbonDrawMode = chimera.Residue.Ribbon_Round
    residue.ribbonDisplay = True
    if getattr(r, attribute, None) is None:
        rad = 0.05
    else:
        rad = r.attribute/2.0 + 1.0
        residue.ribbonStyle = chimera.RibbonStyleWorm([rad])

def alignment2cigar(ref, qry):
    """Reconstruct a CIGAR string from a pair of pairwise-aligned sequence strings

       Function modified from @lh3 on StackOverflow
    """
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

class AssociationContainer(object):
    """A class to hold associations between alignments, structures, and sequences.
    """
    def __init__(self, model):
        self.alignments = []
        self.model = model

# Read in the values for the attributes

# Read data from file of OrderedDicts as literals
data = []
with open(sys.argv[1], 'r') as fh:
    data = [eval(line) for line in fh]
# Consider changing to glob a folder of files once one is working

# Split the model

logger.info("Splitting model chains...")
rc("split")
all_associations = []
attributes = set()
# For each model
for model in openModels.list(modelTypes=Molecule):
    logger.info("Operating on model: {}.{} {}".format(model.id, model.subid, model.name))
    # For each sequence in each model (though this should just be one)
    for modelseq in model.sequences():
        # For each entry in the inputfile
        logger.info("Got model sequence: {}".format(modelseq))
        for entry in data:
            attributes.add(entry["AttributeName"])
            association = AssociationContainer(model)
            logger.info("Processing entry: {} {}...{}".format(entry["AttributeName"],
                                                              entry["Sequence"][0:10],
                                                              entry["Sequence"][-10:]))
            association.entry = entry
            # First test for exact sequence matches for speed
            if str(entry["Sequence"]) == str(modelseq):
                # and add an entry to the associations object
                logger.info("Input sequence is an exact match for model sequence...")
                all_associations.append(association)
            else:
                # If model sequence is some contiguous subsequence of the query seq
                # Figure out where the subsequence starts (reciprocally)
                # and assign values from that point until the end of the vals
        		# try:
        		#     offest = str(modelseq).index(entry["Sequence"])
        		# except ValueError:
                #     try:
            	# 	    offset = str(entry["Sequence"].index(modelseq))
                #     except ValueError():
                #         continue
                # if offset:
                #     for i, val in enumerate(entry["Score"]):
                #         r = seq.residues[i+offset]
                #         setattr(r, entry["AttributeName"], val)

                # Might be able to do away with the block above if the alignment strategy is solid
                logger.info("Exact full length or substring matches not detected. Switching to alignment based approach...")

                # Else, align the sequence of the MODEL to the input seq
                # Treat the entry seq as reference (since it will more often than
                # not be the longer of the 2), model seq as query
                alignment = pairwise2.align.globalds(str(entry["Sequence"]), str(modelseq), blosum62, -10, -0.5)
                logger.info("Alignment:\n {}".format(pairwise2.format_alignment(*alignment[0])))
                # Gather all top alignments to then filter out the best for application
                association.alignments.append(alignment[0])
            all_associations.append(association)

# Now associations are assigned, modify scores and apply accordingly
    # Now, take the best of all of the alignments as the correct one for that protein
    # Sort the alignments in descending order to take the top one based on the alignment score
logger.info("Associations completed...")
for container in all_associations:
    if len(container.alignments) > 0:
        logger.info("Alignments detected, manipulating scores to match aligned positions...")
        container.alignments.sort(key=lambda tup: float(tup[2]))
    # Next, construct a CIGAR string for the alignment to apply the same set of 'moves'
    # to the scores arrays.
        logger.info("Creating CIGARs...")

        container.cigar = alignment2cigar(container.alignments[0][0], container.alignments[0][1]) # Uses the best alignment
        temp_scores = cigar_guided_score_alignment(container.cigar, container.entry["Score"])
        # Finally, drop all the deleted regions
        container.entry["StructureScore"] = list(filter(lambda a: a != "-", temp_scores))
        # For debugging:
        # print("\n".join([alignment[0][0],expand_cigar(cigar),
        #                 alignment[0][1], "".join([str(x) for x in new_scores])]))
        # Doesn't play nice with floats
    else:
        container.entry["StructureScore"] = container.entry["Score"]

    for score, res in zip(container.entry["StructureScore"], container.model.residues):
        setattr(res, container.entry["AttributeName"], score)

logger.info("Final assignments:")
for association in set(all_associations): # set() needed because associations are duplicated for some reason...
    print("Associated entry: {} {}...{}".format(association.entry["AttributeName"],
                                               association.entry["Sequence"][0:10],
                                               association.entry["Sequence"][-10:]))
    print("With model: {}.{} {}".format(association.model.id,
                                        association.model.subid,
                                        association.model.name))

logger.info("Attributes set. Cycling render views...")
for attribute in attributes:
    rc("rangecol {} min blue mid white max red novalue #097b097b097b".format(attribute))
    #[worms(res, attribute) for res in model.residues for model in [model for model in openModels.list(modelTypes=Molecule)]
    sleep(1)
