from __future__ import print_function
import sys
from collections import OrderedDict


import chimera
from chimera import openModels
from chimera import Molecule
from chimera import runCommand as Rc


# Read in the values for the attributes
# glob a folder of attribute files?

# Split the model
rc("split")


# Iterate each model, match sequence against it
data = []
with open(sys.argv[1], 'r') as fh:
    for line in fh:
    data.append(eval(line))
# OrderedDict([('AttributeName', 'parker'),
#              ('Score', [2.867, 2.9219999999999997,...
#              ('Sequence', 'MSTSTSQIAVEYPIPVYRFIVSVGDEKIPFNSVSGLDISYDTIEYRDGVGNWFKMPGQSQSTNITLRKGVFPGKTELFDWINSIQLNQVEKKDITISLTNDAGTELLMTWNVSNAFPTSLTSPSFDATSNDIAVQEITLMADRVIMQAV'),
#              ('Position', [5, 6, 7, 8, 9, 10, 11, 12, 13...

# TODO:
# Return to this block and alter it to enable non-exact sequence matching.
# if matching subsequences, the following could be used:
# for target_seq, vals in data...:
#     try:
#        i = str(seq).index(subseq)
#     except ValueError:
#        continue
#     for offset, val in enumerate(vals)
#        r = seq.residues[i+offest]
#        if r:  # catch missing structural residues
#            setattr(r, attribute, vale)

for m in chimera.openModels.list(modelTypes=[chimera.Molecule]):
    for entry in data:
        if m.sequences()[0] == entry['Sequence']:
            for i, r in enumerate(m.residues):
                if i != entry['Position'][0]
                    continue
                else:
                    if r:
                        setattr(r, entry['AttributeName'], entry['Score'][i])

rc("rangecol {} min white mid white max red".format(entry['AttributeName'])

def worms(residue, attribute):
    residue.ribbonDrawMode = chimera.Residue.Ribbon_Round
    residue.ribbonDisplay = True
    if getattr(r, attribute, None) is None:
        rad = 0.05
    else:
        rad = r.attribute/2.0 + 1.0
        residue.ribbonStyle = chimera.RibbonStyleWorm([rad])
