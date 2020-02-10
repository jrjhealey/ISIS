import chimera
from chimera import openModels
from chimera import Molecule
from chimera import runCommand as rc

rc("split")

for model.sequences()[0] in openModels.list(modelTypes=Molecule):
    for attribute

# Quick notes

 for each model.seq:
    # In future, use a local version of Bcell and run directly here?
    for each attribute:
        if no significant peptides (None):
            next
        else:
            try:
                index the string (tuple number 2): # 2 is peptide seq
                for residueindex to residueindex+len(peptide):
                    # Could include assert that index == tuple 1 (the start position)
                    setattr(attributename, model.residue, tuple 3)
            except IndexError (no subsequence found):
                next
