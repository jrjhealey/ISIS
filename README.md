# ISIS
Nothing to do with the Islamic State; a repository to hold data and workflows for In Silico Immunogenicity Studies based on the tools from IEDB


# Todo
- Interface with eCIS database for orthologue clustering analyses
- Interface with chimera for 3D mapping to protein structures
- Alignment engine for orthologue alignments and scores on a per amino basis
- Codon redundancy back translation for optimised sequences.

# Installing

First, install UCSF Chimera according to its standard instructions.

Install `biopython` in to Chimera:
 - First, find Chimera's own python binary (`which chimera`). Should look something like `/Applications/Chimera.app/Contents/Resources/bin/` (on Mac).
 - Run `/Applications/Chimera.app/Contents/Resources/bin/python2.7 -m pip install biopython`
  - If `pip` isn't present, try `/Applications/Chimera.app/Contents/Resources/bin/python2.7 -m ensurepip`. Refer to section 3b here for more info: http://www.cgl.ucsf.edu/chimera/docs/ProgrammersGuide/faq.html#q3b.
 - It may try to install `numpy`. This can cause issues with chimera's older module, so if issues are encountered, downgrade numpy.


# Workflow

## Linear B-Cell epitope prediction
`http://tools.iedb.org/bcell/help/#Bepipred`

HMMs combined with a propensity scale method that allows for thresholding and scores to be assigned to the antigenicity of a particular sequence.
