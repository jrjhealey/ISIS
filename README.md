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
 - It may try to install `numpy`. This can cause issues with chimera's older module, so if issues are encountered, downgrade numpy, or install without dependencies (`/Applications/Chimera.app/Contents/Resources/bin/python2.7 -m pip install biopython --no-deps`)


# Workflow

## Linear B-Cell epitope prediction

HMMs combined with a propensity scale method that allows for thresholding and scores to be assigned to the antigenicity of a particular sequence.

This tool implements several methods:

## Chou-Fasman
 - Chou PY, Fasman GD. Prediction of the secondary structure of proteins from their amino acid sequence. Adv Enzymol Relat Areas Mol Biol. 1978;47:45-148.
 - http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=364941

Description: The rationale for predicting turns to predict antibody epitopes is based on the paper by Pellequer et al, Immunology Letters, 36 (1993) 83-99. Instead of implementing the turn scale of that paper which has some non-standard properties, we decided to use the Chou and Fasman scale which is commonly used to predict beta turns as described in the reference link above.
