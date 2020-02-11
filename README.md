# ISIS
Nothing to do with the Islamic State; a repository to hold data and workflows for In Silico Immunogenicity Studies based on the tools from IEDB


# Todo
- Interface with eCIS database for orthologue clustering analyses
- Interface with chimera for 3D mapping to protein structures
- Alignment engine for orthologue alignments and scores on a per amino basis
- Codon redundancy back translation for optimised sequences.

# Installing

First, install UCSF Chimera according to its standard instructions.

**IMPORTANT: all of what follows is PYTHON2 ONLY.**

Install `biopython` in to Chimera:
 - First, find Chimera's own python binary (`which chimera`). Should look something like `/Applications/Chimera.app/Contents/Resources/bin/` (on Mac).
 - Run `/Applications/Chimera.app/Contents/Resources/bin/python2.7 -m pip install biopython`
  - If `pip` isn't present, try `/Applications/Chimera.app/Contents/Resources/bin/python2.7 -m ensurepip`. Refer to section 3b here for more info: http://www.cgl.ucsf.edu/chimera/docs/ProgrammersGuide/faq.html#q3b.
 - It may try to install `numpy`. This can cause issues with chimera's older module, so if issues are encountered, downgrade numpy, or install without dependencies (`/Applications/Chimera.app/Contents/Resources/bin/python2.7 -m pip install biopython --no-deps`)


# Workflow
The suite contains 2 approaches:

 1. identification of significant peptides and rendering accordingy.
 2. Calculation of all scores and rendering accordingly.

(2) is still a W.I.P. and accuracy isn't guaranteed. It's also very slow. For the direct detection of significant epitopes however, implementation is reasonably solid now.

-----
# Running

*Still a work in progress, much of this will be combined later for simplicity.*

#### Generating sequences
Since PDBs with missing residues make sequence analysis hard, it is advised to open the model in Chimera first, and export the model chains as they appear in the structure like so, from chimera's commandline:

    runscript /path/to/Chains2Fasta.py /path/to/output.fasta

At the moment, the code can't handle multifastas, so they need to be split. This oneliner can do so (requires BioPython):

    python -c "import sys; from Bio import SeqIO;[SeqIO.write(r,r.description.replace(' ','')+'.fasta', 'fasta') for r in SeqIO.parse(sys.argv[1], 'fasta')];" /path/to/input.fasta

It should convert 'in place' (but won't affect the input file).

*#TODO Split chains to individual files directly.*

#### Generating predictions
Run the newly created and split files through the predictor.

For one file:

    python2.7 Epitopes.py -i single.fasta -v 0 > single.epi

    # -w is optional (a large w (around 9), will lead to fewer peptides predicted)
    # -v must be set to 0 when redirecting to a file so as not to contaminate the output (to be fixed in a later version)

For all split files:

    for file in /path/to/files/*.fasta ; do
        python2.7 Epitopes.py -i ${file} -v 0 > "${file%.*}".epi
    done

or, using `parallel`:

    ls /path/to/files/*.fasta | parallel 'python2.7 Epitopes.py -i {} -v 0 > {.}.epi'

Both of these approaches will write `.epi` files (which are just lines of python dictionaries) in to the same location as the input file.

# Rendering in UCSF chimera
First open the structure file you care about (the one that originated the sequences in the first step).

In chimera's commandline:

    runscript /path/to/EpitopeRenderer.py /path/to/file.epi

The code should then run (note this will take a while if you have a lot of chains since its a combinatorial effect).

It will cycle render colours at the end of its execution so you'll know it finished.

# Resource information

## Linear B-Cell epitope prediction

HMMs combined with a propensity scale method that allows for thresholding and scores to be assigned to the antigenicity of a particular sequence.

This tool implements several methods:

## Chou-Fasman
 - Chou PY, Fasman GD. Prediction of the secondary structure of proteins from their amino acid sequence. Adv Enzymol Relat Areas Mol Biol. 1978;47:45-148.
 - http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=364941

Description: The rationale for predicting turns to predict antibody epitopes is based on the paper by Pellequer et al, Immunology Letters, 36 (1993) 83-99. Instead of implementing the turn scale of that paper which has some non-standard properties, we decided to use the Chou and Fasman scale which is commonly used to predict beta turns as described in the reference link above.
