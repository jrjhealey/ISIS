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

Scale:

| A | C	| D |	E |	F	| G |	H |	I |	K |	L |	M |	N	| P	| Q	| R	| S	| T	| V	| W	| Y |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|0.66|1.19|1.46|0.74|0.6|1.56|0.95|0.47|1.01|0.59|0.6|1.56|1.52|0.98|0.95|1.43|0.96|0.5|0.96|1.14|

Description: The rationale for predicting turns to predict antibody epitopes is based on the paper by Pellequer et al, Immunology Letters, 36 (1993) 83-99. Instead of implementing the turn scale of that paper which has some non-standard properties, we decided to use the Chou and Fasman scale which is commonly used to predict beta turns as described in the reference link above.

The Chou-Fasman method approximately predicts segments of a protein which are likely to create Beta-turn motifs.

---

## Emini
 - Emini EA, Hughes JV, Perlow DS, Boger J. Induction of hepatitis A virus-neutralizing antibody by a virus-specific synthetic peptide. J Virol. 1985 Sep;55(3):836-9.
 - http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=2991600

Scale:

| A |	C |	D |	E |	F |	G |	H |	I |	K |	L |	M |	N |	P |	Q |	R |	S |	T |	V |	W |	Y |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|0.49|0.26|0.81|0.84|0.42|0.48|0.66|0.34|0.97|0.4|0.48|0.78|0.75|0.84|0.95|0.65|0.7|0.36|0.51|0.76|

Description: The calculation was based on surface accessibility scale on a product instead of an addition within the window. The accessibility profile was obtained using the formulae Sn = (dn+4+i ) (0.37)-6 where Sn is the surface probability, dn is the fractional surface probability value, and i vary from 1 to 6. A hexapeptide sequence with Sn greater than 1.0 indicates an increased probability for being found on the surface.

The Emini method predicts likely immunogenic peptides by virtue of their solvent accessibility/surface exposure.

---

## Karplus-Schulz

 - Reference: Karplus PA, Schulz GE. Prediction of Chain Flexibility in Proteins - A tool for the Selection of Peptide Antigens. Naturwissenschafren 1985; 72:212-3.
 - https://link.springer.com/article/10.1007/BF01195768
 -
Description: In this method, flexibility scale based on mobility of protein segments on the basis of the known temperature B factors of the a-carbons of 31 proteins of known structure was constructed. The calculation based on a flexibility scale is similar to classical calculation, except that the center is the first amino acid of the six amino acids window length, and there are three scales for describing flexibility instead of a single one.

The Karplus-Schulz method selects potential antigens based on their proposed B-factor flexibilities. Theoretically, more flexible regions are more likely to be solvent exposed and amenable to immune recognition.

## Kolaskar-Tongaonkar

 - Kolaskar AS, Tongaonkar PC. A semi-empirical method for prediction of antigenic determinants on protein antigens. FEBS Lett. 1990 Dec 10;276(1-2):172-4.
 - http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=2430611

Description: A semi-empirical method which makes use of physicochemical properties of amino acid residues and their frequencies of occurrence in experimentally known segmental epitopes was developed to predict antigenic determinants on proteins. Application of this method to a large number of proteins has shown by the authors that the method can predict antigenic determinants with about 75% accuracy which is better than most of the known methods.

Scale:

| A | C | D | E | F | G | H | I | K | L | M | N | P | Q | R | S | T | V | W | Y |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|1.064|1.412|0.866|0.851|1.091|0.874|1.105|1.152|0.93|1.25|0.826|0.776|1.064|1.015|0.873|1.012|0.909|1.383|0.893|1.161|
