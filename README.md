# DIRpred

Hello fellow reader. Please give the software a try with the bare minimum description that I put here.
If you get really desperate, do not despair and contact me at either stefano.pascarelli1@oist.jp or pasca1989@hotmail.it

If you just want to replicate the paper analysis, that should be easy! You just need to use the test data, test input file, test args, and perform the test run.

## Synopsis
DIRpred is a software to extract conservation and coevolution measures from multiple multiple sequence alignments (paralogs, orthologs and paralogs' orthologs).
These are combined to obtain a prediction score on how important a residue is for paralog functional divergence.
More info on how it works can be found on: https://www.biorxiv.org/content/10.1101/677393v3

## Basic flow
The least headaches way to run the script is to prepare an ARGS file and an INPUT file.
ARGS file contains the parameters of the analysis (i.e. the type of conservation measure, the path to the input file, etc.). More details will come later.
INPUT file has the path to each multiple sequence alignment and a tag to identify the type of MSA. Again, more details later.

## Test run:
> cd src
> python3 dirpred.py @../testargs/testARGS_egfDIRpred_MIblos4444.arg

This will generate a test_output using the data from the previously cited paper (EGF-EGFR ligand receptor system). Let's check the output!

## Test output
- **DIRpred.csv**  
This file contains several info for each reference position. In order of columns:
1. Reference position
2. Reference position with reference amino acid
3. MSA, paralogs multiple sequence alignment conservation
4. MSTA, paralogs multiple structure alignment conservation
5. EVO cons MSA, Evolutionary conservation using the paralogs MSA to match positions #
6. EVO cons MSTA, same as before but using MSTA
7. COEVO receptor MSA, highest co-evolution score of this site and the receptor sites #
8. COEVO receptor MSTA, same as before but using MSTA #
9. COEVO ligands MSA, highest internal co-evolution score, obtained by combining orthologs alignment to the reference using MSA alignment
10. COEVO ligands MSTA, same as before but using MSTA
11. avg DIR score, the DIRpred score, average of MSA and MSTA
12. partial scores, the relative contribution of the four scores to the MSA DIRpred score (first tuple) and the MSTA DIRpred score (second tuple). Used for plotting

*# note: this is an average in all paralogs aligned by MSA or MSTA*

- **DIRpred_sorted.csv** 
Same as previous file, but sorted by the 11th column, the final DIRpred score.

- **ligand_table.csv** 
This file shows how the positions of paralogs align relatively to the reference. When multiple positions are shown to correspond to one reference position, then the aligned position is the last of the array, with all the previous listed position being an unalignable insertion. The only exception is to the last reference position. In that match, the first of the array is the aligning position while the rest are a trailing insertion.
This file also contains the individual conservation scores of the aligning positions.

