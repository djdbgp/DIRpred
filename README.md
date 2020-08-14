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

## Test run
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
This file shows how the positions of paralogs align relatively to the reference for both MSA and MSTA alignments.  
When multiple positions are shown to correspond to one reference position, the aligned position is the last of the list, with all the previously listed positions being an insertion. Only exception to this is if the reference position is the last. In that case, the matching residue is the first of the array, while the rest are a trailing insertion.
This file also contains the individual conservation scores of the aligning positions.

- **plots/** 
This folder contains several DIRpred score plots, and the individual score distributions

- **coevolution/** 
This folder contains the partial results of the co-evolution analysis.

## Input file
The input file is a tab separated list of inputs, with the following columns
1. type: can be LIGAND, MSTA, MSA, RECEPTOR
2. path: the PATH to the file
3. name: the NAME of LIGAND
4. pdb: the PDB ID of LIGAND

An example can be found in testdata/input_files

## ARGS file
The ARGS file contains in a cleaner way the arguments required when calling the python script. Remember to use the @ symbol in front of the path to indicate that the file contains the arguments of the script. An example arg file is found in testargs/testARGS_egfDIRpred_MIblos4444.arg
A list of the required and optional parameters is found just after this text. What are you waiting for, go and read!

## Parameters list and description

| Param | Description| Required | Default | Example |
|:-----:|:----------:|:--------:|:-------:|:-------:|
|--input-file|File containing the multiple MSA and labels|yes|-|--input-file testdata/input_files/egf_inputfile.txt|
|--output-path|Path to the output folder, created if non-existant|yes|-|--output-path ./test_out/|
|--test-type|Type of test, for this version only DIRP|no|DIRP|--test-type DIRP|
|--reference|Reference protein for DIRpred test|yes|-|--reference EGF|
|--conservation-test|Type of conservation measure used for evolutionary alignment, can be: id, blos, jsdw|no|id|--conservation-test blos|
|--alignment-test|Type of conservation measure used for paralogs alignment, can be: id, blos, jsdw. If not specified, will be the same as --conservation-test|no|""|--alignment-test id|
|--coevolution-test|The measure of coevolution, can be: MI,MIp. MIp is APC corrected MI|no|MI|--coevolution-test MIp|
|--abcd|The individual contributions for each of the 4 scores, expressed as in 1 over x. By default is 4,4,4,4 meaning 1 over 4, four times. The four scores in order are: (I)evolutionary conservation (II)paralogs conservation (III)ligand-receptor coevolution (IV) ligand internal coevolution|no|4,4,4,4|--abcd 3,3,6,6|
