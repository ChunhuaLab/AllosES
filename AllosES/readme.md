# **AllosES: an ensemble model for Protein Allosteric site prediction**




The AllosES program is a method for protein allosteric site prediction. This program can be run on Linux system. 

---



## Preparation
Before running the program, the pre-processed functional chain of the target protein is used as an input file. The pre-processing process is to remove water molecules and ligands from the target protein. Also, several software and packages need to be installed to ensure that the program runs properly.

AllosES uses the following dependencies:
- python 3.6
- numpy
- pandas
- joblib
- ProDy

## Example
An allosteric protein with PDB ID "5XZR" is used as an example to show the process. This PDB file is 5XZR.pdb. Only the protein functional chain is preserved.


## How to run
**Step 1:** Run Fpocket program, which can be downloaded from https://fpocket.sourceforge.net/. Follow the instructions to compile the executables.

**Step 2:** Run the DSSP program to obtain the secondary structure file (.dssp) for the target protein.

**Step 3:** Go to the website (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) and generate an (.asn) file.

**Step 4:** Put the PDB file, secondary structure file, the (.asn) file and AllosESmain.py in the same directory with other codes.

**Step 5:** Use the following bash command to run AllosES.


```
python AllosESmain.py --PDBID [pdbid] --CHAIN [chain]
```

> for example:
```
python AllosESmain.py --PDBID 5XZR --CHAIN A
```

> The parameter [pdbid] is the PDB file name of the allosteric protein.

> The parameter [chain] is the functional chain of the target protein.

Then, AllosES program will perform the feature extraction and prediction process, which will take some time, with the calculation of the net transfer entropy feature taking a little more time.

The predicted results will be shown in the output of the AllosESmain.py program, and the pocket ranking can be seen on the bash command interface.


---

## Help
For any questions, please contact us at chunhuali@bjut.edu.cn.
