# PDBMD2CD Standalone

This is the standalone version of [PDBMD2CD](https://pdbmd2cd.cryst.bbk.ac.uk/).

It is a command line program to predict the CD spectra from protein atomic coordinates.

## Installation

The directory containing the program should look like this:

```
pdbmd2cd_standalone/
├── README.md
├── config.py
├── pdbmd2cd_sa.py
├── ref_basis.txt
├── ref_cd.txt
├── ref_lincomb_ss.txt
├── ref_lstsq_ss.txt
├── ref_names.txt
├── requirements.txt
└── tests
    ├── 2m6q.pdb
    ├── 2y65.pdb
    ├── 329d.pdb
    ├── 5h04.pdb
    ├── check
    │   ├── 2m6q_0.pdb_cd.txt
    │   ├── 2m6q_1.pdb_cd.txt
    │   ├── 2m6q_2.pdb_cd.txt
    │   ├── 2m6q_3.pdb_cd.txt
    │   ├── 2m6q_4.pdb_cd.txt
    │   ├── 2m6q_5.pdb_cd.txt
    │   ├── 2m6q_6.pdb_cd.txt
    │   ├── 2m6q_7.pdb_cd.txt
    │   ├── 2m6q_8.pdb_cd.txt
    │   ├── 2m6q_9.pdb_cd.txt
    │   ├── 2y65.pdb_cd.txt
    │   ├── 329d.pdb_cd.txt
    │   ├── 5h04.pdb_cd.txt
    │   └── average_cd.txt
    └── out
```

You will need to install DSSP somewhere on your system. For this package the [conda](https://anaconda.org/salilab/dssp) version is used - this can be installed with `conda install -c salilab dssp`. Other versions of DSSP >= 3.0 should work.

You will need Python 3 installed, as well as the following dependencies:

scipy==1.3.1
numpy==1.19.5
networkx==2.3
biopython==1.79

You can install these using `pip install -r requirements.txt` on the command line, or individually using pip or conda.

You will need to set two variables found in `config.py`. These are the path you want the log file to be written to (LOGPATH), and the path/command for the DSSP executable (DSSPPATH).

## Usage

This version of PDBMD2CD only works with PDB files, not mmCIF files.

Multi-model files will work - each model will be analysed seperately. A PDB file will be created for each model with a numeric suffix e.g. XXXX.pdb  -> XXXX_0.pdb, XXXX_1.pdb... etc.

Put all your PDB files in a directory. You can optionally create an output directory as well.

To make the predictions, run this on the command line:

`python /path/to/pdbmd2cd_sa.py inputdir outputdir`

where `inputdir` is the directory containing all your input PDBs, and `outputdir` is the directory you want the prediction output files to be written to.

DSSP will be run on the all the input files first, then the predictions will be made using these .dssp files.

The output files have a header, with each header line starting with #. The predicted CD spectra is below the header from 180-260nm, and all CD signal values are in units of delta epsilon.

If the analysis fails for a file, an output file will still be created but it will not have a predicted spectra, instead it will say "Calculation Failed".  

The log file will be written to the LOGPATH and will show the date/time of the analysis, the number of failed predictions and the number of successful predictions.

If you want to check everything is working correctly, go to the directory where `pdbmd2cd_sa.py` is and run:

`python pdbmd2cd_sa.py ./tests/ ./tests/out`

and compare the output in `./tests/out` to that in `./tests/check`. You will notice that one of the files will fail - this is expected, and happens because it is purely DNA with no protein.

# IMPORTANT

We ask that you periodically send the log file produced by this program to r.w.janes@qmul.ac.uk

## Reporting issues

If you have problems with this software, please email e.drew@qmul.ac.uk or r.w.janes@qmul.ac.uk

## Licensing:

For non-commercial use, this software can be used according to the terms of the TODO

For commercial use, TODO