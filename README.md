# Flib: A Fragment Library Generation Software

Saulo de Oliveira - February - 2015
Current Version: 1.01

## 1. INSTALLATION

Change to the directory where you have extracted the contents of "Flib.tgz"
and type:

```sh
make
```

This will generate the executable *Flib*.


## 2. DEPENDENCIES

Flib does not require any additional software to be executed. 

However, Flib requires four input files and accepts a fifth optional file:

1- PDB\_ID.fasta.txt : The fasta sequence of the target.
2- PDB\_ID.fasta.ss  : The predicted secondary structure of the target as output by PSIPRED
3- PDB\_ID.spXout    : The predicted torsion angles of the target as output by SPINE-X.
4- PDB\_ID.hhr       : Threading hits as generated by HHBlits.
5- PDB\_ID.homol     : A text file listing all the homologs to the target (one PDB\_ID per line).

An example of each of these files for protein 1AIU can be found in the folder *examples*.

On top of the required input files, Flib also requires a local version of the Protein Data Bank (PDB). All the ".pdb" files should be grouped under the
same folder and their IDs should be in lower case. 

Local versions of the following software can be incorporated into Flib's pipeline:

1- PSIPRED
2- HHSearch/HHBlits
3- SPINE-X
4- BLAST

The script *install\_dependencies.py* can be used to install Flib's dependencies.

To generate a SAINT2-compliant library, the Python script "process\_new.py" should be executed. This script requires Python 2.6 or higher and the Biopython Module.

## 3.CONFIGURATION

In order to configure Flib, alter the path names on the file "pipeline.sh".

Make sure to provide the correct paths to Flib and to the local version of the PDB.

If using any local versions of the software described in the previous section,
make sure to alter the respective flags to "true" as instructed by the comments 
on the script "pipeline.sh".

## 4.RUNNING FLIB

To generate fragment libraries using Flib, you should use the script "pipeline.sh" as follows:

```sh
./pipeline 1AIU
```

The following example should generate the fragment library file "1AIU.lib".

## 5.TROUBLESHOOTING

Contact sauloho@gmail.com for any problems you may experience using Flib.

