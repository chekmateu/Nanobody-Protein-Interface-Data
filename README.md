# Nanobody-Protein-Interface-Data

The code to generate the data files are in DataPreperation.ipynb. Just run all the cells which should fetch the data from
SabDab and download all the valid PDB files for nanobody antigen crystal structures. Then it generates the data files named ScanNetData*.txt.

Most python versions above 3 should work but the one I used was 3.10.11. Just run - pip install -r requirements.txt.

In chance that the metadata.csv files aren't installed on the first run just rerun the createMetadata() from getPDBs.py.