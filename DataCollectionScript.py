from Bio.PDB import *
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import pandas as pd
import warnings
import csv
import os

from SAbDab_downloader import *

#Computes Binding Region Residues

def getBindingRegionResidues(chain1, chain2):
    # Convert the chain objects to lists of atoms
    atoms1 = list(chain1.get_atoms())
    atoms2 = list(chain2.get_atoms())
    
    # Create a NeighborSearch object for chain1
    ns1 = NeighborSearch(atoms1)
    
    # Find the residues in chain2 that are within 4 angstroms of chain1
    nearby_residues2 = []
    residue = None
    neighbors = None
    for atom in atoms2:
        neighbors = ns1.search(atom.coord, 4, level = 'A')
        for neighbor in neighbors:
            residue = neighbor.get_parent()
            if residue not in nearby_residues2:
                nearby_residues2.append(residue)
    
    # Create a NeighborSearch object for chain2
    ns2 = NeighborSearch(atoms2)

    # Find the residues in chain1 that are within 4 angstroms of chain2
    nearby_residues1 = []
    residue = None
    neighbors = None
    for atom in atoms1:
        neighbors = ns2.search(atom.coord, 4, level = 'A')
        for neighbor in neighbors:
            residue = neighbor.get_parent()
            if residue not in nearby_residues1:
                nearby_residues1.append(residue)

    #nearby_residues are on opposite chains
    output = {chain1.get_id():nearby_residues2, 
              chain2.get_id():nearby_residues1}
    
    return output

def getEpitopeParatopeInteractions(pdbFile, chainPairings, header, path = None, outputFile = None):
    warnings.simplefilter('ignore', PDBConstructionWarning)
    parser = PDBParser()
    structure = parser.get_structure(id = pdbFile, file = pdbFile)
    '''
    Calculate all chain interactions in the file also store chains for later 

    Antigen Chain are everything not in chain Pairings and vice versa for antibodies/nanobodies.
    Then make pairs based on just Fv antigen interactions 
    '''
    allChains = [c.get_id() for c in structure.get_chains()]
    Fvs = chainPairings
    antigenCs = [c for c in allChains if c not in chainPairings]

    pairs = []
    for c1 in Fvs:
        c1 = [c for c in structure.get_chains() if c.get_id() == c1][0]
        for c2 in antigenCs:
            c2 = [c for c in structure.get_chains() if c.get_id() == c2][0]
            pairs.append([
                c1, c2
            ])

    ResidueInContact = []
    for pair in pairs:
        ResidueInContact.append(getBindingRegionResidues(pair[0], pair[1]))
    
    # Remove empty entries
    ResidueInContact = [item for item in ResidueInContact if len(item[list(item.keys())[0]]) != 0]

    # Create csv file

    # Consolidate Binding residue sequence IDs with their chain.
    BindingDict = {}
    for pair in ResidueInContact:
        for k, v in pair.items():
            if k not in list(BindingDict.keys()):
                BindingDict[k] = []
            for residue in v:
                BindingDict[k].append(residue.get_id()[1])

    # Check to see if all chains are in BindingDict to avoid KeyError later on. Create an empty entry as that chain has no binding region.
    for chain in allChains:
        if chain not in list(BindingDict.keys()):
            BindingDict[chain] = []

    ResolvedSeq = {}
    for model in structure:
        for chain in model:
            ResolvedSeq[chain.get_id()] = []
            for residue in chain:
                ResolvedSeq[chain.get_id()].append(residue)

    data = []
    for chain, residues in ResolvedSeq.items():
        for residue in residues:
            label = 1 if residue.get_id()[1] in BindingDict[chain] else 0
            data.append([
                chain, residue.get_id()[1], residue.get_resname(), label
            ])
    
    # Remove water
    data = [row for row in data if row[2] != 'HOH']

    # Write file
    if outputFile != None:
        path = os.path.join(path, outputFile)
        with open(path, 'w', newline = '') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow([header])
            csvwriter.writerows(data)

    return data
    
def compileData(pdb_id, chainPairings, header, output = False, removePDB = False, type = 'Nanobody'):
    '''
    Code to compile the data.
    '''
    try:
        getpdb(pdb_entry=pdb_id, out_path = 'pdbFiles')
        if output:
            getEpitopeParatopeInteractions(
                pdbFile = f'pdbFiles/{pdb_id}.pdb',
                chainPairings = chainPairings,
                header = header,
                path = f"Data/{type}",
                outputFile = f'{pdb_id}_binding.csv'
                )
        else:
            getEpitopeParatopeInteractions(
                pdbFile = f'pdbFiles/{pdb_id}.pdb',
                header = header,
                chainPairings = chainPairings
                )

    except Exception as error:
        print(f"PDB ID - {pdb_id} errored with - {error}")

    if removePDB:
        os.remove(f'pdbFiles/{pdb_id}.pdb')

if __name__ == "__main__":

    '''
    Current version only works on nanobodies since the code can't handle different VH and VL regions
    ''' 

    #Check if directories exist
    if os.path.exists("Data") == False:
        os.makedirs("Data/Antibody")
        os.mkdir('Data/Nanobody')

        from getPDBs import *
        createMetadata()

    if os.path.exists("pdbFiles") == False:
        os.mkdir('pdbFiles/')

    assert(os.path.exists("Data"))
    assert(os.path.exists("pdbFiles"))
    
    nanobody_metadata = pd.read_csv("Data/Nanobody/Nanobody_metadata.csv")
    pdbs = nanobody_metadata['PDB'].to_list()
    for pdb in pdbs:

        # Store the full chain pairings data
        row = nanobody_metadata.query(f"PDB == '{pdb}'")['Chain Pairings']

        # Get just the chain pairings
        chainPairings = row.str.findall(r"VH: (\w)").item()

        compileData(pdb, chainPairings = chainPairings, header = row.item(), output = True, removePDB = True)
    
'''
To Use this file compile data will need:
 - The pdb file containing the antibody antigen.
 - chainPairings is the nanobody chains
 - header can be left as any string
 - output bool
 - removePDB bool
'''