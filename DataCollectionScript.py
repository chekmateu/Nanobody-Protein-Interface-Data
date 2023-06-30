from Bio.PDB import *
from Bio.PDB.PDBExceptions import PDBConstructionWarning
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

def getEpitopeParatopeInteractions(pdbFile, outputFile = None):
    warnings.simplefilter('ignore', PDBConstructionWarning)
    parser = PDBParser()
    structure = parser.get_structure(id = pdbFile, file = pdbFile)

    # Calculate all chain interactions in the file also store chains for later
    pairs = []   
    allChains = []
    for model in structure:
        for chain in model:
            allChains.append(chain.get_id())
            for pair in model:
                if chain.get_id() == pair.get_id():
                    continue
                else:
                    pairs.append([
                        chain, pair
                    ])

    ResidueInContact = []
    for pair in pairs:
        ResidueInContact.append(getBindingRegionResidues(pair[0], pair[1]))
    
    # Remove empty entries
    ResidueInContact = [item for item in ResidueInContact if len(item[list(item.keys())[0]]) != 0]

    # Remove intermolecular chain interactions
    compounds = structure.header['compound']
    chains = {k: v['chain'].replace(',', '').upper().split() for k, v in compounds.items()}
    InterChain = []
    for k, v in chains.items():
        for chain in v:
            for copy in v:
                if copy == chain:
                    continue
                else:
                    InterChain.append([
                        chain, copy
                    ])
    ResidueInContact = [item for item in ResidueInContact if list(item.keys()) not in InterChain]

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
        fields = ['chain', 'AA#', 'AA', 'Binding']
        path = os.path.join('Data', outputFile)
        with open(path, 'w', newline = '') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(fields)
            csvwriter.writerows(data)

    return data
    
def compileData(pdb_id, output):
    getpdb(pdb_entry=pdb_id, out_path = 'pdbFiles')
    if output:
        getEpitopeParatopeInteractions(pdbFile = f'pdbFiles/{pdb_id}.pdb', outputFile = f'{pdb_id}_binding.csv')
    else:
        getEpitopeParatopeInteractions(pdbFile = f'pdbFiles/{pdb_id}.pdb')

if __name__ == "__main__":
    compileData('8jlz', output = True)