from Bio.PDB import *
from Bio.PDB import Structure, Model, Chain
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.ResidueDepth import ResidueDepth, get_surface, min_dist
import warnings
import numpy as np
import pandas as pd
import csv
import torch
import os

'''
Data obj for a pdb interaction file

-- Improvement #1: The receptor data only needs to be shown once for a given protein.

'''
class Interactions():

    '''
    Split for receptor and ligand.
    pdb = pdb file path,
    VH = variable heavy chain id, 
    PSSM (optional as can be defined later) = Position Specific Scoring Matrix

    The pdb file should only have 1 model
    '''
    def __init__(self, pdb, VH, pssm = None):
        warnings.simplefilter('ignore', PDBConstructionWarning)

        parser = PDBParser()
        self.model = parser.get_structure(id = pdb, file = pdb)[0]

        self.pdbID = pdb[-8:-4]

        self.receptor = parser.get_structure(id = pdb, file = pdb)[0]
        # Number of interactions is dependent on amount of ligands  
        self.ligands = parser.get_structure(id = pdb, file = pdb)[0]      
        for c in VH:
            self.receptor.detach_child(c)

        for chain in self.model:
            if chain.get_id() not in VH:
                self.ligands.detach_child(chain.get_id())

        self.receptorAA = {c.get_id(): [r for r in c] for c in self.receptor.get_chains()}
        self.ligandsAA = {c.get_id(): [r for r in c] for c in self.ligands.get_chains()}
        '''
        Structured data is a dictionary of receptor chain data, ligand chain data.
        Ligand data is structured with VH and VL keys if VL keys specified. Also a new key entry will be
        created per ligand in file.
        '''
        # Total length of feature = 67 per amino acid in chain (refered to as n). Keys are chain IDs
        self.interactions = {}

    # Simple utility function to get a chain from the structure
    def getChain(self, id):
        for c in self.model:
            if id == c.get_id():
                return c
        return 'Chain not in model'

    '''
    Generates a list of xyz positions for residues based on Carbon alpha position.
    There will be an interaction per ligand. 
    '''
    def calcPercentAAKNearestNeighborhood(structure, residue):
        return

    def generateLabels(self, ScanNet = False, ReturnTensorDict = False):
        '''
        This function generates a 1 for binding residue and 0 for non-binding per ineraction on file
        '''
        # 3 to 1 Amino acid translations
        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

        def getBindingRegionResidues(chain1, chain2, cutoff_distance = 4.0):
            # Create NeighborSearch objects for each chain
            ns1 = NeighborSearch(list(chain1.get_atoms()))

            # Find the residues in chain2 that are within cutoff_distance angstroms of chain1
            nearby_residues2 = set()
            for atom in chain2.get_atoms():
                neighbors = ns1.search(atom.coord, cutoff_distance, level='A')
                for neighbor in neighbors:
                    residue = neighbor.get_parent()
                    nearby_residues2.add(residue.get_id()[1])

            ns2 = NeighborSearch(list(chain2.get_atoms()))

            # Find the residues in chain1 that are within cutoff_distance angstroms of chain2
            nearby_residues1 = set()
            for atom in chain1.get_atoms():
                neighbors = ns2.search(atom.coord, cutoff_distance, level='A')
                for neighbor in neighbors:
                    residue = neighbor.get_parent()
                    nearby_residues1.add(residue.get_id()[1])

            # Combine the results and return as a dictionary
            output = {chain1.get_id(): list(nearby_residues2), chain2.get_id(): list(nearby_residues1)}
            return output
        
        labels = {}

        # Each chain in self.ligands in a nanobody. First chain is on ligand. Second is on Receptor.
        for cl in self.ligands:
            for cr in self.receptor:
                BindingRegionResidues = getBindingRegionResidues(cl, cr)
                # Remove Empty Entries
                if len(BindingRegionResidues[cl.get_id()]) != 0 or len(BindingRegionResidues[cr.get_id()]) != 0:
                    for k, v in BindingRegionResidues.items():
                        if k not in list(labels.keys()):
                            labels[k] = []
                        for resid in v:
                            labels[k].append(resid)
        '''
        Data Generator for ScanNet PPI Interface Predictor
        - Generates a new entry for relevent chains in the residue
        '''
        if ScanNet:
            with open('Data/ScanNetData.txt', 'a') as data:
                for chainID, resids in labels.items():
                    data.writelines(f'>{self.pdbID}_{chainID}\n')
                    for r in self.getChain(chainID):
                        label = 1 if r.get_id()[1] in resids else 0
                        if r.get_resname() in d.keys():
                            data.writelines(
                                f'{chainID} {r.get_id()[1]} {d[r.get_resname()]} {label}\n'
                            )
                        else:
                            print(f'Omited from {self.pdbID} - {chainID} {r.get_id()[1]} {r.get_resname()} {label}')

        return labels

    '''
    Consider replacing self.interactions for returning the data so that PPIDatabase can use that to package the data.
    '''
    def CompileData(self):
        # Kyte Doolittle Hydrophobicity Scale
        kd = {'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
              'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
              'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
              'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2}
        
        # Grantham Scale for Amino Acid Polarity
        GranthamScale = {'ALA': 8.100, 'ARG': 10.500, 'ASN': 11.600, 'ASP': 13.000, 'CYS': 5.500,
                         'GLN': 10.500, 'GLU': 12.300, 'GLY': 9.000, 'HIS': 10.400, 'ILE': 5.200,
                         'LEU': 4.900, 'LYS': 11.300, 'MET': 5.700, 'PHE': 5.200, 'PRO': 8.000,
                         'SER': 9.200, 'THR': 8.600, 'TRP': 5.400, 'TYR': 6.200, 'VAL': 5.900}  
        
        # 3 to 1 Amino acid translations
        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
            
        # Compute surface
        receptorSurface = get_surface(self.receptor)

        receptorData = []
        for c, rs in self.receptorAA.items():
            for i, r in enumerate(rs):  
                if "CA" in r:
                    xyz = np.array(r["CA"].get_coord())
                    ResDepth = np.array([min_dist(tuple(r["CA"].get_coord().tolist()), receptorSurface)])
                    Hydrophobicity = np.array([kd[d[r.get_resname()]]])
                    Polarity = np.array([GranthamScale[r.get_resname()]])
                    receptorData.append(np.concatenate([xyz, ResDepth, Hydrophobicity, Polarity], axis = 0))
                # Reset Params
                xyz, ResDepth, Hydrophobicity, Polarity = None, None, None, None
                    
        # Iterate through ligands
        for i, c in enumerate(self.ligands):
            # Get ligand Surface
            ligandSurface = get_surface(c)

            ligandData = []
            for r in c:
                if "CA" in r:
                    xyz = np.array(r["CA"].get_coord())
                    ResDepth = np.array([min_dist(tuple(r["CA"].get_coord().tolist()), ligandSurface)])
                    Hydrophobicity = np.array([kd[d[r.get_resname()]]])
                    Polarity = np.array([GranthamScale[r.get_resname()]])
                    ligandData.append(np.concatenate([xyz, ResDepth, Hydrophobicity, Polarity], axis = 0))    
                xyz, ResDepth, Hydrophobicity, Polarity = None, None, None, None       

            #Create new Entry
            self.interactions[f'{i}: Ligand Chain - {c.get_id()}'] = {
                'Receptor Data': receptorData,
                'Ligand Data': [ligandData],
                'Interface': []
            }
            ligandData = []