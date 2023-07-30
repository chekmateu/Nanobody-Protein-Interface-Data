from Bio.PDB import *
from Bio.PDB import Structure, Model, Chain

import glob as glob
import numpy as np
import pandas as pd
import random
import asyncio
import aiohttp
import os

import Interactions

'''
Create Database out of list of pdbIDs. Filtering is applied later.
'''
class PPIDatabase(Interactions.Interactions):
    def __init__(self, metadata, pdbPath = 'pdbFiles'):

        self.pdbIDs = list(metadata['PDB'])
        self.metadata = metadata

        self.pdbPath = pdbPath

    '''
    Fetches pdb files asynchronously and cleans pdbFiles removing all header content and extraneous content,
    leaving only atomic information.

    The out path for getpdbs() will default to pdbFiles/

    Run with await.
    '''
    async def getpdbs(self, out_path = 'pdbFiles'):
        async with aiohttp.ClientSession() as session:
            for pdb_entry in self.pdbIDs:
                url = "https://opig.stats.ox.ac.uk/webapps/abdb/entries/%s/structure/%s.pdb"%(pdb_entry,pdb_entry)
                out_file = os.path.join(out_path, "%s.pdb"%pdb_entry)
                async with session.get(url) as response:
                    assert response.status == 200
                    with open(out_file, 'wb') as f:
                        while True:
                            line = await response.content.readline()
                            if not line:
                                break
                            if line.startswith(b'ATOM'):
                                if line.split(b' ')[11] != b'H':
                                    f.write(line)
    '''
    Generates ScanNet Data.
    Requires pdbFiles and metadata for chain pairings. Split - [Train%, Test%]
    Example - 
    >8dt8_A
    {Chain ID} {Res ID} {1 Letter AA Code} {Is Binding or not}
    '''
    def generatetxtLabelFile(self, outdir = 'Data', split = [0.9, 0.1]):

        assert(split[0]+split[1] == 1)

        pdbfiles = glob.glob(f'{self.pdbPath}/*.pdb')
        
        print('Generating Labels')
        startNewLine = 0
        for pdb in pdbfiles:
            # Store the full chain pairings data
            pdbID = pdb[-8:-4]
            NbChains = self.metadata.query(f"PDB == '{pdbID}'")['Chain Pairings'].item()
            interactionObj = Interactions.Interactions(pdb = pdb, VH = NbChains)
            if startNewLine <= 10:
                print(f'{pdbID};', end = ' ')
                startNewLine += 1
            else:
                print(f'{pdbID};')
                startNewLine = 0
            interactionObj.generateLabels(ScanNet = True)

        '''
        Split ScanNetData, 90% train 10% test
        '''
        print('Splitting ScanNet Data')
        delimiter = '>'
        with open(f'{outdir}/ScanNetData.txt', 'r') as f:
            data = f.read()
            data = [delimiter + i for i in data.split(delimiter) if i]
            random.shuffle(data)

        def split_train_test(data, test_size=0.1, random_state=None):
            if random_state is not None:
                random.seed(random_state)
                
            num_test_samples = round(len(data) * test_size)
            test_indices = random.sample(range(len(data)), num_test_samples)

            train_list = [data[i] for i in range(len(data)) if i not in test_indices]
            test_list = [data[i] for i in test_indices]

            return train_list, test_list

        train, test = split_train_test(data = data)

        with open(f'{outdir}/ScanNetData_train.txt', 'a') as f:
            f.writelines(train)

        with open(f'{outdir}/ScanNetData_test.txt', 'a') as f:
            f.writelines(test)

        return train, test
