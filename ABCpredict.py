#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np
import sys


logo = r"""

                                                                                                                         
  ,---.  ,-----.   ,-----.                        ,--.,--.        ,--.   
 /  O  \ |  |) /_ '  .--./ ,---. ,--.--. ,---.  ,-|  |`--' ,---.,-'  '-. 
|  .-.  ||  .-.  \|  |    | .-. ||  .--'| .-. :' .-. |,--.| .--''-.  .-' 
|  | |  ||  '--' /'  '--'\| '-' '|  |   \   --.\ `-' ||  |\ `--.  |  |   
`--' `--'`------'  `-----'|  |-' `--'    `----' `---' `--' `---'  `--'   
                          `--'                                           
                          
University of Pisa - Molecular Modeling and Virtual Screening Laboratory 
https://www.mmvsl.it/wp/
"""

print(logo)

import argparse
parser = argparse.ArgumentParser(description="Machine learning models to predict Caco-2 permeability, BBB permeability and plasma protein binding (expressed as fraction unbound).")
parser.add_argument("--in_csv", "-in", help="Input csv path. Must contain SMILES column.")
parser.add_argument("--endpoint", "-e", help="Endpoint to predict. Available options: [caco2, bbb, fu]")
parser.add_argument("--out_csv", "-out", help="Path of output csv containing predictions.")
args = parser.parse_args()

in_csv, out_csv = args.in_csv, args.out_csv
endpoint = args.endpoint


import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, DataStructs
import numpy as np
import joblib
import sys
import PubChemFingerprints

def fp_as_array(mol, fp_choice, to_numpy=True):
    if mol is None:
        return None
    if fp_choice == "Morgan":
        generator = AllChem.GetMorganGenerator(radius=2, fpSize=1024)
        fp = generator.GetFingerprint(mol)
    elif fp_choice == "RDKit":
        fp = Chem.RDKFingerprint(mol, fpSize=1024)
    elif fp_choice == "PubChem":
        fp = PubChemFingerprints.calcPubChemFingerAll(mol)
        bitstring = "".join([str(x) for x in fp])
        fp = DataStructs.CreateFromBitString(bitstring)
    if to_numpy:
        arr = np.zeros((1,), dtype=int)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    else:
        return fp


df = pd.read_csv(in_csv)
if endpoint == "caco2":
    m = joblib.load("models/caco2_GP_hybrid.dump")
    s = joblib.load("models/caco2_scaler.dump")
    X_fp = np.stack([fp_as_array(Chem.MolFromSmiles(smi), "Morgan", to_numpy=True) for smi in df["SMILES"]])
    desc_list = joblib.load("rdkit_208-desc_list.dump") 
    unscaled_desc = []
    for smi in df["SMILES"]:
        mol = Chem.MolFromSmiles(smi)
        res = {}
        for nm,fn in Descriptors._descList: 
            try:
                val = fn(mol)
            except:
                val = 0
            res[nm] = val
        unscaled_desc.append(np.array([res[key] for key in desc_list]))
    X_scaled_desc = s.transform(np.vstack(unscaled_desc))
    X = np.hstack([X_fp, X_scaled_desc]) 
    df["caco2_predictions"] = m.predict(X)

elif endpoint == "bbb":
    m = joblib.load("models/classBBB_RF_PubChem-FP.dump")
    X = np.stack([fp_as_array(Chem.MolFromSmiles(smi), "PubChem", to_numpy=True) for smi in df["SMILES"]])
    df["bbb_predictions"] = m.predict(X)
    
elif endpoint == "fu":
    m = joblib.load("models/fu_GP_RDKit-desc.dump")
    s = joblib.load("models/fu_scaler.dump")
    desc_list = joblib.load("rdkit_208-desc_list.dump") 
    unscaled_desc = []
    for smi in df["SMILES"]:
        mol = Chem.MolFromSmiles(smi)
        res = {}
        for nm,fn in Descriptors._descList: 
            try:
                val = fn(mol)
            except:
                val = 0
            res[nm] = val
        unscaled_desc.append(np.array([res[key] for key in desc_list]))
    X = s.transform(np.vstack(unscaled_desc))
    df["fu_predictions"] = m.predict(X)
else:
    print(f"Please select one available endpoint: {endpoint} is not among available options. Available endpoints are: [caco2, bbb, fu].")
    sys.exit()

df.to_csv(out_csv, index=False)
print("Results saved as ", out_csv)
