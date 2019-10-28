import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
#dictionary of dictionaries, each key is the smile

#add a key with the tanimoto coefficient
def add_tanimoto_coefficients(smiles):
    #figure out how to get rdkit in here
    for smile in smiles:
        print("\n")
        print(smiles[smile]['murcko'])
        similarities = []
        for othersmile in smiles:
            similarities['othersmile'] = compare_two_smiles(smiles[smile]['murcko'], smiles[othersmile]['murcko'])})
        smiles[smile]['similarities'] = "similarities"    
    return smiles

#need to get the comparison for each smile it would seem, what is the best way to store
#THAT information
#smile 1 and smile two, are they murcko smiles? or 
def compare_two_smiles(smile1, smile2):
    smile1Ms = Chem.MolFromSmiles(smile1)
    smile2Ms = Chem.MolFromSmiles(smile2)
    fps1 = FingerprintMols.FingerprintMol(smile1Ms)
    fps2 = FingerprintMols.FingerprintMol(smile2Ms)

    DataStructs.FingerprintSimilarity(fps1, fps2)
    

