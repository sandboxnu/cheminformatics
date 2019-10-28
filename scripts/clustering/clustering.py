import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
#dictionary of dictionaries, each key is the smile

#add a key with the tanimoto coefficient
def add_tanimoto_coefficients(smiles):
    for smile in smiles:
        smileInfo = smiles[smile]
        similarities = {}
        for othersmile in smiles:
            othersmileInfo = smiles[othersmile]
            similarities[othersmile] = compare_two_smiles(smileInfo['murcko'], othersmileInfo['murcko'])

            print(compare_two_smiles(smileInfo['murcko'], othersmileInfo['murcko']))
        smileInfo['similarities'] = similarities
    return smiles

#need to get the comparison for each smile it would seem, what is the best way to store
#THAT information
#smile 1 and smile two, are they murcko smiles? or 
def compare_two_smiles(smile1, smile2):
    smile1Ms = Chem.MolFromSmiles(smile1)
    smile2Ms = Chem.MolFromSmiles(smile2)
    fps1 = FingerprintMols.FingerprintMol(smile1Ms)
    fps2 = FingerprintMols.FingerprintMol(smile2Ms)

    return DataStructs.FingerprintSimilarity(fps1, fps2)
    

