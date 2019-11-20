import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from app.data.smiles import  convert
from rdkit.ML.Cluster import Butina

#add a key with the tanimoto coefficient
def add_tanimoto_coefficients(smiles):
    for smile in smiles:
        smileInfo = smiles[smile]
        similarities = {}
        for othersmile in smiles:
            othersmileInfo = smiles[othersmile]
            similarities[othersmile] = compare_two_smiles(smileInfo['murcko'], othersmileInfo['murcko'])
        smileInfo['similarities'] = similarities
    return smiles

def compare_two_smiles(smile1, smile2):
    smile1Ms = Chem.MolFromSmiles(smile1)
    smile2Ms = Chem.MolFromSmiles(smile2)
    fps1 = FingerprintMols.FingerprintMol(smile1Ms)
    fps2 = FingerprintMols.FingerprintMol(smile2Ms)

    return DataStructs.FingerprintSimilarity(fps1, fps2)
    

def cluster(smile_keys, cutoff=0.5):
    npts = len(smile_keys)
    print(npts)

    result = Butina.ClusterData(smile_keys, npts, cutoff, False, compare_murcko_smiles, False)
    print(result)




def compare_murcko_smiles(smile1, smile2):
    murcko1 = convert(smile1)
    murcko2 = convert(smile2)
    smile1Ms = Chem.MolFromSmiles(murcko1)
    smile2Ms = Chem.MolFromSmiles( murcko2)
    fps1 = FingerprintMols.FingerprintMol(smile1Ms)
    fps2 = FingerprintMols.FingerprintMol(smile2Ms)

    return DataStructs.FingerprintSimilarity(fps1, fps2)