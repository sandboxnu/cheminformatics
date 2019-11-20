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
    

def cluster(smile_keys, cutoff=0.50):
    #note: it seems cutoff is one - similarity coefficient
    nfps = len(smile_keys)
    dists = []
    #turn into murckos??
    data = [None] * nfps
    for i in range(0, nfps):
        murcko = convert(smile_keys[i])
        mols = Chem.MolFromSmiles(murcko)
        fps = FingerprintMols.FingerprintMol(mols)
        data[i] = fps
    

    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(data[i],data[:i])
        dists.extend([1-x for x in sims])

    print(dists)

    result = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    print(result)




def compare_murcko_smiles(smile1, smile2):
    murcko1 = convert(smile1)
    murcko2 = convert(smile2)
    smile1Ms = Chem.MolFromSmiles(murcko1)
    smile2Ms = Chem.MolFromSmiles( murcko2)
    fps1 = FingerprintMols.FingerprintMol(smile1Ms)
    fps2 = FingerprintMols.FingerprintMol(smile2Ms)

    return DataStructs.FingerprintSimilarity(fps1, fps2)