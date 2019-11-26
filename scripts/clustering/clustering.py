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
    
#return cluster of smile_keys
def cluster(smile_keys, cutoff=0.15):
    #note: it seems cutoff is one - similarity coefficient, it's euclidean distance I think??
    nfps = len(smile_keys)
    dists = []
    
    data = [None] * nfps
    for i in range(0, nfps):
        murcko = convert(smile_keys[i])
        mols = Chem.MolFromSmiles(murcko)
        fps = FingerprintMols.FingerprintMol(mols)
        data[i] = fps
    

    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(data[i],data[:i])
        dists.extend([1-x for x in sims])

    result = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    clusters = []
    for tuple in result:
        set = []
        #elements are index of smile in original array
        for element in tuple:
            corresponding_smile = smile_keys[element]
            set.append(corresponding_smile)
        clusters.append(set)
    return clusters   

def in_same_cluster(s1, s2, clusters):
    for clust in clusters:
        if s1 in clust and s2 in clust:
            return True
        
    return False            

def get_tanimoto_coeffient_by_cluster(smiles, clusters):
    for clust in clusters:
        for smile in clust:
            similarities = {}
            for othersmile in clust:
                if smile != othersmile:
                    similarity = compare_two_smiles(smiles[smile]['murcko'], smiles[othersmile]['murcko'])
                    similarities[othersmile] = similarity
            smiles[smile]['similarities'] = similarities
    return smiles        