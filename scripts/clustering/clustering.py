import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem
from app.data.smiles import  convert
from rdkit.ML.Cluster import Butina
from rdkit.Chem.AtomPairs import Pairs

radius = 4
def compare_two_smiles(smile1, smile2, fp_type):
    smile1Ms = Chem.MolFromSmiles(smile1)
    smile2Ms = Chem.MolFromSmiles(smile2)

    if fp_type == "atom-pair":
        fps1 = Pairs.GetAtomPairFingerprintAsBitVect(smile1Ms)
        fps2 = Pairs.GetAtomPairFingerprintAsBitVect(smile2Ms)
    elif fp_type == "fpm":
        fps1 = FingerprintMols.FingerprintMol(smile1Ms)
        fps2 = FingerprintMols.FingerprintMol(smile2Ms)
    else:
        fps1 = AllChem.GetMorganFingerprintAsBitVect(smile1Ms, radius, nBits=1024)
        fps2 = AllChem.GetMorganFingerprintAsBitVect(smile2Ms, radius, nBits=1024)

    return DataStructs.TanimotoSimilarity(fps1, fps2)
    
#return cluster of smile_keys
def cluster(smile_keys, fp_type, cutoff=0.15):
    #note: it seems cutoff is one - similarity coefficient, it's euclidean distance I think??
    nfps = len(smile_keys)
    dists = []
    
    data = [None] * nfps
    for i in range(0, nfps):
        murcko = convert(smile_keys[i])
        mol = Chem.MolFromSmiles(murcko)
        if fp_type == "atom-pair":
            fps = Pairs.GetAtomPairFingerprintAsBitVect(mol)
        elif fp_type == "fpm":
            fps = FingerprintMols.FingerprintMol(mol)
        else:
            fps = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=1024)
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

def get_tanimoto_coeffient_by_cluster(smiles, clusters, fp_type):
    for clust in clusters:
        index = 0
        for smile in clust:
            similarities = {}
            for othersmile in clust:
                if smile != othersmile:
                    similarity = compare_two_smiles(smiles[smile]['murcko'], smiles[othersmile]['murcko'], fp_type)
                    similarities[othersmile] = similarity
            smiles[smile]['similarities'] = similarities
            smiles[smile]['isCentroid'] = True if index == 0 else False
            if not 'isReclustered' in smiles[smile]:
                smiles[smile]['isReclustered'] =  False # to apply to every smile
            index += 1
    return smiles        

def recluster_singletons(smiles, clusters, recluster_coefficient, fp_type):
    singletons = []
    realClusters = []
    realSingletons = []
    for clust in clusters: 
        if(len(clust) == 1):
            singletons.append(clust[0])
        else:
            realClusters.append(clust)

    if len(singletons) == len(clusters):
        return [smiles, clusters]

    for singleton in singletons:
        centroid_similarities = []
        for cluster in realClusters:
            centroid = cluster[0] 
            similarity = compare_two_smiles(singleton, centroid, fp_type)
            centroid_similarities.append(similarity)
        
        max_sim_index = max(centroid_similarities)
        if max_sim_index < 0:
            continue
        max_sim = centroid_similarities[max_sim_index]
        if(max_sim >= recluster_coefficient):
            smiles[singleton]['isReclustered'] = True
            smiles[singleton]['isCentroid'] = False
            realClusters[max_sim_index].append(singleton)
        else:
            realSingletons.append(singleton)
    
    for singleton in realSingletons:
        realClusters.append([singleton])
            
    return [smiles, realClusters]

def max(array):
    currMax = -1
    index = -1
    for x in array:
        index = index + 1
        if x > currMax:
            currMax = x
    return index        

def get_row_of_similarities(smile, row, smiles, fp_type):
    similarities = {}
    for r in row:
        similarities[r] = compare_two_smiles(smiles[smile]['murcko'], smiles[r]['murcko'], fp_type)
    return similarities
    

