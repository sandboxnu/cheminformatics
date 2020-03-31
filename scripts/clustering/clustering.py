import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from app.data.smiles import  convert
from rdkit.ML.Cluster import Butina

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
    print(clusters)
    for clust in clusters:
        index = 0
        for smile in clust:
            similarities = {}
            for othersmile in clust:
                if smile != othersmile:
                    similarity = compare_two_smiles(smiles[smile]['murcko'], smiles[othersmile]['murcko'])
                    similarities[othersmile] = similarity
            smiles[smile]['similarities'] = similarities
            smiles[smile]['isCentroid'] = True if index == 0 else False
            
            if not 'isReclustered' in smiles[smile]:
                smiles[smile]['isReclustered'] =  False # to apply to every smile
            index += 1

    return smiles        

def recluster_singletons(smiles, clusters, recluster_coefficient):
    singletons = []
    realClusters = []
    realSingletons = []

    for clust in clusters: 
        if(len(clust) == 1):
            singletons.append(clust[0])
        else:
            realClusters.append(clust)
    
    for singleton in singletons:
        centroid_similarities = []
        for cluster in realClusters:
            centroid = cluster[0] 
            similarity = compare_two_smiles(singleton, centroid)
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
            print(singleton)
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

def get_row_of_similarities(smile, row, smiles):
    similarities = {}
    for r in row:
        similarities[r] = compare_two_smiles(smiles[smile]['murcko'], smiles[r]['murcko'])
    return similarities
    

