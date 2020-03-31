import json
import os
import scripts.clustering.clustering as clustering


def get_smiles_json(smiles, cutoff, clusters, include_mpo, mpo_color1=[255,0,0], mpo_color2=[0, 255, 0]):
  nodes = []
  edges = []

  for smile_name, smile_data in smiles.items():

    label = smile_data['label']
    
    if include_mpo:
      label = smile_data['label'] + '\nmpo: ' + str(smile_data['mpo'])

    smile_edges = []
    smile_node = {}
    smile_node['murcko'] = smile_data['murcko']
    smile_node['data'] = {}
    smile_node['data'] = {'id': smile_name, 'label': label,
      'mpo': smile_data['mpo'], 'centroid': smile_data['isCentroid'] , 'type': 'node'}    
    
    
    nodes.append(smile_node)
    for sim, similarity_coefficient in smile_data['similarities'].items():

      if similarity_coefficient >= cutoff and similarity_coefficient != 1:
        new_edge = {}
        new_edge['data'] = {}
        new_edge['data']['source'] = smile_name
        new_edge['data']['target'] = sim 
        new_edge['data']['weight'] = similarity_coefficient
        new_edge['data']['type'] = "controls-state-change-of"
        
        smile_edges.append(new_edge)
    
    edges+= smile_edges 

  result = {'nodes': nodes, 'edges': edges, 'clusterInfo': clusters, 'color1': mpo_color1, 'color2': mpo_color2}

  os.chdir(os.path.abspath(os.path.dirname(__file__)))
  
  with open('../static/js/data.json', 'w+') as outfile:
    json.dump(result, outfile, indent=4, sort_keys=True)



