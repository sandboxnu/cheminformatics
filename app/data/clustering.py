import json
import os


def get_smiles_json(smiles, cutoff):
  nodes = []
  edges = []

  for smile_name, smile_data in smiles.items():
    smile_edges = []

    smile_node = {}
    smile_node['murcko'] = smile_data['murcko']
  
    smile_node['data'] = {}
    smile_node['data'] = {'id': smile_name, 'label': smile_data['label'] + '\nmpo: ' + str(smile_data['mpo'])}
  
    nodes.append(smile_node)
    for sim, similarity_coefficient in smile_data['similarities'].items():

      # TODO: use actual similarity cutoff
      if similarity_coefficient >= cutoff and similarity_coefficient != 1:
        new_edge = {}
        new_edge['data'] = {}
        new_edge['data']['source'] = smile_name
        new_edge['data']['target'] = sim 
        new_edge['data']['weight'] = similarity_coefficient
        new_edge['data']['type'] = "controls-state-change-of"
        smile_edges.append(new_edge)
    
    edges+= smile_edges

  result = {'nodes': nodes, 'edges': edges}

  os.chdir(os.path.abspath(os.path.dirname(__file__)))
  
  with open('../static/js/data.json', 'w+') as outfile:
    json.dump(result, outfile, indent=4, sort_keys=True)



