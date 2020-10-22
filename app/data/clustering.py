import json
import os
import scripts.clustering.clustering as clustering


def get_smiles_json(smiles, cutoff, clusters, include_property, prop_color1=[255,0,0], prop_color2=[0, 255, 0], shouldRecluster=False):
  nodes = []
  edges = []

  for smile_name, smile_data in smiles.items():
    prop_name = smile_data.get('property_name', "")
    prop_val = smile_data['property']
    if include_property:
      label = smile_data['label'] + "\n" + prop_name + ' ' + str(prop_val)
    else:
      label = smile_data['label']

    smile_edges = []
    smile_node = {}
    smile_node['murcko'] = smile_data['murcko']
    smile_node['data'] = {}
    
    if shouldRecluster:
      smile_node['data'] = {'id': smile_name, 'label': label,
        'mpo': prop_val, 'reclustered': smile_data['isReclustered'],'centroid': smile_data['isCentroid'] , 'type': 'node'}
    else:
      smile_node['data'] = {'id': smile_name, 'label': label,
        'mpo': prop_val,'centroid': smile_data['isCentroid'] , 'type': 'node'}
  
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

  result = {'nodes': nodes, 'edges': edges, 'clusterInfo': clusters, 'color1': prop_color1, 'color2': prop_color2}

  os.chdir(os.path.abspath(os.path.dirname(__file__)))
  
  with open('../static/js/data.json', 'w+') as outfile:
    json.dump(result, outfile, indent=4, sort_keys=True)

  return result



