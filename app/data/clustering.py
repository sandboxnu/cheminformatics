import json
import os
import scripts.clustering.clustering as clustering
import pandas as pd


def get_smiles_json(smiles, cutoff, clusters, include_property, lowest_val, highest_val, prop_color1=[255,0,0], prop_color2=[0, 255, 0], shouldRecluster=False):
  nodes = []
  edges = []

  for smile_name, smile_data in smiles.items():
    prop_name = smile_data.get('property_name', "mpo")
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
        'prop_name': prop_name, 'prop_val': prop_val, 'reclustered': smile_data['isReclustered'],'centroid': smile_data['isCentroid'] , 'type': 'node'}
    else:
      smile_node['data'] = {'id': smile_name, 'label': label,
        'prop_name': prop_name, 'prop_val': prop_val, 'centroid': smile_data['isCentroid'] , 'type': 'node'}
  
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

  result = {'nodes': nodes, 'edges': edges, 'clusterInfo': clusters, 'color1': prop_color1, 'color2': prop_color2, 'lowest_val': lowest_val, 'highest_val': highest_val}

  if include_property is None:
    result['lowest_val'] = 0
    result['highest_val'] = 1

  os.chdir(os.path.abspath(os.path.dirname(__file__)))
  
  with open('../static/js/data.json', 'w+') as outfile:
    json.dump(result, outfile, indent=4, sort_keys=True)

  return result

def create_dataframe(result, include_property):
  nodes_murcko = []
  nodes_type = []
  nodes_id = []
  nodes_label = []
  nodes_property = []
  nodes_reclustered = []

  for node in result['nodes']:
    nodes_murcko.append(node['murcko'])
    nodes_id.append(node['data']['id'])
    
    if node['data']['centroid']:
      nodes_type.append('centroid')
    else:
      nodes_type.append('non centroid')

    nodes_reclustered.append(node['data'].get('reclustered', ''))

    if include_property:
      nodes_property.append(node['data']['prop_val'])
    
    nodes_label.append(node['data']['label'].split('\n')[0])

  nodes_cluster = [0] * len(nodes_id)
  
  for group, cluster in enumerate(result['clusterInfo']):
    for index, id in enumerate(nodes_id):
      if id in cluster:
        nodes_cluster[index] = group


  df = pd.DataFrame({'id': nodes_id, 'murcko': nodes_murcko, 'label': nodes_label, 'type': nodes_type, 'cluster': nodes_cluster, 'reclustered':nodes_reclustered})
  
  if include_property:
    df.insert(5, include_property, nodes_property)
  
  return df

