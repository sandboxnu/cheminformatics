from smiles import smiles 


def get_smiles_json():
  nodes = {}
  edges = {}

  for smile_string, smile_data in smiles.values():
    smile_node = {}
    smile_node['data'] = {'id': smile_data['murcko']}
