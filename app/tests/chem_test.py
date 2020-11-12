from app.data.smiles import construct_smiles
from app.data.clustering import get_smiles_json
import scripts.clustering.clustering as clustering
import json
import csv
from os import path

def read_example_data(file_name):
  base = path.dirname(__file__)
  file_path = path.abspath(path.join(base, "..", "..", 'example_data', f'{file_name}.csv'))
  with open(file_path, newline='') as example_file:
    compounds = list(csv.reader(example_file))
  return compounds

def read_smiles_json():
  base = path.dirname(__file__)
  file_path = path.abspath(path.join(base, "..", 'static', 'js', 'data.json'))
  with open(file_path, newline='') as json_file:
    data = json.load(json_file)
  return data

def test_construct_smiles():
  """
  test that construct smiles returns expected results
  """
  compounds = read_example_data("split_data")
  smiles, include_property, highest_val, lowest_val = construct_smiles(compounds)

  assert(include_property == "mpo")
  assert(len(smiles.values()) == 5)
  assert(smiles['c1ccccc1OCC(=O)N/N=C\\2NC(=N\\[H])\\c(c23)cccc3']['property'] == 3.142)
  assert(smiles['c1ccccc1OCC(=O)N/N=C\\2NC(=N\\[H])\\c(c23)cccc3']['property_name'] == "mpo")
  assert(highest_val == 3.761)
  assert(lowest_val == 1.85)

def test_cluster():
  compounds = read_example_data("split_data")
  smiles = construct_smiles(compounds)[0]
  clusters = clustering.cluster(list(smiles.keys()), "atom-pair")
  assert (len(clusters) == 5)

  clusters = clustering.cluster(list(smiles.keys()), "atom-pair", 0.6)
  assert (len(clusters) == 2)

  clusters = clustering.cluster(list(smiles.keys()), "maccs", 0.6)
  assert (len(clusters) == 1)
