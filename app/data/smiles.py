from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import PandasTools

def construct_smiles(csv):

  formattedTitleRow = []

  for title in csv[0]:
    formattedTitleRow.append(title.lower())

  if (formattedTitleRow[:2] != ['\ufeffsmiles', 'label'] and formattedTitleRow[:2] != ['smiles', 'label']):
    raise Exception("Malformed file input")

  if(len(csv[0]) == 3):
    include_property = csv[0][2]
  else:
    include_property = None

  smiles = {}
  csv = csv[1:]

  highest_val = float('-inf')
  lowest_val = float('inf')

  for row in csv:
    smile_string = row[0]
    smiles[smile_string] = {}
    smiles[smile_string]['murcko'] = get_murcko_smile(smile_string)
    smiles[smile_string]['smart'] = convert_from_smart(smile_string)
    smiles[smile_string]['label'] = row[1]
    if include_property is not None:
      prop_val = float(row[2])
      highest_val = max(prop_val, highest_val)
      lowest_val = min(prop_val, lowest_val)
      smiles[smile_string]['property'] = prop_val
      smiles[smile_string]['property_name'] = include_property
    else:
      smiles[smile_string]['property'] = 0

  return smiles, include_property, highest_val, lowest_val


def filter_smiles(good_smiles, smiles):
  return {smi: data for (smi, data) in smiles.items() if data['smart'] in convert_array_of_smarts_to_smiles(good_smiles)}

def get_murcko_smile(mol_smile):
    m = Chem.MolFromSmiles(mol_smile)
    core = MurckoScaffold.GetScaffoldForMol(m)
    return Chem.MolToSmiles(core)

def convert_from_smart(smart):
  conv = Chem.MolFromSmiles(smart)
  back = Chem.MolToSmiles(conv)
  return back



def sanitize(smile):
  try:
    conv = Chem.MolFromSmiles(smile)

    Chem.SanitizeMol(conv)
    return Chem.MolToSmiles(conv)
  except Exception as e:
    return False

def convert_array_of_smarts_to_smiles(smarts):
  return [convert_from_smart(smart) for smart in smarts]

def convert_to_smiles(bad_smiles):
  result = {}

  for (smart, reason) in bad_smiles.items():
    smile = convert_from_smart(smart)
    result[smile] = reason

  return result

def convert_to_smiles_and_labels(bad_smiles, all_smiles):
  result = {}
  label_failures = []

  for (smart, reason) in bad_smiles.items():
    smile = convert_from_smart(smart)
    label = ''
    try:
      label = all_smiles[smile]['label']
    except:
      print('could not find label for ' + smile)
      label_failures.append(smile)
    result[smile] = {'label': label, 'reason': reason}

  return result, label_failures
