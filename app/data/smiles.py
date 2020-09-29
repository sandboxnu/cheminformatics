from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import PandasTools

def construct_smiles(csv):

  formattedTitleRow = []

  for title in csv[0]:
    formattedTitleRow.append(title.lower())

  print(formattedTitleRow)
  if (formattedTitleRow != ['smiles', 'label', 'mpo'] and formattedTitleRow != ['smiles', 'label']
  and formattedTitleRow != ['\ufeffsmiles', 'label', 'mpo'] and formattedTitleRow != ['\ufeffsmiles', 'label']):
    raise Exception("Malformed file input")
  
  if(len(csv[0]) == 3):
    include_property = csv[0][2]
  else:
    include_property = None

  smiles = {}
  
  csv = csv[1:]
  
  for row in csv:
    smile_string = row[0] 
    smiles[smile_string] = {}
    smiles[smile_string]['murcko'] = convert(smile_string)
    smiles[smile_string]['label'] = row[1]
    if include_property is not None:
      smiles[smile_string]['property'] = row[2]
      smiles[smile_string]['property_name'] = include_property
    else:
      smiles[smile_string]['property'] = 0

  return smiles, include_property is not None

def filter_smiles(good_smiles, smiles):
  return {smi: data for (smi, data) in smiles.items() if smi in good_smiles}

def convert(mol_smile):
    m = Chem.MolFromSmiles(mol_smile)
    core = MurckoScaffold.GetScaffoldForMol(m)
    return Chem.MolToSmiles(core)

def convert_from_smart(smart):
  conv = Chem.MolFromSmiles(smart)
  back = Chem.MolToSmiles(conv)
  return back


def convert_to_smiles(bad_smiles):
  result = {}
  
  for (smart, reason) in bad_smiles.items():
    smile = convert_from_smart(smart)
    result[smile] = reason

  return result



