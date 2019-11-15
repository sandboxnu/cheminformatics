from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import PandasTools

smiles = {}

#smiles with murcko smiles example:
smiles_with_murcko = {'COc1cc(OC)cc(C(=O)NS(=O)(=O)c2ccc(CN3CCN(c4ccccc4)CC3)cc2)c1': {'murcko': 'O=C(NS(=O)(=O)c1ccc(CN2CCN(c3ccccc3)CC2)cc1)c1ccccc1'},
'O=C(NS(=O)(=O)c1ccc(CN2CCN(c3ccccc3)CC2)cc1)c1ccc(CCc2ccccn2)cc1': {'murcko': 'O=C(NS(=O)(=O)c1ccc(CN2CCN(c3ccccc3)CC2)cc1)c1ccc(CCc2ccccn2)cc1'},
'CC(=O)Nc1ccc(S(=O)(=O)NC(=O)c2cn(C)c3ccc(-c4cncnc4)cc23)cc1': {'murcko': 'O=C(NS(=O)(=O)c1ccccc1)c1c[nH]c2ccc(-c3cncnc3)cc12'}}


def construct_smiles(csv):
  if (csv[0] == ['snile', 'label', 'mpo']):
    raise Exception("Malformed file input")

  csv = csv[1:]
  for row in csv:
    smile_string = row[0] 
    smiles[smile_string] = {}
    smiles[smile_string]['murcko'] = convert(smile_string)
    smiles[smile_string]['label'] = row[1]
    smiles[smile_string]['mpo'] = row[2]

def filter_smiles(good_smiles):
  return {smi: data for (smi, data) in smiles.items() if smi in good_smiles}

def convert(mol_smile):
    m = Chem.MolFromSmiles(mol_smile)
    core = MurckoScaffold.GetScaffoldForMol(m)
    return Chem.MolToSmiles(core)

def conver_from_smart(smart):
  conv = Chem.MolFromSmiles(smart)
  back = Chem.MolToSmiles(conv)
  return back


def convert_to_smiles(bad_smiles):
  result = {}
  
  for (smart, reason) in bad_smiles.items():
    smile = conver_from_smart(smart)
    result[smile] = reason

  return result
