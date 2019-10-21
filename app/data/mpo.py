from rdkit import Chem
from rdkit.Chem import Descriptors
from app.data.smiles import smiles

def generate_mpo_attributes():
  for smile in smiles.keys():
    molecule = None
    molecule = Chem.MolFromSmiles(smile)
    tpsa = Descriptors.TPSA(molecule)
    logp = Descriptors.MolLogP(molecule)
    mw = Descriptors.MolWt(molecule)
    smiles[smile]['molecule'] = molecule
    smiles[smile]['tpsa'] = tpsa
    smiles[smile]['logp'] = logp
    smiles[smile]['mw'] = mw
  