import chem from rdkit
#dictionary of dictionaries, each key is the smile

#add a key with the tanimoto coefficient
def add_tanimoto_coefficients(smiles):
    #figure out how to get rdkit in here
    for smile in smiles:
        smiles[smile]['tanimoto'] = 0 #some number here

    return smiles