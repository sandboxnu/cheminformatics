##example smile
smile = dict([('smile', ''), ('murcko-smile', '')])

##example lists of smiles
smile_strings = ["COC1=CC2=C(C=C1OC)C(C1=CC=CC=C1)N(CC(=O)NS(=O)(=O)C1=CC=C(CN(C)CCC3=CC=CC=C3)C=C1)CC2",
  "COC1=CC(=CC(=C1)C(=O)NS(=O)(=O)C1=CC=C(CN(C)CCC2=CC=CC=C2)C=C1)OC",
  "CN(CCC1=CC=CC=C1)CC1=CC=C(S(=O)(=O)NC(=O)CCCN2CC(=C(C3=CC=C(F)C=C3F)C2=O)C2=CC=C(F)C=C2F)C=C1",
  "CCOC1=CC=C(C2=CC=C3C(=C2)C(=CN3C)C(=O)NS(=O)(=O)C2=CC=C(CN(C)CCC3=CC=CC=C3)C=C2)C=C1F",
  "CN(CCC1=CC=CC=C1)CC1=CC=C(S(=O)(=O)NC(=O)C2=CC=C(C3=CC=CC=N3)C=C2)C=C1",
  "CN(CCC1=CC=CC=C1)CC1=CC=C(S(=O)(=O)NC(=O)C2=CC=C3C(=C2)C=CN3CC2=CC=C(C(F)(F)F)C=C2C(F)(F)F)C=C1",
  "CN(CCC1=CC=CC=C1)CC1=CC=C(S(=O)(=O)NC(=O)C2=CC=C3OC(=CC3=C2)C2=CC=CC=C2)C=C1",
  "CN(CCC1=CC=CC=C1)CC1=CC=C(S(=O)(=O)NC(=O)C2=CC=C(CN3C(=O)C(=C(C4=CC=C(F)C=C4F)C3=O)C3=CC=C(F)C=C3F)C=C2)C=C1",
  "CN(CCC1=CC=CC=C1)CC1=CC=C(S(=O)(=O)NC(=O)C2=CC=CC(=C2)NCCC2=CC=CC=C2)C=C1",
  "CN(CCC1=CC=CC=C1)CC1=CC=C(S(=O)(=O)NC(=O)C2=CC=C3C(=C2)C(=CN3C)C(=O)C2=CC=CS2)C=C1"]

smiles = dict.fromkeys(smile_strings, {})

#add keys to smiles in the list
for smile in smile_strings: 
    smiles[smile]['murcko'] = ""






