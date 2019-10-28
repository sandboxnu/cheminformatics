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

#smiles with murcko smiles example:
smiles_with_murcko = {'COc1cc(OC)cc(C(=O)NS(=O)(=O)c2ccc(CN3CCN(c4ccccc4)CC3)cc2)c1': {'murcko': 'O=C(NS(=O)(=O)c1ccc(CN2CCN(c3ccccc3)CC2)cc1)c1ccccc1'},
'O=C(NS(=O)(=O)c1ccc(CN2CCN(c3ccccc3)CC2)cc1)c1ccc(CCc2ccccn2)cc1': {'murcko': 'O=C(NS(=O)(=O)c1ccc(CN2CCN(c3ccccc3)CC2)cc1)c1ccc(CCc2ccccn2)cc1'},
'CC(=O)Nc1ccc(S(=O)(=O)NC(=O)c2cn(C)c3ccc(-c4cncnc4)cc23)cc1': {'murcko': 'O=C(NS(=O)(=O)c1ccccc1)c1c[nH]c2ccc(-c3cncnc3)cc12'}}






