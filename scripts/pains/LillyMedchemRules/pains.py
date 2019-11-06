#!/usr/bin/env python

import os, sys


smiles = []
bad_smiles_and_reasons = {}

def run_pains_filter(inputs):

  os.chdir(os.path.abspath(os.path.dirname(__file__)))

  pains = os.system('ruby Lilly_Medchem_Rules.rb input.smi > okmedchem.smi')

  bad_files = list(filter(os.path.isfile, ['bad0.smi', 'bad1.smi', 'bad2.smi', 'bad3.smi']))

  with open('input.smi', "w") as file:
      for input_smile in inputs:
          file.write(input_smile)
          file.write("\n")


  for bad_file in bad_files:
    with open(bad_file) as file:
      bad_smiles = [line.rstrip("\n") for line in file]

      for smile_info in bad_smiles:
        info = smile_info.split(" ")
        smile = info[0]
        reason = " ".join(info[1:])
        bad_smiles_and_reasons[smile] = reason


  good_smiles = [line.rstrip("\n") for line in open('okmedchem.smi')]
  for smile in good_smiles:
    smiles.append(smile)

def get_smiles(inputs):
  run_pains_filter(inputs)
  return smiles

def get_bad_smiles(inputs):
  run_pains_filter(inputs)
  return bad_smiles_and_reasons

def get_bad_smiles_and_reasons(inputs):
  run_pains_filter(inputs)
  return [k + v for k, v in bad_smiles_and_reasons.items()]
