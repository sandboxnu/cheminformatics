#!/usr/bin/env python

import os

smiles = []
bad_smiles_and_reasons = {}

def run_pains_filter(inputs):

  filtered_files = ['bad0.smi', 'bad1.smi', 'bad2.smi', 'bad3.smi']

  if os.path.exists("okmedchem.smi"):
    os.remove("okmedchem.smi")
  if os.path.exists("input.smi"):
    os.remove("input.smi")

  for filtered_file in filtered_files:
    if os.path.exists(filtered_file):
      os.remove(filtered_file)

  with open('input.smi', "w") as file:
      for input_smile in inputs:
          file.write(input_smile)
          file.write("\n")

  os.chdir(os.path.abspath(os.path.dirname(__file__)))
  
  dir_name = os.getcwd()

  dir_name.replace("C:", "/c")

  command = " ".join(["docker", "run", "--rm", "--mount", f"type=bind,destination=/mutable/outside/world,source=\"{dir_name}\"", "--entrypoint", "bash", "ianwatson/lilly_medchem_rules:v1.2", "-c", "\". /etc/profile && rvm use 2.7.1 > /dev/null && cd /mutable/outside/world && /Lilly-Medchem-Rules/Lilly_Medchem_Rules.rb input.smi > okmedchem.smi\""])
  os.system(command)
  bad_files = list(filter(os.path.isfile, filtered_files))

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
    smiles.append(smile.split(" ")[0])

def get_smiles(inputs):
  run_pains_filter(inputs)
  return smiles

def get_bad_smiles(inputs):
  run_pains_filter(inputs)
  return bad_smiles_and_reasons

def get_bad_smiles_and_reasons(inputs):
  run_pains_filter(inputs)
  return [k + v for k, v in bad_smiles_and_reasons.items()]
