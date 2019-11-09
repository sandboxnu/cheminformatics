from app import app
from flask import render_template
import scripts.pains.LillyMedchemRules.pains as pains
from app.data.smiles import smiles
from app.data.smiles import smiles_with_murcko
import scripts.clustering.clustering as clustering
from app.data.clustering import get_smiles_json

global bad_smiles

@app.route('/')
@app.route('/index')
def index():
  return render_template('index.html', title='Cheminformatic Analysis')

@app.route('/')
@app.route('/cluster')
def cluster():

  inputs = smiles.keys()
  
  bad_smiles = pains.get_bad_smiles_and_reasons(inputs)

  tanimoto_smiles = clustering.add_tanimoto_coefficients(smiles)
  
  tanimoto_smile_strings = [" ".join([str(i) for i in v['similarities'].values()]) for v in tanimoto_smiles.values()]

  get_smiles_json(tanimoto_smiles)


  return render_template('cluster.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, tanimoto_smiles=tanimoto_smile_strings)

@app.route('/')
@app.route('/pains_verify_and_coefficient_use')
def pains_verify_and_coefficient_use():
  
  inputs = smiles.keys()
  global bad_smiles #use the global version 
  
  bad_smiles = pains.get_bad_smiles(inputs)
  
  print(bad_smiles)
  return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles)

@app.route('/')
@app.route('/verify_pains')
def verify_pains():
  global bad_smiles

  return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles)
