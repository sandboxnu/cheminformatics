from app import app
from flask import render_template
import scripts.pains.LillyMedchemRules.pains as pains
from app.data.smiles import smiles
import app.data.mpo as mpo

@app.route('/')
@app.route('/index')
def index():

  mpo.generate_mpo_attributes()

  inputs = smiles.keys()
  
  bad_smiles = pains.get_bad_smiles_and_reasons(inputs)

  return render_template('index.html', title='Cheminformatic Analysis', smiles=bad_smiles)