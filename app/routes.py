from app import app
from flask import render_template
import sys
sys.path.append('/scripts/pains/LillyMedchemRules/')
import scripts.pains.LillyMedchemRules.pains as pains
from app.data.smiles import smiles
from app.data.smiles import smiles_with_murcko
import scripts.clustering.clustering as clustering

@app.route('/')
@app.route('/index')
def index():

  inputs = smiles.keys()
  
  bad_smiles = pains.get_bad_smiles_and_reasons(inputs)
  
  tanimoto_smiles = clustering.add_tanimoto_coefficients(smiles_with_murcko)
  

  return render_template('index.html', title='Cheminformatic Analysis', smiles=bad_smiles)