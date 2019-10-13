from app import app
from flask import render_template
import scripts.pains.LillyMedchemRules.pains as pains

@app.route('/')
@app.route('/index')
def index():
  pains.run_pains_filter()
  bad_smiles = pains.get_bad_smiles_and_reasons()

  return render_template('index.html', title='Cheminformatics', smiles=bad_smiles)