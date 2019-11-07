from app import app
from flask import render_template, request, jsonify
import scripts.pains.LillyMedchemRules.pains as pains
from app.data.smiles import smiles, construct_smiles
from app.data.smiles import smiles_with_murcko
import scripts.clustering.clustering as clustering
from app.data.clustering import get_smiles_json

@app.route('/')
@app.route('/index')
def index():     
  return render_template('index.html', title='Cheminformatic Analysis')


@app.route("/cluster", methods=['GET', 'POST'])
def upload():
    if request.method == 'POST':
        data = request.get_array(field_name='file')
        construct_smiles(data)
        print(str(data))
        inputs = smiles.keys()
  
        bad_smiles = pains.get_bad_smiles_and_reasons(inputs)

        tanimoto_smiles = clustering.add_tanimoto_coefficients(smiles)
  
        tanimoto_smile_strings = [" ".join([str(i) for i in v['similarities'].values()]) for v in tanimoto_smiles.values()]

        get_smiles_json(tanimoto_smiles)
        return render_template('cluster.html', title='Cheminformatic Analysis')
