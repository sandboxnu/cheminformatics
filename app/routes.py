from app import app
from flask import render_template, request
import scripts.pains.LillyMedchemRules.pains as pains
from app.data.smiles import smiles, construct_smiles, filter_smiles
import scripts.clustering.clustering as clustering
from app.data.clustering import get_smiles_json

@app.route('/')
@app.route('/index')
def index():     
  return render_template('index.html', title='Cheminformatic Analysis')


@app.route("/cluster", methods=['GET', 'POST'])
def upload():
    if request.method == 'POST':
        try:
            data = request.get_array(field_name='file')
            construct_smiles(data)
        except: 
            return render_template('index.html', title='Cheminformatic Analysis', errors=["Please input a valid file format"])

        inputs = smiles.keys()
        
        good_smiles = filter_smiles(pains.get_smiles(inputs))

        tanimoto_smiles = clustering.add_tanimoto_coefficients(good_smiles)

        get_smiles_json(tanimoto_smiles)

        return render_template('cluster.html', title='Cheminformatic Analysis')
