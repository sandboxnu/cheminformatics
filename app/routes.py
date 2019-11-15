from app import app
from flask import render_template, request
import scripts.pains.LillyMedchemRules.pains as pains
from app.data.smiles import smiles, construct_smiles, filter_smiles, convert, convert_to_smiles
import scripts.clustering.clustering as clustering
from app.data.clustering import get_smiles_json

global bad_smiles 
global good_smiles
global all_smiles

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
        global all_smiles 
        all_smiles = smiles.copy()
        global good_smiles
        good_smiles = convert_to_smiles(filter_smiles(pains.get_smiles(inputs)))
        global bad_smiles
        bad_smiles = convert_to_smiles(pains.get_bad_smiles(inputs))
        
    return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles)

@app.route('/verify_pains', methods=['GET', 'POST'])
def verify_pains():
  global good_smiles
  global bad_smiles
  global all_smiles
  if request.method == 'POST':
    form = request.form
    if form['action'] == 'Drop Selected Compounds':
      for smile in form:
        if(smile != 'action'):
          del bad_smiles[smile]
    elif form['action'] == 'Ignore Errors':
      
      for smile in form:
        if(smile != 'action'):
          del bad_smiles[smile]
          good_smiles[smile] = all_smiles[smile]

  return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles)

@app.route('/final_compounds', methods=['GET', 'POST'])
def final_compounds():
  global good_smiles
  if request.method == 'POST':
    try:
      tanimoto = request.form['tanimoto']
      
      if (float(tanimoto) < 0 or float(tanimoto) > 1):
        return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, errors=["Please input a valid tanimoto coefficient"])
    except:
      return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, errors=["Please input a valid tanimoto coefficient"])
  
  tanimoto_smiles = clustering.add_tanimoto_coefficients(good_smiles)

  get_smiles_json(tanimoto_smiles, float(tanimoto))

  return render_template('cluster.html', title='Cheminformatic Analysis')
