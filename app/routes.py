from app import app
from flask import render_template, request
import scripts.pains.LillyMedchemRules.pains as pains
from app.data.smiles import smiles, construct_smiles, filter_smiles
import scripts.clustering.clustering as clustering
from app.data.clustering import get_smiles_json

global bad_smiles 
global good_smiles

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
        global bad_smiles
        bad_smiles = pains.get_bad_smiles(inputs)
        
    return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles)

@app.route('/')
@app.route('/pains_verify_and_coefficient_use')
def pains_verify_and_coefficient_use():
  inputs = smiles.keys()
  global bad_smiles #use the global version 
  
  bad_smiles = pains.get_bad_smiles(inputs)
  
  print(bad_smiles)
  return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles)


@app.route('/verify_pains', methods=['GET', 'POST'])
def verify_pains():
  global good_smiles
  global bad_smiles
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
          good_smiles[smile] = {}

  return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles)

@app.route('/final_compounds', methods=['GET', 'POST'])
def final_compounds():
  if request.method == 'POST':
    tanimoto = request.form['tanimoto']
  
  #get real final smiles
  inputs = smiles.keys()
  global bad_smiles
  bad_smiles={}

  tanimoto_smiles = clustering.add_tanimoto_coefficients(smiles)
  
  tanimoto_smile_strings = [" ".join([str(i) for i in v['similarities'].values()]) for v in tanimoto_smiles.values()]

  get_smiles_json(tanimoto_smiles, .2)

  return render_template('cluster.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, tanimoto_smiles=tanimoto_smile_strings)
