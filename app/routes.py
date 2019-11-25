from app import app
from flask import render_template, request
import scripts.pains.LillyMedchemRules.pains as pains
from app.data.smiles import smiles, construct_smiles, filter_smiles, convert, convert_to_smiles
import scripts.clustering.clustering as clustering
from app.data.clustering import get_smiles_json
from app.data.color_functions import color_hex_to_array

global bad_smiles 
global good_smiles
global all_smiles
global reasons_for_failure

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
        all_smiles = convert_to_smiles(smiles.copy())
        global good_smiles
        good_smiles = convert_to_smiles(filter_smiles(pains.get_smiles(inputs)))
        global bad_smiles
        bad_smiles = convert_to_smiles(pains.get_bad_smiles(inputs))

        global reasons_for_failure
        reasons_for_failure = []   
        for i in bad_smiles.values():
          if i not in reasons_for_failure:
            reasons_for_failure.append(i) 
        
    return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure)

@app.route('/verify_pains', methods=['GET', 'POST'])
def verify_pains():
  global good_smiles
  global bad_smiles
  global all_smiles
  global reasons_for_failure
  #TODO: check if reasons are still even the list? Does it matter?
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

  return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure)

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
  #TODO: Allow user to input these colors
  #TODO: Add legend for colors in the front end(It is very difficult)
  color1 = '#135476'
  color2 = '#ff0000'
  color1_array = color_hex_to_array(color1)
  color2_array = color_hex_to_array(color2)


  cluster = clustering.cluster(list(good_smiles.keys()), 1 - float(tanimoto))

  get_smiles_json(tanimoto_smiles, float(tanimoto), cluster, color1_array, color2_array)

  return render_template('cluster.html', title='Cheminformatic Analysis', color1=color1, color2=color2)


@app.after_request
def add_header(response):
  response.cache_control.public = True
  response.cache_control.max_age = 0
  return response