from app import app
from flask import render_template, request, session
import scripts.pains.LillyMedchemRules.pains as pains
from app.data.smiles import smiles, construct_smiles, filter_smiles, convert, convert_to_smiles
import scripts.clustering.clustering as clustering
from app.data.clustering import get_smiles_json
from app.data.color_functions import color_hex_to_array

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
        session.clear()
        inputs = smiles.keys()
        #global all_smiles 
        session['all_smiles'] = convert_to_smiles(smiles.copy())
        #global good_smiles
        session['good_smiles'] = convert_to_smiles(filter_smiles(pains.get_smiles(inputs)))
        #global bad_smiles
        bs = convert_to_smiles(pains.get_bad_smiles(inputs)) 
        bad_smiles = bs if isinstance(bs, dict) else {}
        session['bad_smiles'] = bad_smiles
        
        #global reasons_for_failure
        reasons_for_failure = []   
        for i in bad_smiles.values():
          if i not in reasons_for_failure:
            reasons_for_failure.append(i) 
        session['reasons_for_failure'] = reasons_for_failure
        
    return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure)

@app.route('/verify_pains', methods=['GET', 'POST'])
def verify_pains():
  good_smiles = session.get('good_smiles')
  bad_smiles = session.get('bad_smiles')
  all_smiles = session.get('all_smiles')
  reasons_for_failure = session.get('reasons_for_failure')
  #TODO: check if reasons are still even the list? Does it matter?
  if request.method == 'POST':
    form = request.form
    if form['action'] == 'Drop Selected Compounds':
      for smile in form:
        if(smile != 'action'):
          del bad_smiles[smile]
          session['bad_smiles'] = bad_smiles
          session.changed = True

    elif form['action'] == 'Ignore Errors':
      
      for smile in form:
        if(smile != 'action'):
          del bad_smiles[smile]
          session['bad_smiles'] = bad_smiles
          good_smiles[smile] = all_smiles[smile]
          session['good_smiles'][smile] = all_smiles[smile]
          session.changed = True
    
  return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure)

@app.route('/verify_pains_by_error', methods=['GET', 'POST'])
def verify_pains_by_error():
  good_smiles = session.get('good_smiles')
  bad_smiles = session.get('bad_smiles')
  all_smiles = session.get('all_smiles')
  reasons_for_failure = session.get('reasons_for_failure')

  if request.method == 'POST':
    form = request.form
    if form['action'] == 'Drop Selected Errors':
      for reason in form:
        if(reason in reasons_for_failure):
          reasons_for_failure.remove(reason)
          session['reasons_for_failure'] = reasons_for_failure
          smiles_to_remove = []
          for (smile, smile_reason) in bad_smiles.items():
            if(reason == smile_reason):
              smiles_to_remove.append(smile)
          for smile in smiles_to_remove:
            del bad_smiles[smile]
            session['bad_smiles'] = bad_smiles
            session.changed = True    
         
    elif form['action'] == 'Ignore Selected Errors':
      for reason in form:
        if(reason in reasons_for_failure):
          reasons_for_failure.remove(reason)
          session['reasons_for_failure'] = reasons_for_failure
          smiles_to_remove = []
          for (smile, smile_reason) in bad_smiles.items():
            if(reason == smile_reason):
              smiles_to_remove.append(smile)
          for smile in smiles_to_remove:
            del bad_smiles[smile]    
            good_smiles[smile] = all_smiles[smile]
            session['good_smiles'][smile] = all_smiles[smile]
            session.changed = True
            
  return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure)  

@app.route('/final_compounds', methods=['GET', 'POST'])
def final_compounds():
  good_smiles = session.get('good_smiles')
  bad_smiles = session.get('bad_smiles')
  if request.method == 'POST':
    try:
      tanimoto = request.form['tanimoto']
      
      if (float(tanimoto) < 0 or float(tanimoto) > 1):
        return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, errors=["Please input a valid tanimoto coefficient"])
    except:
      return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, errors=["Please input a valid tanimoto coefficient"])
  
  #TODO: Allow user to input these colors. always two colors. one for 0. one for 1.
  #TODO: Add legend for colors in the front end(It is very difficult)
  color1 = request.form['mpo0Color'] 
  color2 = request.form['mpo6Color'] 
  color1_array = color_hex_to_array(color1)
  color2_array = color_hex_to_array(color2)


  cluster = clustering.cluster(list(good_smiles.keys()), 1 - float(tanimoto))
  tanimoto_smiles = clustering.get_tanimoto_coeffient_by_cluster(good_smiles, cluster)

  get_smiles_json(tanimoto_smiles, float(tanimoto), cluster, color1_array, color2_array)
  session.clear()

  return render_template('cluster.html', title='Cheminformatic Analysis', color1=color1, color2=color2)


@app.after_request
def add_header(response):
  response.cache_control.public = True
  response.cache_control.max_age = 0
  return response

