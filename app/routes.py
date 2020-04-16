from app import app
from flask import render_template, request, session
import scripts.pains.LillyMedchemRules.pains as pains
from app.data.smiles import construct_smiles, filter_smiles, convert, convert_to_smiles
import scripts.clustering.clustering as clustering
from app.data.clustering import get_smiles_json
from app.data.color_functions import color_hex_to_array

@app.route('/')
@app.route('/index')
def index():     
  return render_template('index.html', title='Cheminformatic Analysis')


@app.route("/cluster", methods=['GET', 'POST'])
def upload():
    session.clear()
    if request.method == 'POST':
        try:
            data = request.get_array(field_name='file')
            smiles, include_mpo = construct_smiles(data)

        except Exception as e: 
            return render_template('index.html', title='Cheminformatic Analysis', errors=["Please input a valid file format"])
        
        inputs = smiles.keys()

        #global all_smiles 
        session['all_smiles'] = convert_to_smiles(smiles.copy())
        #global good_smiles
        session['good_smiles'] = convert_to_smiles(filter_smiles(pains.get_smiles(inputs), smiles))
        #global bad_smiles
        bs = convert_to_smiles(pains.get_bad_smiles(inputs)) 
        bad_smiles = bs if isinstance(bs, dict) else {}
        session['bad_smiles'] = bad_smiles
        #global include_mpo
        session['include_mpo'] = include_mpo

        
        #global reasons_for_failure
        reasons_for_failure = []   
        for i in bad_smiles.values():
          if i not in reasons_for_failure:
            reasons_for_failure.append(i) 
        session['reasons_for_failure'] = reasons_for_failure
        session.changed = True
    return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure, include_mpo=session['include_mpo'])

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
  return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure, include_mpo=session['include_mpo'])

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
            
  return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure, include_mpo=session['include_mpo'])  

@app.route('/final_compounds', methods=['GET', 'POST'])
def final_compounds():
  good_smiles = session.get('good_smiles')
  bad_smiles = session.get('bad_smiles')
  reasons_for_failure = session.get('reasons_for_failure')
  if request.method == 'POST':
    try:
      tanimoto = request.form['tanimoto']
      reclusterCoefficient = request.form['reclusterCoefficientValue']
      
      if (float(tanimoto) < 0 or float(tanimoto) > 1):
        return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure, errors=["Please input a valid tanimoto coefficient"], include_mpo=session['include_mpo'])
      elif ((reclusterCoefficient != '') and (float(tanimoto) < 0 or float(tanimoto) > 1)):
        return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure, errors=["Please input a valid recluster coefficient"], include_mpo=session['include_mpo'])
    except:
      return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure, errors=["Please input a valid tanimoto coefficient"], include_mpo=session['include_mpo'])

  color1 = request.form['mpo0Color'] 
  color1_array = color_hex_to_array(color1)

  fp_type = request.form['fp_radio'] 

  if(session['include_mpo']):
    color2 = request.form['mpo6Color'] 
    color2_array = color_hex_to_array(color2)
  else:
    color2 = color1
    color2_array = color1_array

  shouldRecluster = reclusterCoefficient != ''
  cluster = clustering.cluster(list(good_smiles.keys()), fp_type, 1 - float(tanimoto))
  if shouldRecluster : 
    recluster_data = clustering.recluster_singletons(good_smiles, cluster, float(reclusterCoefficient))
    recluster_smiles = recluster_data[0]
    recluster_clusters = recluster_data[1]
    tanimoto_smiles = clustering.get_tanimoto_coeffient_by_cluster(recluster_smiles, recluster_clusters, fp_type)
    get_smiles_json(tanimoto_smiles, float(tanimoto), recluster_clusters, session['include_mpo'], color1_array, color2_array, shouldRecluster)
  else :
    tanimoto_smiles = clustering.get_tanimoto_coeffient_by_cluster(good_smiles, cluster, fp_type)
    get_smiles_json(tanimoto_smiles, float(tanimoto), cluster, session['include_mpo'], color1_array, color2_array, shouldRecluster)

  include_mpo = session['include_mpo']
  session.clear()

  return render_template('cluster.html', title='Cheminformatic Analysis', color1=color1, color2=color2, include_mpo=include_mpo)


@app.after_request
def add_header(response):
  response.cache_control.public = True
  response.cache_control.max_age = 0
  return response

