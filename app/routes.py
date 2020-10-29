from flask.helpers import send_file
from app import app
from flask import render_template, request, session
import scripts.pains.LillyMedchemRules.pains as pains
from app.data.smiles import construct_smiles, filter_smiles, convert_to_smiles
import scripts.clustering.clustering as clustering
from app.data.clustering import get_smiles_json
from app.data.color_functions import color_hex_to_array
import pandas as pd

@app.route('/')
@app.route('/welcome')
def welcome():
  return render_template('welcome.html', title='Welcome to Cheminformatics Analysis')


@app.route('/index')
@app.route('/index', methods=['GET', 'POST'])
def index():     
  session.clear()
  if request.method == 'POST':
    try:
      data = request.get_array(field_name='file')
      smiles, include_property = construct_smiles(data)
      session['smiles'] = smiles
      session['include_property'] = include_property

      if include_property:
        unique_compounds = pd.DataFrame(dict((k, [v.get('property', ''), v['label']]) for k, v in smiles.items()), index=['mpo', 'label']).T
      else:
        unique_compounds = pd.DataFrame(dict((k, [v['label']]) for k, v in smiles.items()), index=['label']).T
        
      return render_template('index.html', title='Cheminformatic Analysis', 
        unique_compounds=unique_compounds.to_html(), num_compounds=len(unique_compounds), smiles=smiles, include_property=include_property)

    except Exception as e: 
      return render_template('index.html', title='Cheminformatic Analysis', errors=["Please input a valid file format"])
    
  return render_template('index.html', title='Cheminformatic Analysis')


@app.route("/cluster", methods=['GET'])
def upload():
  smiles = session['smiles']
  include_property = session['include_property']
        
  inputs = smiles.keys()

  #global all_smiles
  session['all_smiles'] = convert_to_smiles(smiles.copy())
  #global good_smiles
  good_smiles = convert_to_smiles(filter_smiles(pains.get_smiles(inputs), smiles))
  session['good_smiles'] = good_smiles
  #global bad_smiles
  bs = convert_to_smiles(pains.get_bad_smiles(inputs))
  bad_smiles = bs if isinstance(bs, dict) else {}
  session['bad_smiles'] = bad_smiles
  #global include_property
  session['include_property'] = include_property

  #global number of compounds
  session["num_remaining"] = len(good_smiles)
  session["num_removed"] = 0

  #global reasons_for_failure
  reasons_for_failure = dict.fromkeys(set(bad_smiles.values()), 0)
  for smile in bad_smiles.values():
    reasons_for_failure[smile] += 1

  session['reasons_for_failure'] = reasons_for_failure
  session.changed = True
    
  return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, num_remaining=session["num_remaining"], num_removed=session["num_removed"], reasons_for_failure=reasons_for_failure, include_property=session['include_property'])

@app.route('/verify_pains', methods=['GET', 'POST'])
def verify_pains():
  good_smiles = session.get('good_smiles')
  bad_smiles = session.get('bad_smiles')
  all_smiles = session.get('all_smiles')
  reasons_for_failure = session.get('reasons_for_failure')
  #TODO: check if reasons are still even the list? Does it matter?
  if request.method == 'POST':
    form = request.form
    if form['action'] == 'Remove Selected Compounds':
      for smile in form:
        if(smile != 'action'):
          del bad_smiles[smile]
          session['bad_smiles'] = bad_smiles
          session["num_removed"]+=1
          session.changed = True

    elif form['action'] == 'Keep Errors':
      for smile in form:
        if(smile != 'action'):
          try:
            del bad_smiles[smile]
            session['bad_smiles'] = bad_smiles
            good_smiles[smile] = all_smiles[smile]
            session['good_smiles'][smile] = all_smiles[smile]
            session['num_remaining']+= 1
            session.changed = True
          except Exception as e:
            print('Could not cluster {smile}')

  return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, num_remaining=session["num_remaining"], num_removed=session["num_removed"], reasons_for_failure=reasons_for_failure, include_property=session['include_property'])

@app.route('/verify_pains_by_error', methods=['GET', 'POST'])
def verify_pains_by_error():
  good_smiles = session.get('good_smiles')
  bad_smiles = session.get('bad_smiles')
  all_smiles = session.get('all_smiles')
  reasons_for_failure = session.get('reasons_for_failure')

  if request.method == 'POST':
    form = request.form
    if form['action'] == 'Remove Selected Errors':
      for reason in form:
        if(reason in reasons_for_failure.keys()):
          del reasons_for_failure[reason]
          session['reasons_for_failure'] = reasons_for_failure
          smiles_to_remove = []
          for (smile, smile_reason) in bad_smiles.items():
            if(reason == smile_reason):
              smiles_to_remove.append(smile)
          for smile in smiles_to_remove:
            del bad_smiles[smile]
            session["num_removed"]+=1
            session['bad_smiles'] = bad_smiles
            session.changed = True

    elif form['action'] == 'Keep Selected Errors':
      for reason in form:
        if(reason in reasons_for_failure.keys()):
          del reasons_for_failure[reason]
          session['reasons_for_failure'] = reasons_for_failure
          smiles_to_remove = []
          for (smile, smile_reason) in bad_smiles.items():
            if(reason == smile_reason):
              smiles_to_remove.append(smile)
          for smile in smiles_to_remove:
            try:
              del bad_smiles[smile]
              good_smiles[smile] = all_smiles[smile]
              session['good_smiles'][smile] = all_smiles[smile]
              session['num_remaining']+= 1
              session.changed = True
            except Exception as e:
              print('Could not cluster {smile}'.format(smile=smile))
            
  return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, num_remaining=session["num_remaining"], num_removed=session["num_removed"], reasons_for_failure=reasons_for_failure, include_property=session['include_property'])  

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
        return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure, errors=["Please input a valid tanimoto coefficient"], include_property=session['include_property'])
      elif ((reclusterCoefficient != '') and (float(tanimoto) < 0 or float(tanimoto) > 1)):
        return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure, errors=["Please input a valid recluster coefficient"], include_property=session['include_property'])
    except:
      return render_template('pains_verify_and_coefficient_use.html', title='Cheminformatic Analysis', bad_smiles=bad_smiles, reasons_for_failure=reasons_for_failure, errors=["Please input a valid tanimoto coefficient"], include_property=session['include_property'])

  color1 = request.form['mpo0Color']
  color1_array = color_hex_to_array(color1)

  fp_type = request.form['fp_radio']

  if(session['include_property']):
    color2 = request.form['mpo6Color'] 

    color2_array = color_hex_to_array(color2)
  else:
    color2 = color1
    color2_array = color1_array

  shouldRecluster = reclusterCoefficient != ''
  cluster = clustering.cluster(list(good_smiles.keys()), fp_type, 1 - float(tanimoto))
  if shouldRecluster :
    recluster_data = clustering.recluster_singletons(good_smiles, cluster, float(reclusterCoefficient), fp_type)
    recluster_smiles = recluster_data[0]
    recluster_clusters = recluster_data[1]
    tanimoto_smiles = clustering.get_tanimoto_coeffient_by_cluster(recluster_smiles, recluster_clusters, fp_type)
    result = get_smiles_json(tanimoto_smiles, float(tanimoto), recluster_clusters, session['include_property'], color1_array, color2_array, shouldRecluster)
  else :
    tanimoto_smiles = clustering.get_tanimoto_coeffient_by_cluster(good_smiles, cluster, fp_type)
    result = get_smiles_json(tanimoto_smiles, float(tanimoto), cluster, session['include_property'], color1_array, color2_array, shouldRecluster)

  include_property = session['include_property']
  session.clear()

  nodes_murcko = []
  nodes_type = []
  nodes_id = []
  nodes_label = []
  nodes_mpo = []
  nodes_reclustered = []

  for node in result['nodes']:
    nodes_murcko.append(node['murcko'])
    nodes_id.append(node['data']['id'])
    
    if node['data']['centroid']:
      nodes_type.append('centroid')
    else:
      nodes_type.append('non centroid')

    nodes_reclustered.append(node['data'].get('reclustered', ''))

    if node['data'].get('mpo'):
      nodes_mpo.append(node['data']['mpo'])
    
    nodes_label.append(node['data']['label'].split('\n')[0])

  nodes_cluster = [0] * len(nodes_id)
  
  for group, cluster in enumerate(result['clusterInfo']):
    for index, id in enumerate(nodes_id):
      if id in cluster:
        nodes_cluster[index] = group


  df = pd.DataFrame({'id': nodes_id, 'murcko': nodes_murcko, 'label': nodes_label, 'type': nodes_type, 'cluster': nodes_cluster, 'reclustered':nodes_reclustered})
  if nodes_mpo:
    df.insert(5, 'mpo', nodes_mpo)
  df.to_csv('export_data/output.csv')
  return render_template('cluster.html', title='Cheminformatic Analysis', color1=color1, color2=color2, include_property=include_property)

@app.route('/getPlotCSV')
def plot_csv():
    return send_file('data/export_data/output.csv',
                     mimetype='text/csv',
                     attachment_filename='output.csv',
                     as_attachment=True)

@app.after_request
def add_header(response):
  response.cache_control.public = True
  response.cache_control.max_age = 0
  return response

