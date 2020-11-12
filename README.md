# cheminformatics

-------------------

Initial Setup:

- clone this repo
- create a conda environment: `conda create -c rdkit -n chem rdkit` (on Windows, this will all be in the Anaconda Prompt)
- activate your environment with `conda activate chem`
- install flask: `pip install flask`
- install flask-excel `pip install Flask_Excel`
- install flask-session `pip install Flask_Session`
- `export FLASK_APP=cheminformatics.py`


Add PAINs script
- `git clone https://github.com/IanAWatson/Lilly-Medchem-Rules`
- move the contents of the folder to LillyMedchemRules and delete empty folder 'Lilly-Medchem-Rules'
- install docker and ensure that it is running
- pull the pains image: `docker pull ianwatson/lilly_medchem_rules:v1.2`

Set up Export Data Folder
- add a folder called export_data under app/data directory. This is where exported data csv files will be stored to access for export

And to run the site locally:
- `flask run`

To run tests:
- install pytest
- with the conda environment activated, `python -m pytest`