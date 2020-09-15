# cheminformatics

-------------------

Initial Setup:

- clone this repo
- create a conda environment: `conda create -c rdkit -n chem rdkit`
- activate your environment with `conda activate chem`
- install flask: `pip install flask`
- install flask-excel `pip install Flask_Excel`
- install flask-session `pip install Flask_Session`
- `export FLASK_APP=cheminformatics.py`


Add PAINs script
- `git clone https://github.com/IanAWatson/Lilly-Medchem-Rules`
- move the contents of the folder to LillyMedchemRules and delete empty folder 'Lilly-Medchem-Rules'
- install docker 
- pull the image: `docker pull ianwatson/lilly-medchem-rules:v1.2`

And to run the site locally:
- `flask run`