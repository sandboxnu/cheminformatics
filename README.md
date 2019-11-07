# cheminformatics

-------------------

Initial Setup:

- clone this repo
- create a conda environment: `conda create -n chem`
- add rdkit `conda install -c conda-forge rdkit`
- install flask: `conda install -c conda-forge flask`
- install flask-excel `conda install -c conda-forge flask-excel`
- `export FLASK_APP=cheminformatics.py`

Add PAINs script
- `git clone https://github.com/IanAWatson/Lilly-Medchem-Rules`
- move the contents of the folder to LillyMedchemRules and delete empty folder 'Lilly-Medchem-Rules'
- make sure ruby is installed on your machine
- run `make` to build the project

And to run the site locally:
- `flask run`