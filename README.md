# cheminformatics

-------------------

Initial Setup:

- clone this repo
- create a conda environment: `conda create -n chem`
- add rdkit `conda install -c conda-forge rdkit`
- install flask: `conda install -c anaconda flask`
- `export FLASK_APP=cheminformatics.py`

Add PAINs script
- `git clone https://github.com/IanAWatson/Lilly-Medchem-Rules`
- move the contents of the folder to LillyMedchemRules and delete empty folder 'Lilly-Medchem-Rules'
- make sure ruby is installed on your machine
- run `make` to build the project

Add cytoscape
- in the main project directory, mkdir `static` 
- `cd static`
- `git clone https://github.com/cytoscape/cytoscape.js`
- now you have cytoscape!

And to run the site locally:
- `flask run`