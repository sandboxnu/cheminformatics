# cheminformatics

-------------------

Initial Setup:

- clone this repo
- create a conda environment with rdkit: `conda create -c rdkit -n my-rdkit-env rdkit`
- install flask: `conda install flask`
- `export FLASK_APP=cheminformatcis.py`

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