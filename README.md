# cheminformatics

-------------------

Initial Setup:

- clone this repo
- create a conda environment with rdkit: `conda create -c rdkit -n my-rdkit-env rdkit`
- install flask: `conda install flask`

Add PAINs script
- first, clone this repo in this directory: https://github.com/IanAWatson/Lilly-Medchem-Rules/blob/master/README
- rename the folder to be LillyMedchemrules (to remove dashes)
- enter the directory and place `pains.py` in it
- make sure ruby is installed on your machine
- run `make` to build the project
- run `./pains.py`
