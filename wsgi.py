from app import app
import os



if __name__ == '__main__':
    extra_files = ['./app/static/js/data.json',]
    json_file = 'app/static/js/data.json'
    pains_files = ['bad0.smi', 'bad1.smi', 'bad2.smi', 'bad3.smi', "input.smi", "okmedchem.smi"]
    
    for pains in pains_files:
        dir_name = './scripts/pains/LillyMedchemRules/'
        file_name = os.path.join(dir_name, pains)
        print(file_name)
        if os.path.isfile(file_name):
            extra_files.append(file_name)
    print(os.getcwd())


    app.run(extra_files=extra_files)

