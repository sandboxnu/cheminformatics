from app import app

from os import path, walk




if __name__ == '__main__':
    extra_files = []
    json_file = 'app/static/js/data.json'
    extra_files.append(json_file)
    pains_files = ['bad0.smi', 'bad1.smi', 'bad2.smi', 'bad3.smi', "input.smi", "okmedchem.smi"]
    
    for pains in pains_files:
        dir_name = '/scripts/pains/LillyMedchemRules/'
        extra_files.append(dir_name + pains_files)

    """
    for extra_dir in extra_dirs:
        for dirname, dirs, files in walk(extra_dir):
            for filename in files:
                filename = path.join(dirname, filename)
                if path.isfile(filename):
                    extra_files.append(filename)
    """
    app.run(extra_files=extra_files)
    app.run()

