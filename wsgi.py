from app import app

from os import path, walk




if __name__ == '__main__':
    extra_dirs = ['directory/to/watch',]
    extra_files = extra_dirs[:]
    for extra_dir in extra_dirs:
        for dirname, dirs, files in walk(extra_dir):
            for filename in files:
                filename = path.join(dirname, filename)
                if path.isfile(filename):
                    extra_files.append(filename)

    app.run(extra_files=extra_files)
    app.run()

