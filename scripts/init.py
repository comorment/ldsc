"""init.py

After creating this project from the containers_template repository,
run this interactive script once to initialize your personal repository
"""

import os
from copy import copy

if __name__ == '__main__':
    name = input(
        'What is your GitHub organization or username hosting this project? ')
    print(f'your GitHub organization/username is: {name}')

    project = input('What is the name of this project? ')
    print(f'your GitHub organization/username is: {project}')
    print(
        f'The project should be hosted at https://github.com/{name}/{project}')

    # old variables which will be renamed
    OLDNAME = 'comorment'
    OLDPROJECT = 'ldsc'

    while True:
        response = input('Is this correct (yes/no)? ')
        if response in ['Y', 'y', 'yes', 'Yes', 'N', 'n', 'no', 'No']:
            break
        else:
            print(f'{response} is not a valid response. Try again.\n')

    if response in ['N', 'n', 'no', 'No']:
        print('Exiting. No file changes were applied.\n')
    else:
        print('Converting repository....\n')

        # walk files and replace occurrences of `OLDPROJECT` by `project`` ID,
        # and `OLDNAME` by `name` (github org/user) as supplied by the user.
        forbiddendirs = ['.git', '.pytest_cache']
        exclude = set(['__pycache__', '_build', '_static', '_template'])
        forbiddenfiles = []
        for root, dirs, files in os.walk('.', topdown=True):
            dirs[:] = [d for d in dirs if d not in exclude]
            # iterate over files
            for filename in files:
                if root != '.':
                    if root.split(os.path.sep)[1] in forbiddendirs:
                        continue

                # modify file contents and rewrite:
                with open(os.path.join(root, filename),
                          'r', encoding="utf8") as f:
                    filedata = f.read()

                newfiledata = copy(filedata)
                newfiledata = filedata.replace(OLDNAME, name)
                newfiledata = newfiledata.replace(OLDPROJECT, project)

                if newfiledata != filedata:
                    print(f'rewriting {os.path.join(root, filename)}')
                    with open(os.path.join(root, filename),
                              'w', encoding="utf8") as f:
                        f.write(newfiledata)

                # modify file names:
                if filename.rfind(OLDPROJECT) >= 0:
                    newfilename = filename.replace(OLDPROJECT, project)
                    newfilepath = os.path.join(root, newfilename)
                    oldfilepath = os.path.join(root, filename)
                    print(
                        f'renaming {oldfilepath} {newfilepath}')
                    os.rename(os.path.join(root, filename),
                              os.path.join(root, newfilename))

            # iterate over directories and rename
            if root != '.':
                for directory in dirs:
                    if root.split(os.path.sep)[1] in forbiddendirs:
                        print(f'skipping {os.path.join(root, directory)}')
                        continue

                    if directory.rfind(OLDPROJECT) >= 0:
                        newdir = directory.replace(OLDPROJECT, project)
                        fullnewdir = os.path.join(root, newdir)
                        fullolddir = os.path.join(root, directory)
                        print(
                            f'renaming {fullolddir} {fullnewdir}')
                        os.rename(fullolddir, fullnewdir)

        # copy ./scripts/PROJECT_README.md over ./README.md
        os.remove('README.md')
        os.rename(os.path.join('scripts', 'PROJECT_README.md'), 'README.md')

        print('The repository has been converted.\n',
              'Commit and push all changes to the remote by issuing \n',
              '``git commit -a -m "initial setup"; git push``')
