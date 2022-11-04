# Project ``container_template``

This is a template repository for singularity/dockerfile containers and build scripts

## Build status

[![License](http://img.shields.io/:license-GPLv3+-green.svg)](http://www.gnu.org/licenses/gpl-3.0.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7244728.svg)](https://doi.org/10.5281/zenodo.7244728)
[![Documentation Status](https://readthedocs.org/projects/container-template/badge/?version=latest)](https://container-template.readthedocs.io/en/latest/?badge=latest)
[![Flake8 lint](https://github.com/precimed/container_template/actions/workflows/python.yml/badge.svg)](https://github.com/precimed/container_template/actions/workflows/python.yml)
[![Dockerfile lint](https://github.com/precimed/container_template/actions/workflows/docker.yml/badge.svg)](https://github.com/precimed/container_template/actions/workflows/docker.yml)

## Citation

If you use this material for scientific publications, please cite it as follows:

    Espen Hagen, Bayram Cevdet Akdeniz, Alexandr Frei. (2022). precimed/container_template: container_template-0.1.0rc0 (v0.1.0rc0). Zenodo. https://doi.org/10.5281/zenodo.7244728

## Getting started

To use these codes, **do not** clone or fork this repository, but click the green
[Use this template](https://github.com/precimed/container_template/generate)
button above.

You may then create a new repository on a GitHub organization or user account (`<user>`) under a new project name (`<project>`),
that can be cloned by issuing in a terminal application

```
git clone git@github.com:<user>/<project>.git
cd <project>
```

The codes may then be adapted further to suit the requirements of the project, as described next

## Setup-requirements

After cloning the template repository as described above a few requirements must be in place.
Some scripts rely on [Python](https://www.python.org) with some packages available in the Python environment.

One convenient way to set up a Python environment is [Conda](https://docs.conda.io/en/latest/), available via [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Miniforge](https://github.com/conda-forge/miniforge).
Assuming one of these Python distributions is installed, issue the following to create a new Python environment using `conda`:

```
conda create -n <container_template> python=3 pip  # create environment named <container_template>
conda activate <container_template>  # activates enviromnent <container_template>
pip install -r setup-requirements.txt
```

## Initial setup

To set up your own project from this template, issue in the terminal:

```
python scripts/init.py
```

The script is interactive, and will prompt the user for some info.
It will modify certain files in the `<project>` directory and also replace this `README` file.
Some files may also be renamed.
Remember to commit and push the changes after running the setup script (`scripts/init.py`):

```
git add <file1> <file2> ...
git commit -m "initial setup"
git push
```

The remaining codes may then be added to and be modified further to suit the requirements of the `<project>`.
See the updated project README.md file for further instructions.
