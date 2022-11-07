# Container `ldsc.sif`

Specific information about the `ldsc.sif` file goes here.
Modify as needed.

You may use the ``ldsc.sif`` container file to familirize yourself with [Singularity](<https://sylabs.io/docs/>),
and the way it works on your secure HPC environment (TSD, Bianca, Computerome, or similar).
This singularity container is indented as a demo.
It only contains Python 3.10.

## Getting Started

* Download ``ldsc.sif`` from [here](https://github.com/comorment/ldsc/tree/main/containers)
* Import these files to your secure HPC environment
* Run ``singularity exec --no-home ldsc.sif python3 --help``, to validate that you can run singularity. This command is expected to produce the standard python help message, starting like this:

```
usage: /usr/local/bin/python [option] ... [-c cmd | -m mod | file | -] [arg] ...
Options and arguments (and corresponding environment variables):
-b     : issue warnings about str(bytes_instance), str(bytearray_instance)
         and comparing bytes/bytearray with str. (-bb: issue errors)
...
```

## Helpful links to singularity documentation

It's good idea to familirize with basics of Singularity, such as these:

* ["singularity shell" options](https://sylabs.io/guides/3.2/user-guide/cli/singularity_shell.html#options)
* [Bind paths and mounts](https://sylabs.io/guides/3.2/user-guide/bind_paths_and_mounts.html).

## Installing Docker and Singularity on your local machine

While you're getting up to speed with singularity, it might be reasonable to have it install on your local machine (laptop or desktop),
and try out containers locally before importing them to your HPC environment.

To install singularity on Ubuntu follow steps described here: <https://sylabs.io/guides/3.7/user-guide/quick_start.html>
Note that ``sudo apt-get`` can give a very old version of singularity, which isn't sufficient.
Therefore it's best to build singularity locally.  Note that singularity depends on GO, so it must be installed first.
If you discovered more speciifc instructions, please submit an issue or pull request to update this documentation.

## Mapping your data to singularity containers

There are several ways to give singularity container access to your data. Here are few examples (extend as required by the project):

1. The command
   
   ```
   singularity exec --home $PWD:/home ldsc.sif python -c "import os; print(f'Hello World from {os.getcwd()}')"
   ```
   
   will map your present workding directory (`$PWD`) into ``/home`` directory within the container, and set it as active working directory.
   In this way in your python command you can refer to the files as if they are in your local folder without specifying the path.

2. Now, let's assume that instead of downloading the ``ldsc.sif`` container file you've cloned the entire github repo
   (``git clone git@github.com:comorment/ldsc.git``), have built and pushed the container, and have transfered it to your HPC environment.
   Then change your folder to the root of the ``ldsc`` repository, and run these commands:

   ```
   singularity exec --bind reference/:/ref:ro,out_dir:/out:rw singularity/ldsc.sif python ...
   ```

   Note that input paths are relative to the current folder. Also, we specified ``ro`` and ``rw`` access, to have reference data as read-only,
   but explicitly allow the container to write into ``/out`` folder (mapped to ``out_dir`` on the host).

3. Run 

   ```
   singularity shell --home $PWD:/home -B $(pwd)/data:/data ldsc.sif
   ``` 
   
   to use singularity in interactive mode.
   In this mode you can interactively run python commands.
   Note that it will consume resources of the machine where you currently run the singularity  comand
   (i.e., most likely, the login node of your HPC cluster). 
   It is usually possible to use the compute nodes in interactive nodes (cf. your system documentation).

## Running as SLURM job

* Run singularity container within SLURM job scheduler, by creating a ``ldsc_job.sh`` file (by adjusting the example below), and running ``sbatch ldsc_job.sh``:

  ```
  #!/bin/bash
  #SBATCH --job-name=ldsc
  #SBATCH --account=$PROEJCT
  #SBATCH --time=00:10:00
  #SBATCH --cpus-per-task=1
  #SBATCH --mem-per-cpu=8000M
  module load singularity/3.7.1
  singularity exec --no-home ldsc.sif python -c "print('Hello World')"
  singularity exec --home $PWD:/home ldsc.sif python -c "print('Hello Moon')"
  ```

Please [let us know](https://github.com/comorment/ldsc/issues/new) if you face any problems.
