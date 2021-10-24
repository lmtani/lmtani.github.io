---
title: Dependencies management with Conda
author: Lucas M. Taniguti
comments: true
date: 2021-10-15 11:33:00 +0800
categories: [Bioinformatics, Dependencies management]
tags: [bioinformatics, dependencies]
---

Compiling a bioinformatics program sometimes takes a lot of time and effort. When using conda, you will install a version that someone already has prepared the software, so it's easier to get it running on your machine.

You only need to download the version that satisfies your system specifications, such as Miniconda3 Linux 64-bit, if running on a standard Ubuntu desktop. Access the official Miniconda documentation to get links to other versions:

https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links

The installer will ask if you want conda as your default python environment. Usually, I accept it, so the conda version replaces the default python from my operational system.

> Note: conda does not require administrative permission, so the default python from your system will remain in the machine. Other users will not be affected.

After reloading your terminal, you'll see a new prefix in your prompt. Like this:

```shell
(base) ✔ ~ $
```
This 'base' informs you that you are using conda's base environment.

## Creating a new environment

Now imagine that you'll start working on one of your projects that have some particular dependencies. The first step is to create a new environment:

```shell
conda create -n name-for-your-new-env python=3.7

#Proceed ([y]/n)?  y
```

Here we are creating an environment called "name-for-your-new-env" that contains python 3.7. You need to activate it every time you want to work on it. Otherwise, you will stay in the base environment.

```sh
conda activate name-for-your-new-env
```

Now you can install all dependencies you need. There is a very useful channel to provide bioinformatics applications called "bioconda". Here is one example of how to search for a package in this channel:

```sh
(name-for-your-new-env) ✔ ~
14:37 $ conda search -c bioconda canu

Loading channels: done
# Name                       Version           Build  Channel
canu                             1.1               0  bioconda
canu                             1.3               0  bioconda
canu                             1.4               0  bioconda
canu                             1.4               1  bioconda
canu                             1.4      pl5.22.0_2  bioconda
canu                             1.5      pl5.22.0_0  bioconda
canu                             1.5      pl5.22.0_1  bioconda
canu                             1.6      pl5.22.0_1  bioconda
canu                             1.7      pl5.22.0_0  bioconda
canu                           1.7.1 pl526h470a237_0  bioconda
canu                             1.8      he1b5a44_1  bioconda
canu                             1.8      he1b5a44_2  bioconda
canu                             1.8 pl526h470a237_0  bioconda
canu                             1.9      he1b5a44_0  bioconda
canu                             1.9      he1b5a44_1  bioconda
canu                             2.0      he1b5a44_0  bioconda
canu                           2.1.1      h1b792b2_2  bioconda
canu                           2.1.1      he1b5a44_0  bioconda
canu                           2.1.1      he1b5a44_1  bioconda
canu                             2.2      ha47f30e_0  bioconda
```

To install a specific version of the [Canu Genome Assembler](https://github.com/marbl/canu) you can write:

```sh
conda install -c bioconda canu
```

And that is it. Now inside your conda environment called "name-for-your-new-env" you´ll have access to the package named canu.

To go back to the base environment, you can run `conda deactivate`.

If you want to remove the environment, you can run `conda env remove -n name-for-your-new-env`

> Note: see [mamba](https://github.com/mamba-org/mamba) if you think conda is too slow. After installing with `conda install mamba -n base -c conda-forge`, you'll be able to replace `conda` to `mamba` in the examples provided in this post.

## Using a provided environment

For the post [Playing with bioinformatics I]({% post_url 2021-10-13-playing-with-bioinformatics-I %}) I've prepared a conda environment and tests for my script. Now you can reproduce the whole process in your machine without too much effort in installing dependencies.

Follow these steps:

```sh
# 1. Clone my repository, so you'll have the script and the test files.
git clone https://github.com/lmtani/lmtani.github.io.git

# 2. Go to the example directory
cd lmtani.github.io/_code/playing-with-bioinformatics-files/

# 3. Install the conda environment. You could use mamba instead
conda env create -f environment.yml

# 4. Activate it
conda activate bioinformatics-files

# 5. Install pytest-workflow to run automated tests
pip install pytest-workflow

# 6. Run pytest
cd ../../
pytest --kwd .
```

With these steps, you've just run the example from the script [playing-with-bioinformatics-files.sh](https://github.com/lmtani/lmtani.github.io/blob/main/_code/playing-with-bioinformatics-files/environment.yml). All dependencies are described in the [environment.yml](https://github.com/lmtani/lmtani.github.io/blob/main/_code/playing-with-bioinformatics-files/environment.yml) file.


Hope you liked this material. In the next post I plan to show how to use the pytest-workflow package.
