# LS_RemoTeC

## Installation
This version of RemoTeC handles dependencies using conda or miniconda. Whenever you want to compile RemoTeC, you need to enter the `remotec-env` conda environment.

### Step 1: Install Miniconda
If conda is installed (`conda --version`) and allows you to create your own environment, skip ahead to step 2.
Else, install miniconda as follows:
```
mkdir -p ~/.local
cd ~/.local
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
When prompted for installation location, enter:
```
$HOME/.local/miniconda3
```
Configure your shell to use Miniconda
```
echo 'export PATH="$HOME/.local/miniconda3/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### Step 2: Create Conda Environment
```
conda env create -f remotec-env.yml
```

## Building
RemoTeC needs to be built from the `remotec-env` conda environment. Activate it using
```
conda activate remotec-env
```
You are now ready to work on RemoTeC. To compile it, run `make`.
