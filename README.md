# LS_RemoTeC
## Installation
This version of RemoTeC handles dependencies using conda or miniconda. Whenever you want to recompile RemoTeC, you need to enter the `remotec-env` conda environment. If conda is installed (`conda --version`), go directory to step 2. If you do this and it turns out the setup you are working on does not allow you to create a conda environment, go back to step 1 and repeat.

### Step 1: Install Miniconda
```
mkdir -p ~/.local
cd ~/.local
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
when prompted for installation location, enter:
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

### Step 3: Activate Conda Environment
```
conda activate remotec-env
```
You are now ready to work on and recompile RemoTeC. Do this using `make clean; make`
