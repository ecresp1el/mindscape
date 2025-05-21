##  Environment Setup

This project uses a Conda-based environment to ensure reproducibility. Setup instructions vary slightly depending on whether you're working on a local machine or the University of Michigan Great Lakes HPC cluster.

### Step 1: Clone the repository

Run this from your terminal:
git clone https://github.com/ecresp1el/MindScape.git cd MindScape


---

### Step 2: Choose your setup environment

#### Option A: Personal machine (e.g., laptop or workstation)

Use this if you have **Anaconda or Miniconda installed** on your system.

Run:

bash scripts/setup_env_local.sh

This script will:
- Check that `conda` is available
- Create the `mindscape-env` environment from `conda-environments/MINDSCAPE.yml`
- Activate the environment
- Verify the installation of required Python packages (e.g., `scanpy`)

Once setup completes, activate the environment manually if needed:
conda activate mindscape-env

---

#### Option B: UMich Great Lakes HPC cluster

Use this if you're logged into the Great Lakes system.

Run:
bash scripts/setup_env_greatlakes.sh

This script will:
- Load required UMich modules: `python3.11-anaconda`, `Bioinformatics`, and `cellranger`
- Create the same `mindscape-env` environment from `conda-environments/MINDSCAPE.yml`
- Activate the environment
- Verify Python libraries and Cell Ranger availability

> Note: If you're using Great Lakes but forget to load the Anaconda module, the local script will attempt to load it for you automatically.

---

### What these scripts do

- Remove any existing `mindscape-env` environment to ensure a clean setup
- Recreate it from the exact specifications defined in the YAML file
- Perform basic tests to confirm successful setup
- They do **not** perform any data analysis or modify project files

---
## Re-entering the Environment (e.g. in VS Code Terminal)

If you're using VS Code or any new terminal session on the UMich Great Lakes cluster:

Run the following script each time you open a terminal:

```bash
bash scripts/load_env_interactive.sh

This will:

    Load the correct Anaconda module

    Activate the mindscape-env Conda environment


---

##  Daily Use: Re-Entering Your Environment on Great Lakes

Each time you open a new terminal (e.g. in VS Code), you must **re-activate** the environment.

To do that, run:

```bash

source scripts/load_env_interactive.sh




