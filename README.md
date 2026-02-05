# EE492 Senior Design II — 5G NR PDSCH + HARQ Optimization

MATLAB-based link-level simulation framework for modeling the 5G New Radio (NR) downlink physical layer with emphasis on the Physical Downlink Shared Channel (PDSCH), Hybrid Automatic Repeat reQuest (HARQ), and machine-learning–enhanced retransmission strategies.

This project evaluates whether data-driven HARQ decision logic reliability, throughput, and latency can be improved by applying machine learning algorithms. 

---

## Repository Structure

.
├─ scripts/        Entry-point simulations and experiments  
├─ functions/      Reusable MATLAB functions  
├─ tests/          Verification/unit tests  
├─ data/           Local datasets (ignored by git)  
├─ results/        Generated plots and outputs  
├─ docs/           Notes, reports, and figures  
├─ README.md  
└─ .gitignore  

---

## Requirements

- MATLAB R2023a or newer (recommended)
- 5G Toolbox
- Communications Toolbox
- Statistics and Machine Learning Toolbox
- Git

Optional:
- WSL/Linux for development workflow

---

## Setup Instructions

### 1. Clone the repository

SSH:
git clone git@github.com:stewrad/ee492_5gNR_sim.git

HTTPS:
git clone https://github.com/stewrad/ee492_5gNR_sim.git

---

### 2. Open in MATLAB

Open MATLAB and change directory:

cd path/to/ee492_5gNR_sim

---

### 3. Add project paths

Run once per session:

addpath(genpath(pwd));

(Optional permanent)
savepath

---

## Running Simulations

Baseline HARQ:
run("scripts/nr_dlsch_tx_rx.m")

<!-- ML-enhanced HARQ:
run("scripts/run_harq_ml.m") -->

<!-- Quick demo:
run("scripts/demo.m") -->

<!-- Outputs are saved in:
results/ -->

---

## Development Workflow

### Branching (required)

Do NOT commit directly to main.

Create a feature branch:

git checkout -b feature/your-feature-name  
git push -u origin feature/your-feature-name  

Open a Pull Request to merge.

---

## Coding Guidelines

- Place reusable logic inside /functions
- Keep scripts short and modular
- Document all functions
- Avoid hardcoded file paths
- Store large datasets outside the repository
- Keep commits small and descriptive

---

## Data Policy

Large files should NOT be committed to git.

Ignored examples:
- .mat
- .bin
- .h5
- raw IQ captures
- generated results

Use cloud storage or shared drives for large datasets.

---

## WSL Users (optional)

You may store the repository inside WSL:

\\wsl.localhost\Ubuntu\home\<user>\repos\

Then open normally in MATLAB.

Note: heavy simulations often run faster when the repo is located on a native Windows path (C:\).

---

## Metrics Evaluated

- Block Error Rate (BLER)
- Throughput
- HARQ retransmission count
- Latency
- Spectral efficiency
- ML decision accuracy

---

## Contributors

- Austin Seward
- Liz Ramirez
- Arnold Kalala 

---

## License

Educational and research use only.
