# Vehicle Routing Problem Solver

This Python code solves the Vehicle Routing Problem (VRP) using Tabu Search, aiming to find an optimal solution for routing a fleet of vehicles to service customers within certain constraints.

## Overview

The project consists of several classes and methods:

- **`Customer`**: Represents a customer with specific attributes such as ID, coordinates, demand, ready time, due time, service time, and distance to the depot.
- **`Vehicle`**: Represents a vehicle with attributes like ID, a list of customers it serves, capacity, time, distance, and time over.
- **`Solution`**: Manages the solution composed of vehicles and customers, handling operations like adding/removing customers from vehicles, swapping customers, and generating neighbors.
- **`Problem`**: Handles the VRP problem instance, calculates fitness, checks feasibility, and manages parameters for the tabu search.

## Usage

### Installation

1. Clone the repository or download the source code.

2. Ensure you have Python installed.

3. Open the terminal or command prompt and navigate to the directory containing the code.

### Running the Code

Run the script providing the input data file:

```bash
python main.py data/i5.txt
Replace data/i5.txt with the path to your input data file.

Results
The output will be written to a file named according to the runtime in the results/ directory.
```

# HMO Validator

* install and use Python >= 3.7

## INSTALL
```bash
pipenv install
source <virtualenv_activate_path> # (reported by virtualenv)
```
or
```bash
# create virtualenv manually
pip install
source <virtualenv_activate_path> # (created manually)
```
or
```bash
# globally
pip3 install click==7.0
pip3 install numpy==1.18.0 # most other versions should be fine
# and now run with global python3.7
```

## USAGE
```bash
python3.7 validator.py -i ../instances/i1.txt -o ../out/res-1m-i1.txt
```
or see `help`
```bash
python3.7 validator.py --help
```
