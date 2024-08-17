# Using ASE-Python: Equation of State

This repository contains a Python script that performs an Equation of State (EOS) fitting on energy-volume data for Nickel Oxide (NiO). The fitting process helps in deriving important mechanical properties, such as bulk modulus, from the energy-volume relationship. The script uses the `ASE` (Atomic Simulation Environment) library's `EquationOfState` class to fit the data and visualize the results.

## Getting Started

### Prerequisites

Before running the script, ensure you have the following Python packages installed:

- `ase`
- `matplotlib`
- `pandas`
- `numpy`

You can install these packages using pip:

```bash
pip install ase matplotlib pandas numpy
```

### Running the Script
1. Load your data
- The data must be ib two-column style, contains the volume and energy. 
- Use any python llibrary to load, numpy or even pandas should be fine. 
2. Lets assume you are working with jupyter notebook, then just simply follow this code

```bash
from ase.eos import EquationOfState
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
plt.rcParams.update({'font.family':'Times New Roman'})

# Load energy-volume data from Excel
dataeos = pd.read_excel('NiO_EOS_New.xlsx')
volumes = list(dataeos['Volume (a.u.)'])
energies = list(dataeos['Calculated'] * 13.6056980659)

# Shift the energy data so that the minimum energy corresponds to zero
emin = np.min(energies)
energies -= emin

# Fit the equation of state
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()

# Plot the fit with the shifted energy data
plt.plot(volumes, energies, 'o', markersize=10)
eos.plot()

# Add y-axis label with shifted energy units
plt.xlabel('Volume (a.u.)$^3$', fontsize=14, fontweight='bold')
plt.ylabel('E - E$_{min}$ (eV)', fontsize=14, fontweight='bold')
plt.tick_params(axis='x', which='major', labelsize=14)
plt.tick_params(axis='y', which='major', labelsize=14)
plt.show()
```
