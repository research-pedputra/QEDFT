# Quantum Confinement Effect: Band gap energy fitting using Brus Equation

This repository contains a Python script that fits the quantum confinement effect on the band gap of nanostructured SnO$_2$. The script uses a custom function to model the relationship between the thickness of the nanostructure and its direct allowed transition energy gap.

## Getting Started

### Prerequisites

Before running the script, ensure you have the following Python packages installed:

- `numpy`
- `matplotlib`
- `scipy`

You can install these packages using pip:

```bash
pip install numpy matplotlib scipy
```
### Running the Script

1. **Prepare your data**: The script uses predefined data for the thickness of nanostructures and their corresponding energy gaps. You can modify the `R` and `y_data` variables in the script to use your own data.

2. **Run the script**: The script performs curve fitting to the quantum confinement model, plots the original data along with the fitted curve, and prints the fitting parameters, R-squared, and chi-square values.

```bash
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams['figure.figsize'] = [8, 6]

def gap_qc(R, gap_qc, C1, C2):
    energy_gap = gap_qc + C1 / (R**2) - C2 / R
    return energy_gap

R = np.array([199, 250, 300, 400, 500])  
y_data = np.array([3.82, 3.74, 3.68, 3.64, 3.62])  # Energy gap data in eV

# Initial guess for the fitting parameters
initial_guess = [1.6, 1e-3, 1e-9]

# Perform curve fitting
fit_params, _ = curve_fit(gap_qc, R, y_data, p0=initial_guess)

# Extract the fitting parameters
gap_qc_fit, C1_fit, C2_fit = fit_params

# Generate the fitted curve
R_fit = np.linspace(min(R), max(R), 100)
y_fit = gap_qc(R_fit, gap_qc_fit, C1_fit, C2_fit)

# Plot the original data and the fitted curve
plt.scatter(R, y_data, label='Nanostructured SnO$_2$')
plt.plot(R_fit, y_fit, label='Fitted', color='r')
plt.xlabel('Thickness (nm)', fontsize=14, fontweight='bold')
plt.ylabel('Direct Allowed Transition Energy gap (eV)', fontsize=14, fontweight='bold')
plt.tick_params(axis='y', which='major', labelsize=14)
plt.tick_params(axis='x', which='major', labelsize=14)
plt.legend()
plt.show()

# Print the fitting parameters
print('Fitted parameters:')
print('gap_qc =', gap_qc_fit)
print('C1 =', C1_fit)
print('C2 =', C2_fit)

# Calculate residuals and goodness-of-fit
residuals = y_data - gap_qc(R, gap_qc_fit, C1_fit, C2_fit)
residual_sum_squares = np.sum(residuals ** 2)
total_sum_squares = np.sum((y_data - np.mean(y_data)) ** 2)
R2 = 1 - (residual_sum_squares / total_sum_squares)
chi2 = np.sum(residuals ** 2) / (y_data.size - len(fit_params))

# Print R-squared and chi-square
print('R-squared:', R2)
print('chi-square:', chi2)
```

### Explanation of Results

1. **Fitted Parameters**:
   - `gap_qc`: The quantum confinement gap.
   - `C1` and `C2`: Coefficients related to the thickness dependence.

2. **R-squared**: Indicates the goodness of fit, with values closer to 1 indicating a better fit.

3. **Chi-square**: Provides an additional measure of fit quality.

