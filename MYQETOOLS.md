# myqetools API Guide

Comprehensive guide for the object-oriented helpers defined in `myqetools.py`. The toolkit streamlines post-processing of Quantum ESPRESSO (QE) calculations, covering band structures, effective masses, projected density of states (PDOS), charge densities, and energetic analyses. Use this document to understand each class individually and to compose them into larger workflows.

> **TL;DR workflow**  
> 1. Load raw QE outputs (band `.gnu`, PDOS `.dat`, charge-density `.pp`, SCF logs, etc.).  
> 2. Instantiate the relevant class (for example `BandStructure`, then `EffectiveMassAnalyzer`).  
> 3. Call plotting or analysis helpers to visualise or extract metrics.  
> 4. Combine results in notebooks, scripts, or GitHub automation.

---

## Prerequisites

- **Python**: 3.9 or newer recommended.
- **Dependencies**: `numpy`, `matplotlib`, `scipy`, `ase`.
  - Install via `pip install numpy matplotlib scipy ase`.
- **QE outputs**:
  - Band structure: two-column text (`k`, `energy`) from `bands.x` (`*.bands.gnu`).
  - PDOS: `.dat` files from `projwfc.x`.
  - Charge density: `.pp` data exported to plain text with a 50-line header.
  - SCF or relaxation logs from `pw.x`.
  - QE structure files in `espresso-in` / `espresso-out` format.

Ensure your environment supports Matplotlib rendering (Jupyter, VS Code, or a desktop Python session).

---

## Quick Start Example

```python
import numpy as np
from myqetools import BandStructure, EffectiveMassAnalyzer

# 1) Load band structure data exported by bands.x (.bands.gnu format)
fermi = 5.6925
data = np.loadtxt("CsGeCl.bands.gnu")
high_symmetry = [("Gamma", 0.0), ("M", 0.7071), ("R", 1.2071), ("X", 1.9142)]

# 2) Build the band structure object and draw the plot
bs = BandStructure()
bs.get_bandstructure(data, fermi, high_symmetry)
bs.plot_shifted_band_structure(start_band_number=23, shift_value=1.827)  # optional rigid shift
bs.save_figure("band_structure.png")

# 3) Analyse effective masses around the extrema
analyzer = EffectiveMassAnalyzer(bs)
me, mh, Eg = analyzer.analyze(
    VBM_main=22,
    CBM_main=23,
    VBM_degenerate=[],
    CBM_degenerate=[24, 25],
    range_e=(35, 45),
    range_h=(35, 45),
)
print(f"Electron m* = {me:.3f} m_e, Hole m* = {mh:.3f} m_e, Eg = {Eg:.3f} eV")
```

---

## Class Reference

Each subsection documents one class from `myqetools.py`, including purpose, required inputs, and key methods.

### BandStructure

**Purpose**: Load and manipulate band energies sampled along a k-path.

- `get_bandstructure(data, fermi_energy, kpoints)`: Accepts a two-column array (`k`, `energy`). Internally reshapes to `(n_bands, n_kpoints)` and stores the Fermi level. `kpoints` is a list of `(label, coordinate)` pairs marking high-symmetry points.
- `plot_band_structure()`: Plots all bands shifted by the Fermi energy. X-axis ticks are labelled with high-symmetry points if supplied.
- `plot_band_number(band_number, overlay=True)`: Highlights an individual band. With `overlay=True`, draws on the existing figure; otherwise spawns a standalone plot.
- `shift_band(start_band_number, shift_value)`: Applies a rigid energy shift to bands at or above `start_band_number`. Use before plotting to mimic scissor corrections.
- `plot_shifted_band_structure(start_band_number, shift_value)`: Convenience wrapper that calls `shift_band` then `plot_band_structure`.
- `get_band_data_point(band_number)`: Returns `(k_values, energies_minus_fermi)` for further analysis.
- `save_figure(filename="band_structure.png", dpi=300, transparent=True)`: Writes the most recent figure to disk.

### EffectiveMassAnalyzer

**Purpose**: Compute electron and hole effective masses using curvature fits near the band extrema.

- Instantiate with an existing `BandStructure`: `analyzer = EffectiveMassAnalyzer(bs)`.
- `analyze(VBM_main, CBM_main, VBM_degenerate=None, CBM_degenerate=None, range_e=(35, 45), range_h=(35, 45))`:
  - `VBM_main`, `CBM_main`: 1-indexed band numbers for valence-band maximum and conduction-band minimum.
  - `VBM_degenerate`, `CBM_degenerate`: optional lists of additional band indices that share the extremum.
  - `range_e`, `range_h`: tuples of `(start_idx, end_idx)` in k-space used for curvature averaging.
  - Returns `(electron_mass, hole_mass, band_gap)` and prints human-readable output.
- `get_band_gap(valence_band, conduction_band)`: Static helper that computes the gap between two arrays.
- Internally `_plot_effective_mass(...)` produces scatter/fit plots to validate the curvature region.

### EffectiveMass

**Purpose**: Lightweight average effective mass estimator for selected bands.

- Instantiate with `EffectiveMass(bandstructure)`.
- `calculate_average_mass(band_indices, k_range=None, plot=False)`:
  - Iterates over the provided bands, optionally slices the k-range, and computes second derivatives via `numpy.gradient`.
  - Setting `plot=True` visualises the selected band segments.
- `effmass(E, k)`: Core routine returning the effective mass (in electron masses) for a single `(E, k)` pair of arrays.

### PDOSPlotter

**Purpose**: Plot projected density of states (PDOS) data from `projwfc.x`.

- Instantiate with `PDOSPlotter(fermi_energy, orientation="vertical")`.
- `read_pdos(filepath)`: Returns `(energy_minus_fermi, pdos_array)` for that file.
- `infer_label(filepath)`: Generates a legend label such as `Pb(p)` based on the filename pattern.
- `plot(files, xlim=None, ylim=None, orientation=None, figsize=(8, 6))`:
  - `files` is a dictionary `{key: (color, filepath)}`.
  - Supports spin-collinear PDOS (two columns) by mirroring the spin-down channel when plotting vertically.
  - `orientation="horizontal"` swaps axes, useful when stacking beside band plots.

### OrbitalCenter

**Purpose**: Identify energy positions where a PDOS curve reaches minimal or maximal intensity.

- Requires a global variable `fermi_energy` defined before instantiation (`OrbitalCenter` reads it at module scope).
- `_load_data()`: Internal loader called during `__init__` to fill `energy` (shifted by Fermi level) and `pdos`.
- `get_extreme_positions(x_min=None, x_max=None)`:
  - Optionally restricts the PDOS window.
  - Returns `(energy_at_minimum, energy_at_maximum)` of the total PDOS.

### ChargeDensitySlice

**Purpose**: Visualise and compare three-dimensional charge density grids from QE.

- Instantiate with the FFT grid sizes and lattice constants:

  ```python
  cd = ChargeDensitySlice("charge.pp.txt", x1=120, x2=120, x3=120, a=11.95, b=11.95, c=11.95)
  ```

- The constructor converts charge density units from bohr^-3 to angstrom^-3 and reshapes the array.
- `plot_slice(hkl="001", position=None, interpolation="gaussian")`: Displays a 2D slice perpendicular to the chosen Miller plane. `position` defines the real-space offset along the corresponding lattice vector; defaults to the central plane.
- `plot_line_profile(axis="z", ix=None, iy=None, iz=None)`: Extracts and plots a 1D cut along the selected axis; other indices default to midpoints.
- `get_line_profile(axis="z", index1=None, index2=None)`: Returns `(coordinates, profile)` without plotting.
- `plot_charge_density_difference(file1, file2, x1, x2, x3, a, b, c, axis="z", index1=None, index2=None, label1="file1", label2="file2", diff_label="Difference")`: Class-level helper that compares two charge-density files along the same line cut.

### QeInputFile

**Purpose**: Compute interatomic distances from QE input files using ASE.

- Accepts a single path or list of paths: `qe = QeInputFile(["scf.in", "relax.in"])`.
- `get_atom_indices(atoms, symbol)`: Return indices of atoms matching the element symbol.
- `get_all_distances(atoms, symbol1, symbol2)`: Compute every pairwise distance between the specified species (excluding self-pairs).
- `print_average_distances(symbol1, symbol2)`: For each file, prints the average of the distances computed above.

### Minimaout

**Purpose**: Inspect relaxed geometries from QE output (`espresso-out`) files.

- Instantiate with `Minimaout("scf.out")` or any geometry-containing log.
- `get_atomic_positions(exclude_element="O")`: Returns Cartesian and fractional coordinates for all atoms except the excluded species.
- `find_highest_z_coordinates(exclude_element="O")`: Identifies the top two atoms along the z-direction and reports their coordinates together with the z-separation.

### Minimain

**Purpose**: Analyse initial structures from QE input (`espresso-in`) files.

- Methods mirror `Minimaout`:
  - `get_atomic_positions(exclude_element="O")`.
  - `find_highest_z_coordinates(exclude_element="O")`.
- Use it to compare starting and final geometries within the same workflow.

### TotalEnergy

**Purpose**: Extract the final total energy from QE SCF logs.

- Instantiate with `TotalEnergy("scf.out")`.
- `extract()`: Stream through the file and capture the last line containing `"!    total energy"`, storing the value in Ry.
- `last_TotalEnergy()`: Print the stored energy in both Ry and eV. Requires a prior call to `extract()`.

### SurfaceEnergyCalculator

**Purpose**: Convert slab energies into surface energies.

- Instantiate optionally overriding the Ry-to-J/m² factor: `SurfaceEnergyCalculator(ry_to_jm2=217.9863)`.
- `surf(data, bulk, N, area)`:
  - `data`: iterable of slab energies (Ry).
  - `bulk`: bulk energy per formula unit (Ry).
  - `N`: number of bulk units in the slab.
  - `area`: iterable of surface areas (square angstrom).
  - Returns surface energies in J/m².

### AdsorptionCalculator

**Purpose**: Evaluate adsorption energies and Arrhenius-type desorption times.

- Instantiate with `AdsorptionCalculator(Eadsorbate, ry_to_eV=13.605698066, kBT=0.0258519)`.
- `Eadsorption(Esys, Eadsorbant)`: Returns `(Esys - (Eadsorbant + Eadsorbate)) * ry_to_eV` in eV.
- `desorption_time(t, Esys, Eadsorbant)`: Multiplies the attempt time `t` by `exp(-E_ads / kBT)` using the adsorption energy from `Eadsorption`.

---

## Workflow Examples

1. **Electronic structure**: Use `BandStructure` to load `.bands.gnu`, highlight specific bands with `plot_band_number`, and feed the instance into `EffectiveMassAnalyzer` or `EffectiveMass` for transport metrics.
2. **Density of states**: Instantiate `PDOSPlotter` to overlay orbital contributions, while `OrbitalCenter` pinpoints the energy centroids of selected orbitals.
3. **Real-space analysis**: Load charge-density grids with `ChargeDensitySlice`, examine planar slices, and compare pristine versus doped structures using `plot_charge_density_difference`.
4. **Geometry review**: Combine `QeInputFile`, `Minimaout`, and `Minimain` to inspect interatomic distances and surface atom positions across input and output geometries.
5. **Energetics**: Summarise SCF runs via `TotalEnergy`, estimate surface energies with `SurfaceEnergyCalculator`, and study adsorption thermodynamics with `AdsorptionCalculator`.

---

## Tips and Extensibility

- **Unit discipline**: QE outputs mix Ry and eV. Confirm each method's expected units before combining results.
- **Matplotlib styling**: The module sets `Times New Roman` globally. Override `plt.rcParams` before importing `myqetools` if you prefer a different style.
- **Interpreter availability**: If your system lacks a `python` executable, create a virtual environment (for example `conda create -n qe python=3.10`) and install the required packages there.
- **Custom extensions**: All classes expose NumPy arrays, making it easy to integrate additional fitting, machine learning, or plotting routines on top.

For contributions or questions, open an issue or pull request in the repository and reference this guide.
