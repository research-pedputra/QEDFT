# myqetools API Guide

Comprehensive guide for the object-oriented helpers defined in `myqetools.py`. The toolkit streamlines post-processing of Quantum ESPRESSO (QE) calculations, covering band structures, effective masses, projected density of states (PDOS), charge densities, and energetic analyses. This document explains every public class, the data each method expects, and practical usage patterns so you can integrate the utilities into notebooks, scripts, or reproducible GitHub workflows.

> **TL;DR workflow**  
> 1. Load raw QE outputs (band `.gnu`, PDOS `.dat`, charge-density `.pp`, SCF text, etc.).  
> 2. Instantiate the relevant class (for example `BandStructure` followed by `EffectiveMassAnalyzer`).  
> 3. Call plotting or analysis helpers to visualise or extract metrics.  
> 4. Combine results in your reports or automation scripts.

---

## 1. Prerequisites

- **Python**: 3.9+ recommended.
- **Dependencies**: `numpy`, `matplotlib`, `scipy`, and `ase` (`pip install numpy matplotlib scipy ase`).
- **QE outputs**:
  - Band structure: plain text with two columns (`k`, `energy`) as produced by `bands.x`.
  - PDOS: `.dat` files created by `projwfc.x`.
  - Charge density: `.pp` files (converted to plain text with a 50-line header).
  - SCF/relax output: `pw.x` stdout log.
  - QE input (`*.in`) or output (`*.out`) structures for geometry utilities.

All plotting functions rely on Matplotlib’s current backend; run in an environment that supports figure rendering (Jupyter, VS Code, or a desktop Python session).

---

## 2. Quick Start Example

```python
import numpy as np
from myqetools import BandStructure, EffectiveMassAnalyzer

# 1) Load band structure data exported from bands.x (.bands.gnu format)
fermi = 5.6925
data = np.loadtxt("CsGeCl.bands.gnu")
high_symmetry = [("Γ", 0.0), ("M", 0.7071), ("R", 1.2071), ("X", 1.9142)]

# 2) Prepare the band structure object and draw the plot
bs = BandStructure()
bs.get_bandstructure(data, fermi, high_symmetry)
bs.plot_shifted_band_structure(start_band_number=23, shift_value=1.827)  # optional rigid shift
bs.save_figure("band_structure.png")

# 3) Analyse sub-gap effective masses around the extrema
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

The remaining sections describe every class in detail so you can chain analyses beyond the example above.

---

## 3. Band-Resolved Analysis

### 3.1 `BandStructure`

Purpose: Load band energies sampled along a k-path and provide convenient plotting and slicing utilities.

| Method | Description |
| --- | --- |
| `get_bandstructure(data, fermi_energy, kpoints)` | Accepts a 2-column array (`k`, `energy`) and reshapes it into `(nbands, nk)` for internal use. `kpoints` is a list of `(label, path_coordinate)` pairs marking high-symmetry points. |
| `plot_band_structure()` | Plots every band shifted by the stored Fermi level. Axis limits default to ±5 eV around `E_f`. |
| `plot_band_number(band_number, overlay=True)` | Highlights an individual band, either overlayed on the main figure or in a separate plot. |
| `shift_band(start_band_number, shift_value)` / `plot_shifted_band_structure(...)` | Apply a rigid energy shift to bands at or above `start_band_number` before plotting (useful for scissor corrections). |
| `get_band_data_point(band_number)` | Returns `(k_values, energies)` for the requested band with the Fermi level removed. |
| `save_figure(filename, dpi=300, transparent=True)` | Persist the latest band structure figure. |

> **Tip**: Load GMO/bands data via `numpy.loadtxt`. Ensure the rows are ordered by band index and k-point, exactly as exported by QE (`bands.x`).

### 3.2 `EffectiveMassAnalyzer`

Purpose: Compute electron and hole effective masses by averaging curvature of the selected conduction and valence bands.

- Initialise with an existing `BandStructure` instance: `analyzer = EffectiveMassAnalyzer(bs)`.
- Call `analyze(...)` with the band indices and the k-index ranges that bracket the parabolic region near the extremum. Parameters:
  - `VBM_main`, `CBM_main`: primary band numbers for valence and conduction edges (1-indexed, matching the `.bands.gnu` ordering).
  - `VBM_degenerate`, `CBM_degenerate`: optional lists of additional degenerate bands to average.
  - `range_e`, `range_h`: tuples `(start_index, end_index)` specifying the slice in k-space used for the curvature fit.
- Returns `(m_e*, m_h*, Eg)` and prints a report. Internal helper `_plot_effective_mass` visualises the scatter and spline fit for each side.
- Use `get_band_gap(valence_band, conduction_band)` if you only need the gap from arrays.

### 3.3 `EffectiveMass`

Purpose: Lightweight alternative that estimates an average effective mass across several bands.

- Instantiate with the same `BandStructure` object.
- `calculate_average_mass(band_indices, k_range=None, plot=False)` loops over the selected bands, optionally clips the k-range (`(start, end)` indices), and returns the mean effective mass computed by finite differences (`np.gradient`). Set `plot=True` to inspect the energy slices used for the fit.

---

## 4. Density of States & Orbital Analysis

### 4.1 `PDOSPlotter`

Plot total or spin-resolved projected DOS.

- Construct with a Fermi level in eV: `plotter = PDOSPlotter(fermi_energy, orientation="vertical")`.
- `plot(files, xlim=None, ylim=None, orientation=None, figsize=(8, 6))` expects a dictionary of the form:

  ```python
  files = {
      "Pb_p": ("red", "Pb_p.dat"),
      "Br_p": ("green", "(Br)_(p).dat"),
  }
  plotter.plot(files, xlim=(-5, 5))
  ```

- The helper automatically reads each file, subtracts the Fermi level, and infers labels such as `Pb(p)` for legend entries. If the PDOS file contains two columns after the energy (spin up / spin down), the routine mirrors the down component when plotting vertical orientation.
- `orientation="horizontal"` swaps axes, helpful for stacking with band plots.
- Use `read_pdos(filepath)` if you need the raw `(energy, pdos_array)` pair without plotting.

### 4.2 `OrbitalCenter`

Find the energy positions where a specific orbital PDOS is minimal or maximal, optionally within a sub-window.

```python
fermi_energy = 3.5926  # must exist in the namespace before instantiating!
oc = OrbitalCenter("PDOS/(Pb)_(p).dat")
emin, emax = oc.get_extreme_positions(x_min=0.0, x_max=4.0)
```

> **Important**: The class reads the module-level variable `fermi_energy`. Define it before creating an `OrbitalCenter` instance, or refactor the class to pass the value explicitly.

---

## 5. Charge Density Post-Processing

### 5.1 `ChargeDensitySlice`

Visualise and compare real-space charge densities imported from QE `.pp` files (converted to plain text with a 50-line header).

```python
cd = ChargeDensitySlice(
    file_path="pure.pp.txt",  # plain text data
    x1=120, x2=120, x3=120,   # FFT grid dimensions
    a=11.95, b=11.95, c=11.95 # lattice vectors in Å
)
```

Key methods:

- `plot_slice(hkl='001', position=None, interpolation='gaussian')`: Draws a 2D slice perpendicular to the specified Miller index. `position` is a real-space distance along the lattice vector; omitted value defaults to the mid-plane.
- `plot_line_profile(axis='z', ix=None, iy=None, iz=None)`: Extracts and plots a 1D line cut along the chosen axis, fixing the other coordinates (defaults to midpoints).
- `get_line_profile(axis='z', index1=None, index2=None)`: Returns `(coordinates, profile)` for scripting without plotting.
- `plot_charge_density_difference(...)`: Class-level helper that loads two files, aligns the same line profile, and plots the difference. Call it as `ChargeDensitySlice.plot_charge_density_difference(file1=..., file2=..., ...)`.

All charge densities are converted from Bohr units to Å⁻³ during initialisation.

---

## 6. Geometry Utilities

### 6.1 `QeInputFile`

Parse one or multiple QE input structures (`*.in`) and compute interatomic statistics using ASE.

- Instantiate with either a single file path or a list: `qe = QeInputFile(["scf.in", "relax.in"])`.
- `print_average_distances(symbol1, symbol2)` prints the average distance between every pair of atoms matching the provided symbols across all files.
- `get_all_distances(atoms, symbol1, symbol2)` and `get_atom_indices(atoms, symbol)` are exposed for finer-grained scripting.

### 6.2 `Minimaout` and `Minimain`

Extract atomic positions from QE output (`espresso-out`) or input (`espresso-in`) files respectively.

- `get_atomic_positions(exclude_element='O')` returns tuples of Cartesian and fractional coordinates after removing unwanted elements (defaults to dropping oxygen).
- `find_highest_z_coordinates(exclude_element='O')` identifies the highest and second-highest atoms along z, and reports their Cartesian/Fractional coordinates alongside the vertical separation. Useful for surface adsorption analyses.

Both classes rely on ASE’s parsers; ensure the corresponding `ase.io.read` format support (`espresso-out`, `espresso-in`) is available.

---

## 7. Energetics

### 7.1 `TotalEnergy`

Scan a QE SCF output log and expose the last reported total energy.

```python
te = TotalEnergy("scf.out")
te.extract()
te.last_TotalEnergy()
```

- `.extract()` updates `last_total_energy_ry` when encountering lines containing `"!    total energy"`.
- `.last_TotalEnergy()` prints the stored value in Ry and eV (requires a previous call to `.extract()`).

### 7.2 `SurfaceEnergyCalculator`

Convert slab total energies into surface energies.

```python
sec = SurfaceEnergyCalculator()
surface_E = sec.surf(data=[-358.412, -358.387], bulk=-179.201, N=2, area=[40.12, 40.12])
```

- Parameters: `data` (array of slab energies, Ry), `bulk` (bulk energy per formula unit, Ry), `N` (number of bulk units in the slab), `area` (surface area values, Å²).
- Returns surface energies in J/m² (`ry_to_jm2` factor defaults to 217.9863).

### 7.3 `AdsorptionCalculator`

Estimate adsorption energies and desorption times for adsorbates.

```python
ads = AdsorptionCalculator(Eadsorbate=-20.321)  # Ry
E_ads = ads.Eadsorption(Esys=-380.512, Eadsorbant=-358.941)
tau = ads.desorption_time(t=1e-6, Esys=-380.512, Eadsorbant=-358.941)
```

- `Eadsorption(Esys, Eadsorbant)` returns `(Esys - (Eadsorbant + Eadsorbate)) * 13.6057` in eV.
- `desorption_time(t, Esys, Eadsorbant)` applies an Arrhenius-like factor `t * exp(-E_ads / k_BT)` with `kBT` defaulting to 0.0258519 eV (~300 K).

---

## 8. Putting It Together

Typical post-processing pipeline for a QE project:

1. **Electronic structure**: Use `BandStructure` to load `.bands.gnu`, plot the band diagram, and annotate important bands with `plot_band_number`. Derive transport metrics with `EffectiveMassAnalyzer`.
2. **Density of states**: Instantiate `PDOSPlotter` to overlay orbital contributions and capture images for manuscripts. `OrbitalCenter` helps pinpoint orbital centroids.
3. **Real-space insight**: For surfaces or defects, feed charge-density grids to `ChargeDensitySlice`, check planar slices, and compare pristine vs. doped structures with the difference helper.
4. **Structural context**: `QeInputFile`, `Minimaout`, and `Minimain` give quick access to interatomic distances and atomic heights without leaving Python.
5. **Energetics**: Summarise SCF outputs with `TotalEnergy`, evaluate surface energies with `SurfaceEnergyCalculator`, and study adsorption thermodynamics via `AdsorptionCalculator`.

By composing these objects you can automate full QE post-processing routines directly in scripts or Jupyter notebooks tracked in your GitHub repository.

---

## 9. Tips & Extensibility

- **Consistent units**: Energies are handled in Ry internally where appropriate but reported in eV for plots; always confirm the expected unit before combining results.
- **Matplotlib styling**: The module sets `Times New Roman` globally. Adjust `plt.rcParams` before importing if you need different fonts for publication.
- **Missing Python interpreter**: If your environment lacks a `python` binary, set up a virtual environment (e.g. `conda create -n qe python=3.10`) and install the dependencies above.
- **Custom pipelines**: Because each class exposes Numpy arrays, you can stack additional analysis layers (e.g. fitting, machine learning, custom plots) on top.

For questions or contributions, open an issue or pull request in the repository and reference this guide to keep the documentation up to date.
