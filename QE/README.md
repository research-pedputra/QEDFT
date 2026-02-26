# QEDFT Toolkit (QE API)

This repository exposes a broad Quantum ESPRESSO API under `QE/myqetools`, intended to cover most practical QE capabilities while preserving the original utility modules.

## Location

- API package: `QE/myqetools/`
- Core QE utilities: `QE/*.py`
- XRD helpers consumed by API: `XRD/`

## Dependencies

- For binary execution only (`QERunner`): standard library only.
- For full analysis API (`QuantumEspressoAPI`, plotting helpers): install `numpy`, `matplotlib`, `scipy`, `ase`, and optional `PyQt5`.

## Main import

```python
from QE.myqetools import QERunner
# optional heavy API:
from QE.myqetools import QuantumEspressoAPI
```

## Coverage in this API

The wrapper includes:

1. Input-generation wrappers (`scf`, `nscf`, `bands`, `pdos`).
2. Binary execution for many QE tools through aliases.
3. Grouped capability metadata (`supported_groups`).
4. Runtime binary discovery (`discover_available_binaries`).
5. Multi-step execution (`run_many`) and group execution (`run_group`).

## Quick usage

### 1) Generate inputs + run a binary

```python
from QE.myqetools import QuantumEspressoAPI

qe = QuantumEspressoAPI(qe_bin_dir="/opt/qe/bin", mpiexec="mpirun", mpi_nprocs=16)
qe.scf_or_nscf("sample.cif", "./inputs", calculation_type="scf").get_input()
qe.run("pw", input_file="./inputs/scf_sample.in", output_file="./outputs/scf.out")
```

### 2) Inspect groups and discover available binaries

```python
from QE.myqetools import QERunner

runner = QERunner(qe_bin_dir="/opt/qe/bin")
print(runner.supported_groups())
print(runner.discover_available_binaries())
```

### 3) Run multiple steps

```python
from QE.myqetools import QuantumEspressoAPI

qe = QuantumEspressoAPI(qe_bin_dir="/opt/qe/bin")
steps = [
    {"binary": "pw", "input_file": "scf.in", "output_file": "scf.out"},
    {"binary": "pw", "input_file": "nscf.in", "output_file": "nscf.out"},
    {"binary": "bands", "input_file": "bands_pp.in", "output_file": "bands_pp.out"},
]
qe.run_many(steps)
```

### 4) Run a capability group with alias->input mapping

```python
from QE.myqetools import QERunner

runner = QERunner(qe_bin_dir="/opt/qe/bin")
runner.run_group(
    group="phonon",
    input_map={
        "ph": "ph.in",
        "q2r": "q2r.in",
        "matdyn": "matdyn.in",
    },
    output_dir="./phonon_out",
)
```

## QE capability groups exposed

- `pw`
- `phonon`
- `neb`
- `interfaces`
- `tddfpt`
- `spectroscopy`
- `transport`
- `pseudo`

The aliases can be read directly from `QE_BINARY_ALIASES`, and non-listed binaries can still be run by explicit name:

```python
qe.run("my_custom_qe_binary.x", input_file="input.in", output_file="out.log")
```


## Band unfolding (Python-first)

Use the Python unfolding path (bands_unfold-like) by default:

```python
import numpy as np
from QE.myqetools import QuantumEspressoAPI

qe = QuantumEspressoAPI()

# synthetic or parsed supercell data
k_sc = np.random.rand(200, 3)
E_sc = np.random.normal(loc=0.0, scale=2.0, size=(200, 8))
S = np.diag([2, 2, 1])

unfolded = qe.run_python_unfold(
    kpoints_supercell=k_sc,
    energies_supercell=E_sc,
    supercell_matrix=S,
    sigma_k=0.02,
    sigma_e=0.08,
)

fig, ax, im = qe.unfolding().plot_intensity_notebook(unfolded, fermi_energy=0.0)
qe.unfolding().save_unfolded_python(unfolded, "./out/unfolded_python.dat")
```

Optional legacy path (if desired):

```python
qe.run_bands_unfold(input_file="./in/unfold.in", output_file="./out/unfold.log")
```
