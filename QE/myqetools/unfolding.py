"""Pure-Python band unfolding helpers (bands_unfold-like workflow).

This module provides a QE-independent unfolding approximation implemented in
Python and suitable for notebooks. A legacy QE binary hook is kept optional.
"""

from __future__ import annotations

from math import exp, sqrt, pi
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .qe_binaries import QERunner, QERunResult


class BandUnfoldingAPI:
    """Band unfolding interface.

    Main path: pure Python unfolding via :meth:`run_like_bands_unfold`.
    Optional path: execute QE `bands_unfold.x` via :meth:`run_qe_binary`.
    """

    def __init__(self, runner: Optional[QERunner] = None, **runner_kwargs) -> None:
        self.runner = runner or QERunner(**runner_kwargs)

    @staticmethod
    def write_input_file(
        file_path: str,
        namelist: Dict[str, object],
        card_lines: Optional[Sequence[str]] = None,
    ) -> str:
        """Write a simple QE-style `bands_unfold.x` input file (optional path)."""
        target = Path(file_path)
        target.parent.mkdir(parents=True, exist_ok=True)

        def fmt(value: object) -> str:
            if isinstance(value, str):
                return f"'{value}'"
            if isinstance(value, bool):
                return ".true." if value else ".false."
            return str(value)

        with target.open("w", encoding="utf-8") as f:
            f.write("&inputpp\n")
            for key, value in namelist.items():
                f.write(f"  {key} = {fmt(value)}\n")
            f.write("/\n")
            if card_lines:
                for line in card_lines:
                    f.write(f"{line}\n")

        return str(target)

    def run_qe_binary(
        self,
        input_file: str,
        output_file: Optional[str] = None,
        cwd: Optional[str] = None,
        timeout: Optional[int] = None,
        dry_run: bool = False,
    ) -> QERunResult:
        """Optional QE-native execution via `bands_unfold.x`."""
        return self.runner.run(
            binary="bands_unfold",
            input_file=input_file,
            output_file=output_file,
            cwd=cwd,
            timeout=timeout,
            dry_run=dry_run,
        )

    # backward-compatible alias
    run = run_qe_binary

    @staticmethod
    def _fold_kpoint(k: Sequence[float], supercell_matrix: Sequence[Sequence[float]]) -> Tuple[float, float, float]:
        # simple inverse-transpose mapping for diagonal-like transforms
        # robust enough for common supercell transforms used in scripts
        s = supercell_matrix
        # assume 3x3 and invert numerically (small system, explicit formula avoided)
        import numpy as np

        s_inv_t = np.linalg.inv(np.asarray(s, dtype=float)).T
        k_prim = (np.asarray(k, dtype=float) @ s_inv_t) % 1.0
        return float(k_prim[0]), float(k_prim[1]), float(k_prim[2])

    @staticmethod
    def run_like_bands_unfold(
        kpoints_supercell,
        energies_supercell,
        supercell_matrix,
        weights=None,
        k_path_points: int = 400,
        energy_grid_points: int = 600,
        sigma_k: float = 0.01,
        sigma_e: float = 0.05,
        emin: Optional[float] = None,
        emax: Optional[float] = None,
    ):
        """Pure-Python `bands_unfold`-like unfolding.

        Returns a dict with flattened points and a gridded intensity map.
        """
        import numpy as np

        k_sc = np.asarray(kpoints_supercell, dtype=float)
        e_sc = np.asarray(energies_supercell, dtype=float)
        if e_sc.ndim == 1:
            e_sc = e_sc[:, None]

        if weights is None:
            w_sc = np.ones_like(e_sc)
        else:
            w_sc = np.asarray(weights, dtype=float)
            if w_sc.shape != e_sc.shape:
                raise ValueError("weights must have same shape as energies_supercell")

        s_inv_t = np.linalg.inv(np.asarray(supercell_matrix, dtype=float)).T
        k_prim = (k_sc @ s_inv_t) % 1.0

        # order along pseudo path
        order = np.lexsort((k_prim[:, 2], k_prim[:, 1], k_prim[:, 0]))
        k_sorted = k_prim[order]
        e_sorted = e_sc[order]
        w_sorted = w_sc[order]

        dk = np.linalg.norm(np.diff(k_sorted, axis=0), axis=1)
        k_line_raw = np.concatenate([[0.0], np.cumsum(dk)])
        if k_line_raw[-1] <= 1e-15:
            k_line_raw[-1] = 1.0
        k_line_raw /= k_line_raw[-1]

        unfolded_k = np.repeat(k_line_raw, e_sorted.shape[1])
        unfolded_e = e_sorted.reshape(-1)
        unfolded_w = w_sorted.reshape(-1)

        if emin is None:
            emin = float(np.nanmin(unfolded_e))
        if emax is None:
            emax = float(np.nanmax(unfolded_e))

        k_line = np.linspace(0.0, 1.0, k_path_points)
        e_grid = np.linspace(emin, emax, energy_grid_points)

        k_diff = k_line[None, :] - unfolded_k[:, None]
        k_kernel = np.exp(-0.5 * (k_diff / sigma_k) ** 2)
        e_diff = e_grid[:, None] - unfolded_e[None, :]
        e_kernel = np.exp(-0.5 * (e_diff / sigma_e) ** 2)

        intensity = (e_kernel * unfolded_w[None, :]) @ k_kernel
        norm = 2.0 * pi * sigma_k * sigma_e
        if norm > 0:
            intensity /= norm

        return {
            "k_line": k_line,
            "e_grid": e_grid,
            "intensity": intensity,
            "unfolded_k": unfolded_k,
            "unfolded_e": unfolded_e,
            "unfolded_w": unfolded_w,
        }

    @staticmethod
    def save_unfolded_python(unfolded: Dict[str, object], output_file: str) -> str:
        import numpy as np

        out = Path(output_file)
        out.parent.mkdir(parents=True, exist_ok=True)
        data = np.column_stack(
            [
                np.asarray(unfolded["unfolded_k"]),
                np.asarray(unfolded["unfolded_e"]),
                np.asarray(unfolded["unfolded_w"]),
            ]
        )
        np.savetxt(out, data, header="k E weight")
        return str(out)

    @staticmethod
    def load_unfolded_data(
        data_file: str,
        k_col: int = 0,
        e_col: int = 1,
        w_col: int = 2,
        skiprows: int = 0,
    ) -> Tuple[object, object, object]:
        import numpy as np

        data = np.loadtxt(data_file, skiprows=skiprows)
        return data[:, k_col], data[:, e_col], data[:, w_col]

    @staticmethod
    def plot_notebook(
        data_file: str,
        fermi_energy: float = 0.0,
        k_col: int = 0,
        e_col: int = 1,
        w_col: int = 2,
        skiprows: int = 0,
        cmap: str = "viridis",
        marker_size: float = 10.0,
        alpha: float = 0.8,
        figsize: Tuple[float, float] = (7, 5),
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        show_colorbar: bool = True,
        title: str = "Unfolded band structure",
    ):
        import matplotlib.pyplot as plt

        k, e, w = BandUnfoldingAPI.load_unfolded_data(
            data_file=data_file,
            k_col=k_col,
            e_col=e_col,
            w_col=w_col,
            skiprows=skiprows,
        )

        fig, ax = plt.subplots(figsize=figsize)
        sc = ax.scatter(
            k,
            e - fermi_energy,
            c=w,
            s=marker_size,
            cmap=cmap,
            alpha=alpha,
            vmin=vmin,
            vmax=vmax,
            edgecolors="none",
        )
        ax.set_xlabel("k-path")
        ax.set_ylabel("E - E$_f$ (eV)")
        ax.set_title(title)
        ax.axhline(0.0, color="gray", linestyle="--", linewidth=1)

        if show_colorbar:
            cbar = fig.colorbar(sc, ax=ax)
            cbar.set_label("Unfolding weight")

        return fig, ax, sc

    @staticmethod
    def plot_intensity_notebook(
        unfolded: Dict[str, object],
        fermi_energy: float = 0.0,
        figsize: Tuple[float, float] = (7, 5),
        cmap: str = "magma",
        title: str = "Unfolded spectral intensity",
    ):
        import matplotlib.pyplot as plt
        import numpy as np

        k_line = np.asarray(unfolded["k_line"])
        e_grid = np.asarray(unfolded["e_grid"]) - fermi_energy
        intensity = np.asarray(unfolded["intensity"])

        fig, ax = plt.subplots(figsize=figsize)
        im = ax.imshow(
            intensity,
            origin="lower",
            aspect="auto",
            extent=[k_line.min(), k_line.max(), e_grid.min(), e_grid.max()],
            cmap=cmap,
        )
        ax.set_xlabel("k-path")
        ax.set_ylabel("E - E$_f$ (eV)")
        ax.set_title(title)
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Spectral intensity")

        return fig, ax, im
