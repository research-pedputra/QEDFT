"""High-level API wrappers for the repository tools."""

from __future__ import annotations

from typing import Dict, Iterable, Optional, Sequence, Tuple

from .. import (
    BandEnergy,
    CombinedPlot,
    SCFNSCFGenerator,
    PDOSGenerator,
    BandGenerator,
    BandsGenerator,
    ChargeDensityPlotter,
    qe2xyz,
    StructureVisualizer,
)
from XRD import BackSub, peaks, PlotTxtFiles
from .qe_binaries import QERunner, QERunResult
from .unfolding import BandUnfoldingAPI


class QuantumEspressoAPI:
    """Unified Quantum ESPRESSO API.

    This class groups input generation and binary execution in one interface.
    """

    def __init__(
        self,
        qe_bin_dir: Optional[str] = None,
        mpiexec: Optional[str] = None,
        mpi_nprocs: Optional[int] = None,
    ) -> None:
        self.runner = QERunner(qe_bin_dir=qe_bin_dir, mpiexec=mpiexec, mpi_nprocs=mpi_nprocs)

    @staticmethod
    def supported_binaries() -> Dict[str, str]:
        return QERunner.supported_binaries()

    @staticmethod
    def supported_groups():
        return QERunner.supported_groups()

    def discover_available_binaries(self) -> Dict[str, str]:
        return self.runner.discover_available_binaries()

    @staticmethod
    def scf_or_nscf(cif_file_path: str, output_dir: str, calculation_type: str = "scf") -> SCFNSCFGenerator:
        return SCFNSCFGenerator(cif_file_path, output_dir, calculation_type)

    @staticmethod
    def pdos_input(cif_file_path: str, output_dir: str) -> PDOSGenerator:
        return PDOSGenerator(cif_file_path, output_dir)

    @staticmethod
    def band_input(cif_file_path: str, output_dir: str) -> BandGenerator:
        return BandGenerator(cif_file_path, output_dir)

    @staticmethod
    def bands_post_input(cif_file_path: str, output_dir: str) -> BandsGenerator:
        return BandsGenerator(cif_file_path, output_dir)

    def run(
        self,
        binary: str,
        input_file: Optional[str] = None,
        output_file: Optional[str] = None,
        npool: Optional[int] = None,
        extra_args: Optional[Sequence[str]] = None,
        cwd: Optional[str] = None,
        timeout: Optional[int] = None,
        dry_run: bool = False,
    ) -> QERunResult:
        return self.runner.run(
            binary=binary,
            input_file=input_file,
            output_file=output_file,
            npool=npool,
            extra_args=extra_args,
            cwd=cwd,
            timeout=timeout,
            dry_run=dry_run,
        )

    def run_many(self, steps, cwd: Optional[str] = None):
        return self.runner.run_many(steps=steps, cwd=cwd)

    def run_group(self, group: str, input_map, output_dir: str = ".", cwd: Optional[str] = None):
        return self.runner.run_group(group=group, input_map=input_map, output_dir=output_dir, cwd=cwd)

    def unfolding(self) -> BandUnfoldingAPI:
        """Return a band-unfolding helper bound to the same runner config."""
        return BandUnfoldingAPI(runner=self.runner)

    def run_bands_unfold(self, *args, **kwargs) -> QERunResult:
        """Explicit convenience runner for ``bands_unfold.x``."""
        return self.runner.run_bands_unfold(*args, **kwargs)

    def run_python_unfold(self, *args, **kwargs):
        """Run pure-Python unfolding via :class:`BandUnfoldingAPI`."""
        return self.unfolding().run_like_bands_unfold(*args, **kwargs)


class BandAnalysis:
    """Band structure and combined band+PDOS analysis helpers."""

    def __init__(self) -> None:
        self.engine = BandEnergy()

    def load(self, data, fermi_energy: float, kpoints: Sequence[Tuple[str, float]]) -> "BandAnalysis":
        self.engine.get_bandstructure(data, fermi_energy, kpoints)
        return self

    def plot(self) -> None:
        self.engine.plot_band_structure()

    def plot_band(self, band_number: int, overlay: bool = True) -> None:
        self.engine.plot_band_number(band_number, overlay=overlay)

    def plot_shifted(self, start_band_number: int, shift_value: float) -> None:
        self.engine.plot_shifted_band_structure(start_band_number, shift_value)

    def combine_with_pdos(
        self,
        band_data,
        fermi_energy: float,
        kpoints: Sequence[Tuple[str, float]],
        pdos_files: Dict[str, Tuple[str, str]],
        pdos_fermi_energy: float,
        band_xlim: Optional[Tuple[float, float]] = None,
        band_ylim: Optional[Tuple[float, float]] = None,
        pdos_xlim: Optional[Tuple[float, float]] = None,
        pdos_ylim: Optional[Tuple[float, float]] = None,
    ) -> None:
        CombinedPlot(self.engine).plot_combined(
            band_data=band_data,
            fermi_energy=fermi_energy,
            kpoints=kpoints,
            pdos_files=pdos_files,
            pdos_fermi_energy=pdos_fermi_energy,
            band_xlim=band_xlim,
            band_ylim=band_ylim,
            pdos_xlim=pdos_xlim,
            pdos_ylim=pdos_ylim,
        )


class PDOSAnalysis:
    """PDOS parsing and plotting via BandEnergy helpers."""

    def __init__(self) -> None:
        self.engine = BandEnergy()

    def read(self, files: Dict[str, Tuple[str, str]], fermi_energy: float):
        return self.engine.get_aopdos(files, fermi_energy)

    def plot(self, aopdos_data, xlim=None, ylim=None) -> None:
        self.engine.plot_aopdos(aopdos_data, xlim=xlim, ylim=ylim)


class QEInputFactory:
    """Factory for generating Quantum ESPRESSO input files from CIF structures."""

    @staticmethod
    def scf_or_nscf(cif_file_path: str, output_dir: str, calculation_type: str = "scf") -> SCFNSCFGenerator:
        return SCFNSCFGenerator(cif_file_path, output_dir, calculation_type)

    @staticmethod
    def pdos(cif_file_path: str, output_dir: str) -> PDOSGenerator:
        return PDOSGenerator(cif_file_path, output_dir)

    @staticmethod
    def band(cif_file_path: str, output_dir: str) -> BandGenerator:
        return BandGenerator(cif_file_path, output_dir)

    @staticmethod
    def bands_post(cif_file_path: str, output_dir: str) -> BandsGenerator:
        return BandsGenerator(cif_file_path, output_dir)


class ChargeDensityAnalysis:
    """Charge density difference workflow wrapper."""

    def __init__(self, cube_files: Sequence[str], dimensions: Tuple[int, int, int], scale_factor: float):
        self.engine = ChargeDensityPlotter(cube_files, dimensions, scale_factor)

    def run(self, output_file: str) -> None:
        self.engine.load_cube_data()
        self.engine.compute_charge_density_difference()
        self.engine.reshape_and_sum_charge_density()
        self.engine.plot_charge_density_difference(output_file)


class StructureTools:
    """Structure conversion and visualization helper methods."""

    @staticmethod
    def qe_to_xyz(input_file: str, output_file: str) -> None:
        qe2xyz(input_file).convert(output_file)

    @staticmethod
    def html_view(file_path: str):
        return StructureVisualizer(file_path).Structure()


class XRDAnalysis:
    """XRD background subtraction, peak finding, and batch plotting."""

    @staticmethod
    def subtract_background(x: Iterable[float], y: Iterable[float], tol: float = 1):
        return BackSub(x, y).backsub(tol=tol)

    @staticmethod
    def find_peaks(
        x: Iterable[float],
        y: Iterable[float],
        height=None,
        distance=None,
        prominence=None,
        width=None,
    ):
        return peaks(x, y).find_peaks(
            height=height,
            distance=distance,
            prominence=prominence,
            width=width,
        )

    @staticmethod
    def plot_folder(folder_path: str) -> None:
        PlotTxtFiles(folder_path).plot_files()
