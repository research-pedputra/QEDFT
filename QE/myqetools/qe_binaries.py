"""Quantum ESPRESSO binary runner API.

This module provides a broad, uniform wrapper around Quantum ESPRESSO executables,
including PWscf, PHonon, NEB, TDDFPT, EPW, HP, XSPECTRA, and converter tools.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import os
import shutil
import subprocess
from typing import Dict, List, Optional, Sequence


QE_BINARY_ALIASES: Dict[str, str] = {
    # Core PWscf and post-processing
    "pw": "pw.x",
    "cp": "cp.x",
    "pp": "pp.x",
    "dos": "dos.x",
    "bands": "bands.x",
    "projwfc": "projwfc.x",
    "average": "average.x",
    "plan_avg": "plan_avg.x",
    "plotband": "plotband.x",
    "plotproj": "plotproj.x",
    "plotrho": "plotrho.x",
    "sumpdos": "sumpdos.x",
    "open_grid": "open_grid.x",
    # Phonons
    "ph": "ph.x",
    "q2r": "q2r.x",
    "matdyn": "matdyn.x",
    "dynmat": "dynmat.x",
    "lambda": "lambda.x",
    # Nudged Elastic Band
    "neb": "neb.x",
    "path_interpolation": "path_interpolation.x",
    # Interfaces / converters
    "pw2wannier90": "pw2wannier90.x",
    "wannier_ham": "wannier_ham.x",
    "pw2bgw": "pw2bgw.x",
    "pw2gw": "pw2gw.x",
    "pw2casino": "pw2casino.x",
    "pw2critic": "pw2critic.x",
    "bands_unfold": "bands_unfold.x",
    "kcw": "kcw.x",
    "kcw_pp": "kcw_pp.x",
    # TDDFPT
    "turbo_lanczos": "turbo_lanczos.x",
    "turbo_davidson": "turbo_davidson.x",
    "turbo_eels": "turbo_eels.x",
    "turbo_spectrum": "turbo_spectrum.x",
    # GWW / XSPECTRA
    "gww_fit": "gww_fit.x",
    "gww_pp": "gww_pp.x",
    "epsilon": "epsilon.x",
    "xspectra": "xspectra.x",
    # EPW / HP
    "epw": "epw.x",
    "hp": "hp.x",
    # Atomic / pseudo tools
    "ld1": "ld1.x",
    "atomic": "atomic.x",
    "cppp": "cppp.x",
    "upfconv": "upfconv.x",
    "virtual_v2": "virtual_v2.x",
    # Legacy / ecosystem tools often distributed with QE stacks
    "manycp": "manycp.x",
    "dist": "dist.x",
    "ev": "ev.x",
    "initial_state": "initial_state.x",
    "fermi_velocity": "fermi_velocity.x",
}

QE_BINARY_GROUPS: Dict[str, Sequence[str]] = {
    "pw": ["pw", "cp", "pp", "dos", "bands", "projwfc", "average", "plan_avg", "plotband", "plotproj", "plotrho", "sumpdos", "open_grid"],
    "phonon": ["ph", "q2r", "matdyn", "dynmat", "lambda"],
    "neb": ["neb", "path_interpolation"],
    "interfaces": ["pw2wannier90", "wannier_ham", "pw2bgw", "pw2gw", "pw2casino", "pw2critic", "bands_unfold", "kcw", "kcw_pp"],
    "tddfpt": ["turbo_lanczos", "turbo_davidson", "turbo_eels", "turbo_spectrum"],
    "spectroscopy": ["gww_fit", "gww_pp", "epsilon", "xspectra"],
    "transport": ["epw", "hp"],
    "pseudo": ["ld1", "atomic", "cppp", "upfconv", "virtual_v2"],
}


@dataclass
class QERunResult:
    """Result information from a QE binary execution."""

    command: List[str]
    returncode: int
    stdout: str = ""
    stdout_file: Optional[str] = None
    stderr: str = ""


@dataclass
class QERunner:
    """Unified runner for Quantum ESPRESSO executables."""

    qe_bin_dir: Optional[str] = None
    mpiexec: Optional[str] = None
    mpi_nprocs: Optional[int] = None
    extra_env: Dict[str, str] = field(default_factory=dict)

    @staticmethod
    def supported_binaries() -> Dict[str, str]:
        return dict(QE_BINARY_ALIASES)

    @staticmethod
    def supported_groups() -> Dict[str, Sequence[str]]:
        return dict(QE_BINARY_GROUPS)

    def discover_available_binaries(self) -> Dict[str, str]:
        """Return aliases that are currently discoverable on this machine."""
        found: Dict[str, str] = {}
        for alias, executable in QE_BINARY_ALIASES.items():
            resolved = None
            if self.qe_bin_dir:
                candidate = Path(self.qe_bin_dir) / executable
                if candidate.exists() and os.access(candidate, os.X_OK):
                    resolved = str(candidate)
            if resolved is None:
                resolved = shutil.which(executable)
            if resolved:
                found[alias] = resolved
        return found

    def _resolve_binary(self, binary: str) -> str:
        executable = QE_BINARY_ALIASES.get(binary, binary)
        if self.qe_bin_dir:
            path = Path(self.qe_bin_dir) / executable
            if path.exists():
                return str(path)
        found = shutil.which(executable)
        if found:
            return found
        raise FileNotFoundError(
            f"Could not find QE binary '{executable}'. Provide qe_bin_dir or add it to PATH."
        )

    def _build_command(
        self,
        binary: str,
        input_file: Optional[str] = None,
        npool: Optional[int] = None,
        extra_args: Optional[Sequence[str]] = None,
    ) -> List[str]:
        cmd: List[str] = []
        if self.mpiexec:
            cmd.append(self.mpiexec)
            if self.mpi_nprocs:
                cmd.extend(["-np", str(self.mpi_nprocs)])

        cmd.append(self._resolve_binary(binary))

        if npool is not None:
            cmd.extend(["-nk", str(npool)])
        if input_file is not None:
            cmd.extend(["-in", input_file])
        if extra_args:
            cmd.extend(list(extra_args))
        return cmd

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
        command = self._build_command(binary=binary, input_file=input_file, npool=npool, extra_args=extra_args)
        if dry_run:
            return QERunResult(command=command, returncode=0, stdout_file=output_file)

        env = os.environ.copy()
        env.update(self.extra_env)

        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with output_path.open("w", encoding="utf-8") as out:
                proc = subprocess.run(
                    command,
                    cwd=cwd,
                    env=env,
                    stdout=out,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=timeout,
                    check=False,
                )
                stdout_text = ""
        else:
            proc = subprocess.run(
                command,
                cwd=cwd,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout,
                check=False,
            )
            stdout_text = proc.stdout

        if proc.returncode != 0:
            raise RuntimeError(
                f"QE run failed (exit={proc.returncode}) for command: {' '.join(command)}\n{proc.stderr}"
            )

        return QERunResult(
            command=command,
            returncode=proc.returncode,
            stdout=stdout_text,
            stdout_file=output_file,
            stderr=proc.stderr,
        )

    def run_many(
        self,
        steps: Sequence[Dict[str, object]],
        cwd: Optional[str] = None,
    ) -> List[QERunResult]:
        """Execute a sequence of QE steps.

        Each step supports keys from ``run``: binary, input_file, output_file,
        npool, extra_args, timeout, and dry_run.
        """
        results: List[QERunResult] = []
        for step in steps:
            result = self.run(
                binary=str(step["binary"]),
                input_file=step.get("input_file"),
                output_file=step.get("output_file"),
                npool=step.get("npool"),
                extra_args=step.get("extra_args"),
                cwd=step.get("cwd", cwd),
                timeout=step.get("timeout"),
                dry_run=bool(step.get("dry_run", False)),
            )
            results.append(result)
        return results
