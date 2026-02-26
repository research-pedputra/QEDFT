"""Quantum ESPRESSO binary runner API.

This module provides a uniform wrapper around Quantum ESPRESSO executables.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import os
import shutil
import subprocess
from typing import Dict, List, Optional, Sequence


QE_BINARY_ALIASES: Dict[str, str] = {
    # pw suite
    "pw": "pw.x",
    "cp": "cp.x",
    "pp": "pp.x",
    "dos": "dos.x",
    "bands": "bands.x",
    "projwfc": "projwfc.x",
    "ph": "ph.x",
    "q2r": "q2r.x",
    "matdyn": "matdyn.x",
    "dynmat": "dynmat.x",
    "neb": "neb.x",
    "pw2wannier90": "pw2wannier90.x",
    # converters/tools
    "pw2bgw": "pw2bgw.x",
    "pw2gw": "pw2gw.x",
    "pw2casino": "pw2casino.x",
    "pw2critic": "pw2critic.x",
    "bands_unfold": "bands_unfold.x",
    # tddfpt
    "turbo_lanczos": "turbo_lanczos.x",
    "turbo_davidson": "turbo_davidson.x",
    "turbo_eels": "turbo_eels.x",
    "turbo_spectrum": "turbo_spectrum.x",
    # gww / xspectra
    "gww_fit": "gww_fit.x",
    "epsilon": "epsilon.x",
    "xspectra": "xspectra.x",
    # EPW / hp
    "epw": "epw.x",
    "hp": "hp.x",
    # misc
    "ld1": "ld1.x",
    "atomic": "atomic.x",
    "cppp": "cppp.x",
    "upfconv": "upfconv.x",
}


@dataclass
class QERunResult:
    """Result information from a QE binary execution."""

    command: List[str]
    returncode: int
    stdout_file: Optional[str] = None
    stderr: str = ""


@dataclass
class QERunner:
    """Unified runner for Quantum ESPRESSO executables.

    Parameters
    ----------
    qe_bin_dir:
        Optional directory that contains QE binaries. If omitted, binaries are
        resolved from the system PATH.
    mpiexec:
        Optional MPI launcher command (for example: ``"mpirun"``).
    mpi_nprocs:
        Number of MPI processes when ``mpiexec`` is set.
    """

    qe_bin_dir: Optional[str] = None
    mpiexec: Optional[str] = None
    mpi_nprocs: Optional[int] = None
    extra_env: Dict[str, str] = field(default_factory=dict)

    @staticmethod
    def supported_binaries() -> Dict[str, str]:
        """Return supported shorthand->executable mappings."""
        return dict(QE_BINARY_ALIASES)

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
        """Execute a QE binary with a normalized interface.

        Parameters
        ----------
        binary:
            Either shorthand (for example ``"pw"``) or executable name
            (for example ``"pw.x"``).
        input_file:
            QE input file path passed as ``-in <file>``.
        output_file:
            Optional file path where stdout is written.
        npool:
            Optional QE k-point pool parallelism (``-nk``).
        extra_args:
            Additional command-line arguments appended to the command.
        cwd:
            Working directory for the process.
        timeout:
            Timeout in seconds.
        dry_run:
            If True, command is built and returned but not executed.
        """
        command = self._build_command(
            binary=binary,
            input_file=input_file,
            npool=npool,
            extra_args=extra_args,
        )
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

        if proc.returncode != 0:
            raise RuntimeError(
                f"QE run failed (exit={proc.returncode}) for command: {' '.join(command)}\n{proc.stderr}"
            )

        return QERunResult(
            command=command,
            returncode=proc.returncode,
            stdout_file=output_file,
            stderr=proc.stderr,
        )

    # Convenience methods for commonly used QE binaries
    def run_pw(self, *args, **kwargs) -> QERunResult:
        return self.run("pw", *args, **kwargs)

    def run_bands(self, *args, **kwargs) -> QERunResult:
        return self.run("bands", *args, **kwargs)

    def run_projwfc(self, *args, **kwargs) -> QERunResult:
        return self.run("projwfc", *args, **kwargs)

    def run_dos(self, *args, **kwargs) -> QERunResult:
        return self.run("dos", *args, **kwargs)

    def run_pp(self, *args, **kwargs) -> QERunResult:
        return self.run("pp", *args, **kwargs)

    def run_ph(self, *args, **kwargs) -> QERunResult:
        return self.run("ph", *args, **kwargs)

    def run_q2r(self, *args, **kwargs) -> QERunResult:
        return self.run("q2r", *args, **kwargs)

    def run_matdyn(self, *args, **kwargs) -> QERunResult:
        return self.run("matdyn", *args, **kwargs)

    def run_dynmat(self, *args, **kwargs) -> QERunResult:
        return self.run("dynmat", *args, **kwargs)

    def run_neb(self, *args, **kwargs) -> QERunResult:
        return self.run("neb", *args, **kwargs)
