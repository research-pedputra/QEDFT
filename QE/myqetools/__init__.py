"""Public API for the QE myqetools wrapper.

Notes
-----
This module uses lazy attribute loading so lightweight utilities like
``QERunner`` can be imported without requiring heavy scientific dependencies.
"""

from .qe_binaries import QERunner, QERunResult, QE_BINARY_ALIASES, QE_BINARY_GROUPS

__all__ = [
    "QuantumEspressoAPI",
    "BandAnalysis",
    "PDOSAnalysis",
    "QEInputFactory",
    "ChargeDensityAnalysis",
    "StructureTools",
    "XRDAnalysis",
    "QERunner",
    "QERunResult",
    "QE_BINARY_ALIASES",
    "QE_BINARY_GROUPS",
]


def __getattr__(name):
    if name in {
        "QuantumEspressoAPI",
        "BandAnalysis",
        "PDOSAnalysis",
        "QEInputFactory",
        "ChargeDensityAnalysis",
        "StructureTools",
        "XRDAnalysis",
    }:
        from .api import (
            QuantumEspressoAPI,
            BandAnalysis,
            PDOSAnalysis,
            QEInputFactory,
            ChargeDensityAnalysis,
            StructureTools,
            XRDAnalysis,
        )

        mapping = {
            "QuantumEspressoAPI": QuantumEspressoAPI,
            "BandAnalysis": BandAnalysis,
            "PDOSAnalysis": PDOSAnalysis,
            "QEInputFactory": QEInputFactory,
            "ChargeDensityAnalysis": ChargeDensityAnalysis,
            "StructureTools": StructureTools,
            "XRDAnalysis": XRDAnalysis,
        }
        return mapping[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
