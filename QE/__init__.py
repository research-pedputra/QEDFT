"""Quantum ESPRESSO utilities package.

Exports are loaded lazily so lightweight subpackages can be imported even when
optional scientific dependencies are not installed.
"""

__all__ = [
    "BandEnergy",
    "CombinedPlot",
    "SCFNSCFGenerator",
    "PDOSGenerator",
    "BandGenerator",
    "BandsGenerator",
    "ChargeDensityPlotter",
    "shirley",
    "linear_bg",
    "qe2xyz",
    "StructureVisualizer",
]


def __getattr__(name):
    if name in {"BandEnergy", "CombinedPlot"}:
        from .bandenergy import BandEnergy, CombinedPlot

        return {"BandEnergy": BandEnergy, "CombinedPlot": CombinedPlot}[name]
    if name in {"SCFNSCFGenerator", "PDOSGenerator", "BandGenerator", "BandsGenerator"}:
        from .qeinputmaker import SCFNSCFGenerator, PDOSGenerator, BandGenerator, BandsGenerator

        return {
            "SCFNSCFGenerator": SCFNSCFGenerator,
            "PDOSGenerator": PDOSGenerator,
            "BandGenerator": BandGenerator,
            "BandsGenerator": BandsGenerator,
        }[name]
    if name == "ChargeDensityPlotter":
        from .densityplot import ChargeDensityPlotter

        return ChargeDensityPlotter
    if name in {"shirley", "linear_bg"}:
        from .background import shirley, linear_bg

        return {"shirley": shirley, "linear_bg": linear_bg}[name]
    if name in {"qe2xyz", "StructureVisualizer"}:
        from .strucvis import qe2xyz, StructureVisualizer

        return {"qe2xyz": qe2xyz, "StructureVisualizer": StructureVisualizer}[name]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
