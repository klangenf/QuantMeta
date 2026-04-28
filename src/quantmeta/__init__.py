"""QuantMeta: Absolute abundance quantification using synthetic DNA standards."""

from .detection_threshold.run_sh import main as detection_threshold
from .standard_curve.run_sh import main as standard_curve
from .quant_targets.run_sh import main as quant_targets

__all__ = [
    "detection_threshold",
    "standard_curve",
    "quant_targets"
]
__version__ = "2.0.0"