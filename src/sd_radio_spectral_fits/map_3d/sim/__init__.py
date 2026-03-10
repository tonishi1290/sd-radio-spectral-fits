"""
OTF Simulation Tools
Mock data generators for testing OTF gridding, basket-weaving, and pipelines.
"""

from .otf_simulator import (
    GaussianSource, 
    SkyModel, 
    create_complex_sky_model, 
    simulate_otf_observation
)

__all__ = [
    "GaussianSource",
    "SkyModel",
    "create_complex_sky_model",
    "simulate_otf_observation"
]
