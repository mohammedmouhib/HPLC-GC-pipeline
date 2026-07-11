"""HPLC (MOCCA2) data-processing pipeline.

A reusable, config-driven tool that turns Agilent Chemstation .D injection
folders into integrated peak tables, compound concentrations, replicate
statistics, and publication-quality plots. Point it at any experiment folder;
all settings live in that folder's hplc_config.yaml.
"""

from .config import Config, load_config

__version__ = "0.1.0"
__all__ = ["Config", "load_config", "__version__"]
