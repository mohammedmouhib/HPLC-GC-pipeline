"""
Agilent GC-MS ``.D`` folder discovery and trace loading.

The HPLC path uses MOCCA2 to read DAD ``.D`` folders. GC-MS ``.D`` folders have a
different internal layout (``data.ms`` mass-spectral scans, an ``FID1A.ch``
channel, ``AcqData/sample_info.xml`` metadata) and MOCCA2 does not model them, so
GC uses the ``rainbow`` reader instead.

This module only *reads* data. It exposes, per injection:

* the retention-time axis (minutes),
* the monitored m/z list and the per-ion signal matrix (the SIM channels),
* the total ion chromatogram (sum over ions),
* the FID trace, if present.

Peak detection and integration live in :mod:`gc_processing`.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np


@dataclass
class GCInjection:
    """One GC-MS ``.D`` injection folder and its parsed metadata."""
    folder_path: Path
    folder_name: str          # e.g. "50_s11_A.D"
    sample_name: str          # from AcqData/sample_info.xml, else folder stem
    experiment_name: str      # run subfolder name, or "Main" for a flat layout


@dataclass
class GCTraces:
    """Loaded signal traces for one GC-MS injection.

    ``eic`` is a (n_scans x n_mz) matrix of the monitored SIM ions; column ``i``
    corresponds to ``mz[i]``. ``tic`` is the sum across ions. FID fields are
    ``None`` when the injection has no FID channel.
    """
    sample_name: str
    folder_name: str
    experiment_name: str
    time_min: np.ndarray                 # (n_scans,) retention time, minutes
    mz: np.ndarray                       # (n_mz,) monitored m/z values
    eic: np.ndarray                      # (n_scans, n_mz) per-ion signal
    tic: np.ndarray                      # (n_scans,) total ion current
    fid_time_min: Optional[np.ndarray] = None
    fid_signal: Optional[np.ndarray] = None

    def eic_for(self, target_mz: float, tol: float = 0.5) -> np.ndarray:
        """Return the extracted-ion chromatogram closest to ``target_mz``.

        Raises if no monitored ion falls within ``tol`` of the request, so a
        typo'd quantifier m/z in the compound table fails loudly rather than
        silently integrating the wrong ion.
        """
        if self.mz.size == 0:
            raise ValueError("injection has no monitored ions")
        i = int(np.argmin(np.abs(self.mz - target_mz)))
        if abs(self.mz[i] - target_mz) > tol:
            raise ValueError(
                f"m/z {target_mz} not monitored in {self.sample_name} "
                f"(SIM ions: {[float(m) for m in self.mz]})"
            )
        return self.eic[:, i]


def parse_gc_sample_name(d_folder: Path) -> str:
    """Read the 'Sample Name' field from ``AcqData/sample_info.xml``.

    Falls back to the ``.D`` folder stem when the file or field is missing.
    """
    xml_path = d_folder / "AcqData" / "sample_info.xml"
    stem = d_folder.name[:-2] if d_folder.name.endswith(".D") else d_folder.name
    if not xml_path.exists():
        return stem
    try:
        root = ET.parse(xml_path).getroot()
        for field_el in root.findall("Field"):
            name_el = field_el.find("Name")
            value_el = field_el.find("Value")
            if name_el is not None and (name_el.text or "").strip() == "Sample Name":
                text = (value_el.text or "").strip() if value_el is not None else ""
                return text or stem
    except Exception as exc:  # malformed XML, encoding issues, etc.
        print(f"  Warning: could not parse {xml_path}: {exc}")
    return stem


def _is_gc_d_folder(path: Path) -> bool:
    return path.is_dir() and path.name.endswith(".D") and (path / "data.ms").exists()


def find_gc_injections(data_root: Path) -> list[GCInjection]:
    """Find GC-MS ``.D`` folders (those containing ``data.ms``) anywhere under data_root.

    Fully recursive — works regardless of nesting depth.  The ``experiment_name``
    is set to the first-level subdirectory name below data_root.
    """
    data_root = Path(data_root)
    injections: list[GCInjection] = []

    def scan(directory: Path, batch_name: str = "") -> None:
        try:
            entries = sorted(directory.iterdir())
        except PermissionError:
            return
        for item in entries:
            if not item.is_dir():
                continue
            if _is_gc_d_folder(item):
                injections.append(
                    GCInjection(
                        folder_path=item,
                        folder_name=item.name,
                        sample_name=parse_gc_sample_name(item),
                        experiment_name=batch_name or "Main",
                    )
                )
            elif not item.name.endswith(".D") and not item.name.startswith("."):
                scan(item, batch_name or item.name)

    scan(data_root)
    injections.sort(key=lambda inj: (inj.experiment_name, inj.folder_name))
    return injections


def load_gc_traces(inj: GCInjection) -> GCTraces:
    """Load MS (SIM) and FID traces for one injection via ``rainbow``."""
    import rainbow  # imported lazily so non-GC runs don't need the dependency

    data = rainbow.read(str(inj.folder_path))

    ms = data.get_file("data.ms")
    time_min = np.asarray(ms.xlabels, dtype=float).flatten()
    mz = np.asarray(ms.ylabels, dtype=float).flatten()
    eic = np.asarray(ms.data, dtype=float)
    if eic.ndim == 1:
        eic = eic.reshape(-1, 1)
    tic = eic.sum(axis=1)

    fid_time = fid_signal = None
    fid_files = [f for f in data.datafiles if f.name.lower().endswith(".ch")]
    if fid_files:
        fid = fid_files[0]
        fid_time = np.asarray(fid.xlabels, dtype=float).flatten()
        fid_signal = np.asarray(fid.data, dtype=float).flatten()

    return GCTraces(
        sample_name=inj.sample_name,
        folder_name=inj.folder_name,
        experiment_name=inj.experiment_name,
        time_min=time_min,
        mz=mz,
        eic=eic,
        tic=tic,
        fid_time_min=fid_time,
        fid_signal=fid_signal,
    )
