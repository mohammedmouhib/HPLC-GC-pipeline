"""
Agilent Chemstation .D folder discovery and chromatogram loading.

MOCCA2 reads .D folders natively, so this module only handles:
- locating .D injection folders (flat or one level of run subfolders)
- parsing the sample name out of each sample.xml
- loading each chromatogram at the configured wavelength

Kept deliberately close to the original ``find_injection_folders`` /
``prepare_dataset`` logic; the change is that the experiment directory and
wavelength are passed in rather than hardcoded.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from mocca2 import Chromatogram


@dataclass
class Injection:
    """One .D injection folder and its parsed metadata."""
    folder_path: Path
    folder_name: str        # e.g. "007-0301.D"
    xml_file: Path
    sample_name: str
    experiment_name: str    # run subfolder name, or "Main" for flat layout


def parse_sample_name_from_xml(xml_path: Path) -> str:
    """Extract the sample name from an Agilent sample.xml.

    Falls back to the .D folder name if the <Name> tag can't be read.
    """
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()

        name_element = root.find("Name")
        if name_element is not None and name_element.text:
            return name_element.text.strip()

        for tag in ["SampleName", "Sample_Name", "sample_name"]:
            element = root.find(f".//{tag}")
            if element is not None and element.text:
                return element.text.strip()

        print(f"  Warning: no <Name> tag in {xml_path.name}; using folder name")
        return xml_path.parent.name
    except Exception as exc:  # malformed XML, encoding issues, etc.
        print(f"  Warning: could not parse {xml_path.name}: {exc}")
        return xml_path.parent.name


def find_injection_folders(experiment_dir: Path) -> list[Injection]:
    """Find all .D folders containing a sample.xml.

    Supports a flat layout (``.D`` folders directly under ``experiment_dir``)
    and one level of run subfolders (``experiment_dir/<run>/*.D``).
    """
    experiment_dir = Path(experiment_dir)
    injections: list[Injection] = []

    def scan(directory: Path, parent_name: str = "") -> None:
        for item in sorted(directory.iterdir()):
            if not item.is_dir():
                continue

            if item.name.endswith(".D"):
                xml_file = _find_sample_xml(item)
                if xml_file is None:
                    print(f"  Warning: {item.name} has no sample.xml; skipping")
                    continue
                injections.append(
                    Injection(
                        folder_path=item,
                        folder_name=item.name,
                        xml_file=xml_file,
                        sample_name=parse_sample_name_from_xml(xml_file),
                        experiment_name=parent_name or "Main",
                    )
                )
            elif parent_name == "":
                # Only recurse one level, from the top, to avoid scanning the
                # whole disk. Method/result folders end in .M/.B and are ignored.
                scan(item, item.name)

    scan(experiment_dir)
    injections.sort(key=lambda inj: (inj.experiment_name, inj.folder_name))
    return injections


def _find_sample_xml(d_folder: Path) -> Optional[Path]:
    for xml_name in ["sample.xml", "Sample.xml", "SAMPLE.XML"]:
        candidate = d_folder / xml_name
        if candidate.exists():
            return candidate
    return None


def load_chromatogram(d_folder: Path, wavelength: int) -> Chromatogram:
    """Load a .D folder via MOCCA2 and extract a single wavelength."""
    chrom = Chromatogram(str(d_folder))
    chrom.extract_wavelength(wavelength, wavelength, inplace=True)
    return chrom
