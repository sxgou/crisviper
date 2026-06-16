"""Dataclass definitions for amplicon structure configuration and YAML config loading."""

from dataclasses import dataclass, field
from typing import List, Dict, Optional


@dataclass
class CutsiteRegion:
    """Definition of a single cutsite region on the reference sequence.

    All coordinates are 0-indexed and inclusive (both start and end positions
    are part of the cutsite). The length in bp is end - start + 1.
    """
    name: str               # Display name, e.g. "Target1", "Target2"
    start: int              # Start position on the reference sequence (0-indexed)
    end: int                # End position on the reference sequence (0-indexed, inclusive)


@dataclass
class AmpliconConfig:
    """Amplicon structure configuration — centralized CARLIN-specific parameters.

    All hard-coded parameters specific to a particular multi-target lineage
    tracer amplicon design are collected here. Modify these values to adapt
    the pipeline to different amplicon architectures (different number of
    targets, linker lengths, primer positions, etc.).

    The reference sequence layout is:
      Primer5 | CGCCG-prefix | [Target + Linker] × N_targets | Postfix | Primer3
    """
    primer5_len: int = 23       # 5' primer length in bp
    primer3_len: int = 33       # 3' primer length in bp
    prefix: str = "CGCCG"       # 5' prefix sequence (immediately after Primer5)
    postfix_len: int = 8        # 3' postfix length in bp (before Primer3)

    target_size: int = 20       # Each target region: conserved 13bp + cutsite 7bp
    linker_size: int = 7        # PAM/Linker length between targets
    n_targets: int = 10         # Number of target sites in the amplicon
    cutsite_offset: int = 13    # Cutsite start offset within a target (after conserved region)
    cutsite_len: int = 7        # Cutsite length in bp

    dual_anchor_tolerance: int = 4  # Allowed mismatches for dual primer anchoring

    @property
    def period(self) -> int:
        """Combined period of one Target + one Linker (bp)."""
        return self.target_size + self.linker_size

    @property
    def expected_full_length(self) -> int:
        """Expected full-length amplicon including primers (bp).

        Layout: Primer5 + prefix + N_targets×(target) + (N_targets-1)×(linker) + postfix + Primer3
        """
        return (self.primer5_len + len(self.prefix) +
                self.n_targets * self.target_size +
                (self.n_targets - 1) * self.linker_size +
                self.postfix_len + self.primer3_len)

    @classmethod
    def carlin_standard(cls) -> "AmpliconConfig":
        """Standard CARLIN amplicon (332 bp) configuration with default values."""
        return cls()

    @classmethod
    def from_json(cls, path: str) -> "AmpliconConfig":
        """Load configuration from a JSON file."""
        import json
        with open(path) as f:
            data = json.load(f)
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})

    @classmethod
    def from_dict(cls, data: dict) -> "AmpliconConfig":
        """Load from a dict (e.g. parsed from YAML 'amplicon' section)."""
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


def load_yaml_config(path: str) -> dict:
    """Load a YAML configuration file.

    Args:
        path: Path to the YAML file.

    Returns:
        Dict with the following possible keys:
        - amplicon: AmpliconConfig field dict
        - cutsites: list of cutsite region dicts [{"name": ..., "start": ..., "end": ...}, ...]
        - pipeline: PipelineConfig field dict (defined in models.py)
    """
    try:
        import yaml
    except ImportError:
        raise ImportError("PyYAML is required. Install with: pip install pyyaml")
    try:
        with open(path) as f:
            data = yaml.safe_load(f)
        return data if isinstance(data, dict) else {}
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file not found: {path}")
    except yaml.YAMLError as e:
        raise ValueError(f"YAML format error: {e}")


def cutsites_from_list(data: list) -> List["CutsiteRegion"]:
    """Create a list of CutsiteRegion from a list of dicts (e.g. YAML 'cutsites' section)."""
    result = []
    for i, entry in enumerate(data):
        if not isinstance(entry, dict):
            continue
        result.append(CutsiteRegion(
            name=entry.get("name", f"Target{i+1}"),
            start=entry["start"],
            end=entry["end"],
        ))
    return result


