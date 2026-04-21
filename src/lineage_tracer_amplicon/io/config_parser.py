"""
Configuration parser for amplicon structure.
"""

import json
import yaml
from pathlib import Path
from typing import Union, Dict, Any
from ..core.structures import AmpliconStructure


def load_config(config_path: Union[str, Path]) -> AmpliconStructure:
    """
    Load amplicon structure configuration from JSON or YAML file.
    
    Args:
        config_path: Path to configuration file
        
    Returns:
        AmpliconStructure object
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        if config_path.suffix.lower() in ['.yaml', '.yml']:
            config_data = yaml.safe_load(f)
        else:
            # Assume JSON
            config_data = json.load(f)
    
    return parse_config(config_data)


def parse_config(config_data: Dict[str, Any]) -> AmpliconStructure:
    """
    Parse configuration dictionary to AmpliconStructure.
    
    Args:
        config_data: Configuration dictionary
        
    Returns:
        AmpliconStructure object
    """
    # Validate required fields
    required_fields = ["reference", "features"]
    for field in required_fields:
        if field not in config_data:
            raise ValueError(f"Missing required field: {field}")
    
    # Ensure features is a list
    features = config_data["features"]
    if not isinstance(features, list):
        raise TypeError("features must be a list")
    
    # Validate each feature
    for i, feature in enumerate(features):
        if not isinstance(feature, dict):
            raise TypeError(f"feature {i} must be a dictionary")
        
        required_feature_fields = ["name", "start", "end", "type"]
        for field in required_feature_fields:
            if field not in feature:
                raise ValueError(f"feature {i} missing field: {field}")
    
    return AmpliconStructure(
        reference=config_data["reference"],
        features=features
    )