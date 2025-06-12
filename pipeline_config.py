"""Load default pipeline config."""

import pathlib

import yaml


default_config_file = pathlib.Path(__file__).parent / "config" / "default_config.yml"

with open(default_config_file, "r", encoding="utf-8") as default_config:
    WORKFLOW_DEFAULT_CONF = yaml.safe_load(default_config)

