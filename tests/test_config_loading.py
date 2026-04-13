import os
import tempfile

import pytest
import yaml

from insilicho import parameters, run


class TestYAMLConfigLoading:
    def test_load_from_yaml_file(self):
        config = {
            "parameters": {"K_lys": 0.05, "Ndays": 6},
            "initial_conditions": {"V": 0.025},
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(config, f)
            path = f.name

        try:
            model = run.GrowCHO(
                path,
                feed_fn=lambda t: 0.003,
                temp_fn=lambda t: 36.4,
            )
            assert model.params.Ndays == 6
            assert model.initial_conditions.V == 0.025
        finally:
            os.unlink(path)

    def test_load_with_unit_strings_from_yaml(self):
        config = {
            "parameters": {"K_lys": "0.05 1/h"},
            "initial_conditions": {"V": "50 mL"},
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(config, f)
            path = f.name

        try:
            model = run.GrowCHO(
                path,
                feed_fn=lambda t: 0.003,
                temp_fn=lambda t: 36.4,
            )
            assert model.initial_conditions.V == pytest.approx(0.05)
        finally:
            os.unlink(path)

    def test_invalid_yaml_raises(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            # Write content that yaml.safe_load can parse but gives
            # bad structure — config_parser itself won't raise on valid
            # YAML, but an empty file gives None which triggers downstream
            # errors
            f.write("")
            path = f.name

        try:
            with pytest.raises((TypeError, AttributeError)):
                run.GrowCHO(
                    path,
                    feed_fn=lambda t: 0.003,
                    temp_fn=lambda t: 36.4,
                )
        finally:
            os.unlink(path)

    def test_default_params_when_no_config(self):
        params, ic = run.unpack()
        assert isinstance(params, parameters.InputParameters)
        assert isinstance(ic, parameters.InitialConditions)
        # Should have default values
        assert params.Ndays == 12
        assert ic.V == pytest.approx(40 / 1000)


class TestConfigParser:
    def test_parses_valid_yaml(self):
        config = {"parameters": {"Ndays": 8}}
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(config, f)
            path = f.name

        try:
            data = run.config_parser(path)
            assert data["parameters"]["Ndays"] == 8
        finally:
            os.unlink(path)
