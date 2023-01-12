import pytest

from insilicho import growth_model, parameters


class TestGrowthModel:
    def test_move_away_from_optima_reduces_values(self):
        T_optima = 36.4
        spread = 3.12

        T = 36.4
        assert (
            growth_model.exponential_dependence_around_optima(T, T_optima, spread)
            == 1.0
        )

        T = 36.0
        assert (
            growth_model.exponential_dependence_around_optima(T, T_optima, spread) < 1.0
        )

        T = 37.0
        assert (
            growth_model.exponential_dependence_around_optima(T, T_optima, spread) < 1.0
        )

    def test_missing_fns_raise_errors(self):
        with pytest.raises(IOError):
            assert growth_model.state_vars(1, ([1] * 10), parameters.InputParameters)
