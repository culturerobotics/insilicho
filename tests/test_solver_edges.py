import pytest

from insilicho import parameters, run


class TestSolverFailure:
    def test_integration_failure_raises_runtime_error(self):
        """Extreme parameters that cause the ODE solver to fail."""
        bad_config = {
            "parameters": {
                "mu_max": 1e10,
                "mu_d_max": 1e10,
                "q_glc_max": 1e10,
                "q_gln_max": 1e10,
                "Ndays": 1,
            },
            "initial_conditions": {
                "Xv": 1e20,
                "Xt": 1e20,
                "V": 0.001,
            },
        }
        model = run.GrowCHO(
            bad_config,
            feed_fn=lambda t: 0.003,
            temp_fn=lambda t: 36.4,
            param_rel_stddev=0.0,
        )
        with pytest.raises(RuntimeError, match="Integration failed"):
            model.execute()


class TestMissingInitialConditions:
    def test_raises_without_initial_conditions(self):
        """When unpack returns None for initial_conditions."""
        model = run.GrowCHO(
            {"parameters": {"Ndays": 1}},
            feed_fn=lambda t: 0.003,
            temp_fn=lambda t: 36.4,
        )
        # Force initial_conditions to None to trigger the IOError path
        model.initial_conditions = None
        with pytest.raises(IOError, match="Initial conditions undefined"):
            model.execute()


class TestSolverDefaults:
    def test_solve_with_default_tspan(self):
        """solver.solve with tspan=None uses default linspace."""
        from insilicho import solver

        params = parameters.InputParameters(Ndays=1)
        ic = parameters.InitialConditions()
        state, state_vars, info = solver.solve(
            params,
            ic,
            feed_fn=lambda t: 0.003,
            temp_fn=lambda t: 36.4,
        )
        assert info["message"] == "Integration successful."
        assert state.shape[0] == 10000  # default tspan length

    def test_execute_with_starting_day_offset(self):
        model = run.GrowCHO(
            {"parameters": {"Ndays": 2}, "initial_conditions": {"V": 0.025}},
            feed_fn=lambda t: 0.003,
            temp_fn=lambda t: 36.4,
        )
        result = model.execute(starting_at_day=3)
        # Time should start at day 3 = 72 hours
        assert result["time"][0] == pytest.approx(72.0)
