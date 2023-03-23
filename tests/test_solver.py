import pytest

from insilicho import run


class TestBolusFeed:
    def test_solver_with_bolus_feed(self, bolus_feed: run.GrowCHO):
        bolus_feed.execute(plot=False)

        assert bolus_feed.full_result.info["message"] == "Integration successful."
        # Ending concentrations
        assert (bolus_feed.full_result.state[-1, 0]) * 1e-9 == pytest.approx(
            26.1, rel=0.005
        )  # Xv
        assert (bolus_feed.full_result.state[-1, 1]) * 1e-9 == pytest.approx(
            26.7, rel=0.005
        )  # Xt
        assert (bolus_feed.full_result.state[-1, 2]) == pytest.approx(
            35.7, rel=0.005
        )  # CGlc


class TestConstantFeed:
    def test_solver_with_constant_feed(self, constant_feed: run.GrowCHO):
        constant_feed.execute(plot=False)

        assert constant_feed.full_result.info["message"] == "Integration successful."
        # Ending concentrations
        assert (constant_feed.full_result.state[-1, 0]) * 1e-9 == pytest.approx(
            60.1, rel=0.005
        )  # Xv
        assert (constant_feed.full_result.state[-1, 1]) * 1e-9 == pytest.approx(
            63.9, rel=0.005
        )  # Xt
        assert (constant_feed.full_result.state[-1, 2]) == pytest.approx(
            0.069, rel=0.005
        )  # CGlc

    def test_example(self):
        def T(time):
            """returns temperature in degC"""
            return 36.4

        def F(time):
            """returns flow rate in L/hr"""
            return 0.003

        model = run.GrowCHO(
            {"parameters": {"K_lys": "0.05 1/h"}},
            feed_fn=F,
            temp_fn=T,
        )

        model.execute(plot=False, initial_conditions={"V": "50 mL"})
        assert model.initial_conditions.V == 50 / 1000
        final_V = model.full_result.state[-1, 8]
        assert final_V == pytest.approx(0.914)
