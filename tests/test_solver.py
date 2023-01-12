import pytest

from insilicho import run

CFG_DICT = {"parameters": {"K_lys": "0.05 1/h"}, "initial_conditions": {"V": 0.025}}


@pytest.fixture
def bolus_feed():
    def F(time, bolus_size=0.03, num_bolus=10, bolus_frequency=24):
        def dirac(t):
            # https://stackoverflow.com/questions/63230568/numerical-solution-to-a-differential-equation-containing-a-dirac-delta-function
            def tentfunc(t):
                return max(0, 1 - abs(t))

            N = 10.0
            return sum(
                N * tentfunc(N * (t - (i + 1) * bolus_frequency))
                for i in range(num_bolus)
            )

        return bolus_size * dirac(time)  # L/h

    def T(time):
        return 36.4

    # bolus feed requires manual reduction in step size to capture all steps.
    return run.GrowCHO(
        CFG_DICT,
        feed_fn=F,
        temp_fn=T,
        solver_max_step_size=0.1,
    )


@pytest.fixture
def constant_feed():
    def F(time):
        return 0.003

    def T(time):
        return 36.4

    return run.GrowCHO(
        CFG_DICT,
        feed_fn=F,
        temp_fn=T,
    )


class TestBolusFeed:
    def test_solver_with_bolus_feed(self, bolus_feed: run.GrowCHO):
        bolus_feed.execute(plot=False)

        assert bolus_feed.full_result.info["message"] == "Integration successful."
        # Ending concentrations
        assert (bolus_feed.full_result.state[-1, 0]) * 1e-9 == pytest.approx(
            26.5, rel=0.005
        )  # Xv
        assert (bolus_feed.full_result.state[-1, 1]) * 1e-9 == pytest.approx(
            27.0, rel=0.005
        )  # Xt
        assert (bolus_feed.full_result.state[-1, 2]) == pytest.approx(
            34.266, rel=0.005
        )  # CGlc


class TestConstantFeed:
    def test_solver_with_constant_feed(self, constant_feed: run.GrowCHO):
        constant_feed.execute(plot=False)

        assert constant_feed.full_result.info["message"] == "Integration successful."
        # Ending concentrations
        assert (constant_feed.full_result.state[-1, 0]) * 1e-9 == pytest.approx(
            60.9, rel=0.005
        )  # Xv
        assert (constant_feed.full_result.state[-1, 1]) * 1e-9 == pytest.approx(
            64.9, rel=0.005
        )  # Xt
        assert (constant_feed.full_result.state[-1, 2]) == pytest.approx(
            0.068, rel=0.005
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
