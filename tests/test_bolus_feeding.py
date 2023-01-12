import pytest

from insilicho import run


def F(time, bolus_size=0.03, num_bolus=10, bolus_frequency=24):
    def dirac(t):
        # https://stackoverflow.com/questions/63230568/numerical-solution-to-a-differential-equation-containing-a-dirac-delta-function
        def tentfunc(t):
            return max(0, 1 - abs(t))

        N = 10.0
        return sum(
            N * tentfunc(N * (t - (i + 1) * bolus_frequency)) for i in range(num_bolus)
        )

    # if 0 < V <= reactor.Vessel.volume:
    return bolus_size * dirac(time)  # L/h
    # return 0


def T(time):
    return 37


@pytest.fixture
def model_obj():
    cfg_dict = {"parameters": {"K_lys": "0.05 1/h"}, "initial_conditions": {"V": 0.025}}
    return run.GrowCHO(cfg_dict, feed_fn=F, temp_fn=T, solver_max_step_size=0.1)


class TestBolusFeed:
    def test_solver_with_bolus(self, model_obj: run.GrowCHO):
        model_obj.execute(plot=False)

        assert model_obj.full_result.info["message"] == "Integration successful."
        # Ending concentrations
        assert (model_obj.full_result.state[-1, 0]) == pytest.approx(
            26446818993.382
        )  # Xv
        assert (model_obj.full_result.state[-1, 1]) == pytest.approx(
            27031305161.744
        )  # Xt
        assert (model_obj.full_result.state[-1, 2]) == pytest.approx(34.266)  # CGlc
