import pytest

from insilicho import run


CFG_DICT = {"parameters": {"K_lys": "0.05 1/h"}, "initial_conditions": {"V": 0.025}}


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
