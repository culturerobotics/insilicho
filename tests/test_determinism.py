import copy
import typing
from datetime import datetime

import numpy as np
import pandas as pd
import pytest
from scipy.stats import qmc

from tests import run_doe


def _create_lhc_settings(
    n_runs: int = 24,
    random_seed: int = 42,
) -> typing.List[typing.Any]:
    sampler = qmc.LatinHypercube(d=len(run_doe.ranges.keys()), seed=random_seed)
    l_bounds = [r[0] for r in run_doe.ranges.values()]
    u_bounds = [r[1] for r in run_doe.ranges.values()]
    exps = qmc.scale(sampler.random(n=n_runs), l_bounds, u_bounds)
    return list(exps)


def _run_initial_set(n_runs):
    observations = []

    for e in _create_lhc_settings(n_runs, int(datetime.now().timestamp())):
        condition = {k: v for k, v in zip(run_doe.ranges.keys(), e)}
        results, score = run_doe.run_exp(condition, sampling_stddev=0)
        observations.append(
            {
                "condition": condition,
                "results": results,
                "full_result": copy.deepcopy(run_doe.RealCHO.full_result),
                "score": score,
            }
        )
    return pd.DataFrame(observations)


def eval_param_mse(parameters, condition_idxs, ref_observations):
    def calc_sse(ref_results, results):
        assert set(ref_results.keys()) == set(results.keys())
        ret = 0
        for k in ref_results.keys():
            ret += np.sum((np.array(ref_results[k]) - np.array(results[k])) ** 2)
        return ret

    my_cho = copy.deepcopy(run_doe.RealCHO)
    for k, v in parameters.items():
        setattr(my_cho, k, v)

    sse = 0
    for i in condition_idxs:
        r, _ = run_doe.run_exp(
            ref_observations["condition"][i], model=my_cho, sampling_stddev=0
        )
        sse += calc_sse(ref_observations["results"][i], r)
    return sse / len(condition_idxs)


def get_mismatches(n_runs, ref_observations):
    return [
        i for i in range(n_runs) if (eval_param_mse(dict(), [i], ref_observations) > 0)
    ]


class TestDeterminism:
    @pytest.mark.timeout(60)
    def test_determinism(self) -> None:
        n_runs = 128
        obs = _run_initial_set(n_runs)
        assert obs.shape[0] == n_runs
        assert not get_mismatches(n_runs, obs)
