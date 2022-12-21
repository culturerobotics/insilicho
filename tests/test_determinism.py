import copy
import typing
from datetime import datetime

import numpy as np
import pandas as pd
import pytest
from scipy.stats import qmc

from insilicho import run
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
    @pytest.mark.timeout(90)
    def test_determinism(self) -> None:
        n_runs = 64
        obs = _run_initial_set(n_runs)
        assert obs.shape[0] == n_runs
        assert not get_mismatches(n_runs, obs)


class TestBlips:
    N_ITERS = 3
    KNOWN_BLIP_CONDITIONS = [
        {
            "batch_glc": 50.0,
            "batch_gln": 50.0,
            "batch_pH": 7.085020796463857,
            "feed_glc": 100.0,
            "feed_gln": 100.0,
            "prod_start_eft": 40.57103186367486,
            "batch_temp": 36.36799201861109,
            "prod_temp": 35.188648493781166,
            "day_0_feed": 0,
            "day_1_feed": 0.0021022266574628997,
            "day_2_feed": 0.00012093627022496579,
            "day_3_feed": 0.0033983883685194784,
        },
        {
            "batch_glc": 50.0,
            "batch_gln": 50.0,
            "batch_pH": 6.750062183704625,
            "feed_glc": 100.0,
            "feed_gln": 100.0,
            "prod_start_eft": 30.738714263859187,
            "batch_temp": 36.130589554264475,
            "prod_temp": 36.17424156024803,
            "day_0_feed": 0,
            "day_1_feed": 0.003580592381995898,
            "day_2_feed": 0.000524712751487901,
            "day_3_feed": 0.0017155253782759289,
        },
        {
            "batch_glc": 50.0,
            "batch_gln": 100.0,
            "batch_pH": 7.050000000000001,
            "feed_glc": 500.0,
            "feed_gln": 300.0,
            "prod_start_eft": 72.0,
            "batch_temp": 36.0,
            "prod_temp": 36.0,
            "day_0_feed": 0.0,
            "day_1_feed": 0.0,
            "day_2_feed": 0.005,
            "day_3_feed": 0.005,
        },
    ]

    @pytest.mark.timeout(15)
    def test_blips(self):
        for _ in range(self.N_ITERS):
            for conditions in self.KNOWN_BLIP_CONDITIONS:
                model = run.GrowCHO(
                    {"parameters": {"Ndays": 4, "Nsamples": 2}},
                    None,
                    None,
                )
                run_doe.run_exp(
                    conditions,
                    model=model,
                    plot=False,
                )
                vol = model.full_result.state[:, 8]
                ph = model.full_result.state[:, 9]
                assert np.argwhere(ph <= 0).shape[0] == 0
                assert np.argwhere(vol <= 0).shape[0] == 0
