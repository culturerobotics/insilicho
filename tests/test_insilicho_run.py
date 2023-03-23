import copy
import dataclasses

import numpy as np
import pytest

from insilicho import run


class TestInsilichoRun:
    def test_normal_noise(self):
        assert run.add_relative_normal_noise(2, 0.1)

        assert isinstance(
            run.add_relative_normal_noise(np.array([1, 1.5, 2]), 0.05),
            np.ndarray,
        )

        with pytest.raises(ValueError):
            run.add_relative_normal_noise(np.array([1, 1.5, 2]), -0.05)

    def test_flex2_sampling(self, constant_feed: run.GrowCHO):
        constant_feed.execute(plot=False)
        result_without_noise = run.flex2_sampling(
            constant_feed.full_result.state,
            constant_feed.full_result.state_vars,
            constant_feed.params,
            np.linspace(
                0, 24 * constant_feed.params.Ndays, 1000 * constant_feed.params.Ndays
            ),
            0.0,
        )
        result_with_noise = run.flex2_sampling(
            constant_feed.full_result.state,
            constant_feed.full_result.state_vars,
            constant_feed.params,
            np.linspace(
                0, 24 * constant_feed.params.Ndays, 1000 * constant_feed.params.Ndays
            ),
            0.07,
        )

        allow_list = ["time", "V"]

        for (k, v1), v2 in zip(
            result_without_noise.items(), result_with_noise.values()
        ):
            if k in allow_list:
                assert v1 == v2
            else:
                assert v1 != v2


class TestGrowCHO:
    def test_random_noise(self, constant_feed: run.GrowCHO):
        grow_cho_copy = copy.deepcopy(constant_feed)
        original_params = copy.copy(grow_cho_copy.params)

        grow_cho_copy._randomize_params(0.05)
        new_params = copy.copy(grow_cho_copy.params)

        noisy_params = constant_feed.params_with_noise()

        for f, v1, v2 in zip(
            dataclasses.fields(grow_cho_copy.params),
            original_params.tolist(),
            new_params.tolist(),
        ):
            if f.name in noisy_params:
                assert v1 != v2
