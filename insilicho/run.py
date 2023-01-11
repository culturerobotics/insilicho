import dataclasses
import types
import typing

import numpy as np
import yaml

from insilicho import parameters, plotter, solver, util


class GrowCHO:
    def __init__(
        self,
        config: typing.Union[typing.Dict[str, typing.Any], str],
        feed_fn: typing.Optional[typing.Callable[[typing.Any], typing.Any]],
        temp_fn: typing.Optional[typing.Callable[[typing.Any], typing.Any]],
        random_seed: int = 0,
        param_rel_stddev: float = 0.05,
        solver_max_step_size=np.inf,
    ):
        cfg_dict, cfg_path = None, None
        if type(config) == dict:
            cfg_dict = config
        elif type(config) == str:
            cfg_path = config

        self.params, self.initial_conditions = unpack(cfg_dict, cfg_path)
        self.seed = np.random.seed(random_seed)
        self._randomize_params(param_rel_stddev)

        self.feed_fn = feed_fn
        self.temp_fn = temp_fn
        self.solver_max_step_size = solver_max_step_size

        self._full_result = types.SimpleNamespace(
            state=[], state_vars=[], t=[], info={}
        )

    def _randomize_params(self, rel_stddev):
        float_params = {
            f.name: getattr(self.params, f.name)
            for f in dataclasses.fields(self.params)
            if (
                f.type == typing.Union[float, str]
                and f.name not in ["Cglc_feed", "Cgln_feed"]
            )
        }
        for pname, pval in float_params.items():
            new_val = np.random.normal(loc=pval, scale=rel_stddev * pval)
            setattr(self.params, pname, new_val)

    def execute(
        self,
        initial_conditions: typing.Optional[typing.Dict[str, typing.Any]] = None,
        plot: bool = False,
        sampling_stddev: float = 0.05,
    ):
        if initial_conditions:
            self.initial_conditions = util.DataClassUnpack.instantiate(
                parameters.InitialConditions, initial_conditions
            )

        if not self.initial_conditions:
            raise IOError("Initial conditions undefined for sim")

        tspan = np.linspace(0, 24 * self.params.Ndays, 1000 * self.params.Ndays)
        solver_fn = solver.solve
        if plot:
            solver_fn = plotter.solve_and_plot

        state, state_vars, infodict = solver_fn(
            self.params,
            self.initial_conditions,
            tspan=tspan,
            feed_fn=self.feed_fn,
            temp_fn=self.temp_fn,
            solver_hmax=self.solver_max_step_size,
        )
        self._full_result = types.SimpleNamespace(
            state=state,
            state_vars=state_vars,
            t=tspan,
            info=infodict,
        )

        if infodict["message"] != "Integration successful.":
            raise RuntimeError(
                "Integration failed at specified params and/or initial values."
            )

        return flex2_sampling(
            state,
            state_vars,
            self.params,
            tspan,
            sampling_stddev=sampling_stddev,
        )

    @property
    def full_result(self):
        return self._full_result


def config_parser(cfg_path):
    data = {}
    with open(cfg_path, "r") as f:
        try:
            data = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            raise ValueError(f"Error parsing sim config: {exc}")
    return data


def unpack(cfg_dict=None, cfg_path=None):
    def instantiate_objects_from_dict(data):
        return util.DataClassUnpack.instantiate(
            parameters.InputParameters, data.get("parameters")
        ), util.DataClassUnpack.instantiate(
            parameters.InitialConditions, data.get("initial_conditions")
        )

    if cfg_dict:
        params, ic = instantiate_objects_from_dict(cfg_dict)
    elif cfg_path:
        params, ic = instantiate_objects_from_dict(config_parser(cfg_path))
    else:
        params = parameters.InputParameters()
        ic = parameters.InitialConditions()
    return params, ic


def flex2_sampling(
    state,
    state_vars,
    params,
    tspan,
    sampling_stddev=0.05,
):
    Xv, Xt, Cglc, Cgln, Clac, Camm, Cmab, Coxygen, V, pH = state.transpose()
    Osmolarity = state_vars[:, 9]
    time = tspan.transpose()
    res_map = {
        "time": time,  # hrs
        "Xv": Xv * 1e-9,  # viable cells (millions/mL conversion)
        "Cglc": Cglc,  # mmol
        "Cgln": Cgln,
        "Clac": Clac,
        "Camm": Camm,
        "Cmab": Cmab,
        "Osmolarity": Osmolarity,
    }
    res = {}

    # get idx to sample at
    idx = np.round(
        np.linspace(0, len(Xv) - 1, params.Ndays * params.Nsamples + 1)
    ).astype(int)

    # sample across small time range, add noise
    for k, var in res_map.items():
        var = np.maximum(var, parameters.EPSILON)
        if k == "time":
            res[k] = var[idx].tolist()
        else:
            res[k] = np.maximum(
                np.random.normal(
                    loc=var[idx], scale=float(sampling_stddev), size=len(idx)
                ),
                parameters.EPSILON,
            ).tolist()

    return res
