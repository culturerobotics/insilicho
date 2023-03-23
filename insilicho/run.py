import dataclasses
import types
import typing

import numpy as np
import yaml

from insilicho import growth_model, parameters, plotter, solver, util


def add_relative_normal_noise(
    a: typing.Union[np.ndarray, float, int], relative_std_dev: float
) -> np.ndarray:
    """_summary_

    Args:
        a (typing.Union[np.ndarray, float, int]): Original value(s) that will act as
            means.
        relative_std_dev (float): Relative standard deviation. The resulting noise
            sample will be taken from a normal distribution of mean = 1, standard
            deviation = `relative_std_dev`. The noise sample(s) will then be multiplied
            to value(s) in `a`.

    Returns:
        np.ndarray: an array of values with added noise.
    """
    if isinstance(a, float) or isinstance(a, int):
        array = np.array([a])
    else:
        array = a

    return array * np.random.normal(
        loc=1.0, scale=float(relative_std_dev), size=len(array)
    )


class GrowCHO:
    def __init__(
        self,
        config: typing.Union[typing.Dict[str, typing.Any], str],
        feed_fn: typing.Optional[growth_model.FeedFunctionType],
        temp_fn: typing.Optional[growth_model.TempFunctionType],
        random_seed: int = 0,
        param_rel_stddev: float = 0.05,
        solver_max_step_size: float = np.inf,
    ):
        """Class to simulate CHO growth.

        Args:
            config (typing.Union[typing.Dict[str, typing.Any], str]): A path to yaml
                file or a dictionary with initial conditions and parameter values.
            feed_fn (typing.Optional[growth_model.FeedFunctionType]): A callable
                describing time dependence of feed profile, expected units for feed rate
                are in L/h.
            temp_fn (typing.Optional[growth_model.TempFunctionType]): A callable
                describing time dependence of temperature profile, expected units for
                temp are in degC.
            random_seed (int, optional): random seed to control sampling events.
                Defaults to 0.
            param_rel_stddev (float, optional): Relative std deviation while sampling
                parameter, assumes a normal distribution.. Defaults to 0.05.
            solver_max_step_size (float, optional): Max step size for the odeint solver
                to take. Defaults to np.inf.
        """
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

    def _randomize_params(self, rel_stddev: float):
        """Randomize parameters for the model.

        Args:
            rel_stddev (float): Relative std deviation while sampling parameter, assumes
                a normal distribution.
        """
        noisy_params = self.params_with_noise()
        for pname, pval in noisy_params.items():
            new_val = add_relative_normal_noise(pval, rel_stddev)[0]
            setattr(self.params, pname, new_val)

    def params_with_noise(self):
        """Currently we add noise to a subset of params"""
        return {
            f.name: getattr(self.params, f.name)
            for f in dataclasses.fields(self.params)
            if (
                f.type == typing.Union[float, str]
                and f.name
                not in [
                    "Cglc_feed",
                    "Cgln_feed",
                    "T_optimal",
                    "T_optimal_decay_spread",
                    "pH_optimal",
                    "pH_optimal_decay_spread",
                ]
            )
        }

    def execute(
        self,
        initial_conditions: typing.Optional[typing.Dict[str, typing.Any]] = None,
        plot: bool = False,
        sampling_stddev: float = 0.05,
        starting_at_day: int = 0,
    ) -> typing.Dict[str, typing.Any]:
        """Execute the GrowCHO model object

        Args:
            initial_conditions (typing.Optional[typing.Dict[str, typing.Any]],
                optional): Initial conditions for the solver. Defaults to conditions
                given in insilicho.parameters.
            plot (bool, optional): option to plot data using matplotlib. Defaults to
                False.
            sampling_stddev (float, optional): scale of error in normal distributed
                sampling event, relative to sample magnitude. Defaults to 0.05.
            starting_at_day (int, optional): day at which to start the simulation.
                Defaults to 0.

        Raises:
            IOError: If initial conditions were not supplied.
            RuntimeError: If integration/LSODA solver runs into failures.

        Returns:
            typing.Dict[str, typing.Any]: Sampled metabolite, volume and cell
                concentrations.

        """

        if initial_conditions:
            self.initial_conditions = util.DataClassUnpack.instantiate(
                parameters.InitialConditions, initial_conditions
            )

        if not self.initial_conditions:
            raise IOError("Initial conditions undefined for sim")

        tmin = starting_at_day * 24
        tspan = np.linspace(
            tmin,
            tmin + 24 * self.params.Ndays,
            1000 * self.params.Ndays,
        )

        state, state_vars, infodict = solver.solve(
            self.params,
            self.initial_conditions,
            tspan=tspan,
            feed_fn=self.feed_fn,
            temp_fn=self.temp_fn,
            solver_hmax=self.solver_max_step_size,
        )
        if plot:
            plotter.plot(tspan, state, state_vars)

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
            sampling_rel_stddev=sampling_stddev,
        )

    @property
    def full_result(self):
        """A property to get the full unsampled data."""
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
    state: np.ndarray,
    state_vars: np.ndarray,
    params: parameters.InputParameters,
    tspan: np.ndarray,
    sampling_rel_stddev: float = 0.05,
) -> typing.Dict[str, typing.Any]:
    """Samples datapoints from a simulation output.

    Args:
        state (np.ndarray): Array of state solutions (Xv, Xt, Cglc, Cgln, Clac, Camm,
            Cmab, Coxygen, V, pH) for all points in tspan.
        state_vars (np.ndarray): Array of state variable (F, T, mu, mu_d, q_glc, q_gln,
            q_lac, q_amm, q_mab, Osmolarity) solutions for all points in tspan.
        params (parameters.InputParameters): Input parameters for simulation system.
        tspan (np.ndarray): time array (in hours) over which the system was solved.
        sampling_rel_stddev (float, optional): scale of error in normal distributed
            sampling event, relative to sample magnitude. Defaults to 0.05.

    Returns:
        typing.Dict[str, typing.Any]: Results from sampling i.e., Xv, Xt, Cglc, Cgln,
            Clac, Camm, Cmab, Osmolarity and time.
    """

    Xv, Xt, Cglc, Cgln, Clac, Camm, Cmab, Coxygen, V, pH = state.transpose()
    Osmolarity = state_vars[:, 9]
    time = tspan.transpose()
    res_map = {
        "time": time,  # hrs
        "Xv": Xv * 1e-9,  # viable cells (millions/mL conversion)
        "Xt": Xt * 1e-9,  # total cells (millions/mL conversion)
        "Cglc": Cglc,  # mmol
        "Cgln": Cgln,
        "Clac": Clac,
        "Camm": Camm,
        "Cmab": Cmab,
        "Coxygen": Coxygen,
        "V": V,
        "pH": pH,
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
        if k in ["time", "V"]:
            res[k] = var[idx].tolist()
        else:
            res[k] = np.maximum(
                add_relative_normal_noise(var[idx], sampling_rel_stddev),
                parameters.EPSILON,
            ).tolist()

    return res
