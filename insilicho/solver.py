import typing

import numpy as np
from scipy.integrate import odeint

from insilicho import growth_model, parameters


def solve(
    params: parameters.InputParameters,
    initial_conditions: parameters.InitialConditions,
    model: typing.Any = growth_model.model,
    tspan: typing.Any = None,
    feed_fn: typing.Optional[growth_model.FeedFunctionType] = None,
    temp_fn: typing.Optional[growth_model.TempFunctionType] = None,
    solver_hmax: float = np.inf,
) -> typing.Tuple[np.ndarray, np.ndarray, typing.Any]:
    """Solves the supplied differential equation system using scipy.odeint (LSODA) solver.

    Args:
        params (parameters.InputParameters): Parameters for the model.
        initial_conditions (parameters.InitialConditions): Initial conditions for the
            solver.
        model (growth_model.model, optional): Differential equations to solve. Defaults
            to growth_model.model.
        tspan (List, optional): time array (in hrs) over which to solve the system.
            Defaults to np.linspace(0, 288, 10000).
        feed_fn (growth_model.FeedFunctionType, optional): Callable describing feed
            profile. Defaults to None.
        temp_fn (growth_model.TempFunctionType, optional): Callable describing temp
            profile. Defaults to None.
        solver_hmax (float, optional): max step size solver can take. Defaults to
            np.inf.

    Returns:
        state_model: Array of state solutions for all points in tspan.
        state_model: Array of state solutions for all points in tspan.
        infodict: Dictionary of LSODA solver behavior.
    """
    if tspan is None:
        tspan = np.linspace(0, 288, 10000)

    # Use default InputParameters and InitialConditions if not provided
    if params is None:
        params = parameters.InputParameters()
    if initial_conditions is None:
        initial_conditions = parameters.InitialConditions()

    IC = initial_conditions.tolist()
    args = params.tolist()

    state_model, info = odeint(
        model,
        IC,
        tspan,
        (
            args,
            feed_fn,
            temp_fn,
        ),
        tfirst=True,
        printmessg=False,
        full_output=True,
        hmax=solver_hmax,
    )
    state_vars = []
    for i in range(len(tspan)):
        state_vars.append(
            list(
                growth_model.state_vars(
                    tspan[i], state_model[i], params, feed_fn, temp_fn
                )
            )
        )
    return state_model, np.array(state_vars, dtype=float), info
