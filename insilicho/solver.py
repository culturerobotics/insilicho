import typing
import numpy as np
from scipy.integrate import odeint

from insilicho import growth_model, parameters


def solve(
    params: typing.Optional[parameters.InputParameters],
    initial_conditions: typing.Optional[parameters.InitialConditions],
    model=growth_model.model,
    tspan=np.linspace(0, 288, 10000),
    feed_fn: typing.Optional[growth_model.FeedFunctionType] = None,
    temp_fn: typing.Optional[growth_model.TempFunctionType] = None,
    solver_hmax=np.inf,
):
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
