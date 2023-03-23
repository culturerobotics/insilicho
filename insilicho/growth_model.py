# mypy: disable-error-code=operator

import typing

import numpy as np

from insilicho import parameters
from insilicho.chemistry import Species

FeedFunctionType = typing.Callable[[float], float]
TempFunctionType = typing.Callable[[float], float]


def exponential_dependence_around_optima(
    x: float, optima: float, spread: float = 1.0
) -> float:
    """Function that decays both side of an optimal value.

    Args:
        x (float): value of a variable.
        optima (float): Optimal value for a variable.
        spread (float, optional): Spread around the optimal, sets rate of decay.
            Defaults to 1.0.

    Returns:
        float: Number between 0 and 1 indicating distance from optima.
    """
    return np.exp(-((x - optima) ** 2.0) / spread**2.0)


def state_vars(
    t: float,
    state: np.ndarray,
    params: parameters.InputParameters,
    feed_fn: typing.Optional[FeedFunctionType] = None,
    temp_fn: typing.Optional[TempFunctionType] = None,
):
    """Variables related to state (species and condition) but are not solved for in the
    ODE system.

    Args:
        t (float): Time point to evaluate on.
        state (np.ndarray): Current states (e.g. Xv, Cglc). Must follow order specified
            in parameters.InitialConditions.
        params (parameters.InputParameters): Parameters of GrowCHO model.
        feed_fn (typing.Optional[FeedFunctionType], optional): Callable describing feed
            profile. Defaults to None.
        temp_fn (typing.Optional[TempFunctionType], optional): Callable describing temp
            profile. Defaults to None.

    Raises:
        ValueError: Raised if `feed_fn` or `temp_fn` are not provided.

    Returns:
        typing.Tuple[float]: Returns a tuple of floats of the following variables:
            - F
            - T
            - mu
            - mu_d
            - q_glc
            - q_gln
            - q_lac
            - q_amm
            - q_mab
            - Osmolarity
    """

    # state
    Xv, Xt, Cglc, Cgln, Clac, Camm, Cmab, Coxygen, V, pH = state

    if not feed_fn or not temp_fn:
        raise ValueError("feed/temp model missing")

    F = feed_fn(t)
    T = temp_fn(t)

    Osmolarity = (
        Cglc * Species.Glc.phi
        + Cgln * Species.Gln.phi
        + Clac * Species.Lac.phi
        + Camm * Species.NH3.phi
    )

    # mus
    mu = (
        (
            params.mu_max
            * Cglc
            / (Cglc + params.Ks_glc)
            * Cgln
            / (Cgln + params.Ks_gln)
            * params.Ki_amm
            / (Camm + params.Ki_amm)
        )
        * exponential_dependence_around_optima(
            T, params.T_optimal, params.T_optimal_decay_spread  # type: ignore[arg-type]
        )  # This comes from Carcano et al.
        * exponential_dependence_around_optima(
            pH, params.pH_optimal, params.pH_optimal_decay_spread  # type: ignore[arg-type]
        )  # This is arbitrary
    )
    # TODO: make these optima and spread parameters

    mu_d = params.mu_d_min + (
        params.mu_d_max
        * params.Ks_glc
        / (Cglc + params.Ks_glc)
        * params.Ks_gln
        / (Cgln + params.Ks_gln)
        * Camm
        / (Camm + params.Ki_amm)
    )

    # qs
    q_glc = (
        params.q_glc_max
        * Cglc
        / (Cglc + params.k_glc)
        * (mu / (mu + params.mu_max) + 0.5)
    )
    q_gln = params.q_gln_max * Cgln / (Cgln + params.k_gln)

    if Cglc < 0.5:
        q_lac_uptake = float(params.q_lac_max)
    else:
        q_lac_uptake = 0.0
    q_lac = (
        params.Y_lac_glc * Cglc / (Clac + parameters.SMALL_CONC) * q_glc - q_lac_uptake
    )
    q_amm = params.Y_amm_gln * q_gln
    q_mab = params.q_mab
    if Camm > params.Ki_amm:
        q_mab = 0

    return (
        F,
        T,
        mu,
        mu_d,
        q_glc,
        q_gln,
        q_lac,
        q_amm,
        q_mab,
        Osmolarity,
    )


def model(
    t,
    state,
    args,
    feed_fn: typing.Optional[FeedFunctionType],
    temp_fn: typing.Optional[TempFunctionType],
):
    # params repacking
    params = parameters.InputParameters(*args)
    # state
    Xv, Xt, Cglc, Cgln, Clac, Camm, Cmab, Coxygen, V, pH = state
    F, T, mu, mu_d, q_glc, q_gln, q_lac, q_amm, q_mab, Osmolarity = state_vars(
        t, state, params, feed_fn, temp_fn
    )

    # diffeqs
    dXv = (mu - mu_d - F / V) * Xv
    dXt = mu * Xv - params.K_lys * (Xt - Xv) - F / V * Xt
    dCglc = -q_glc * Xv + F * (params.Cglc_feed - Cglc) / V
    dCgln = -q_gln * Xv + F * (params.Cgln_feed - Cgln) / V
    dClac = q_lac * Xv - F * Clac / V
    dCamm = q_amm * Xv - F * Camm / V
    dCmab = q_mab * Xv - F * Cmab / V

    # Process
    dCoxygen = 0
    dV = F
    dpH = 0

    return [
        dXv,
        dXt,
        dCglc,
        dCgln,
        dClac,
        dCamm,
        dCmab,
        dCoxygen,
        dV,
        dpH,
    ]
