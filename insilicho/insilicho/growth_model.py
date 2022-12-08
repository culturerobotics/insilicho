# mypy: disable-error-code=operator


import numpy as np

from insilicho import chemistry, parameters
from insilicho.chemistry import Species
from insilicho.chemistry import Thermodynamics


class ProcessDependencies:
    @staticmethod
    def GrowthRatepHSensitivity(pH, pHoptimal=7.05, deltapH_sensitivity=0.89):
        return np.exp(-((pH - pHoptimal) ** 2) / deltapH_sensitivity**2)

    @staticmethod
    def OsmolarityDependence(Osmolarity):
        isotonic_osmolality = 320  # mOSm/kg
        critical_osmolality = 400  # beyond this no production of mAbs
        return np.exp(
            -1
            * max(0, Thermodynamics.molar_to_molal(Osmolarity, 1) - isotonic_osmolality)
            / (critical_osmolality - isotonic_osmolality)
        )

    @staticmethod
    def GrowthRateOxygenDependence(Coxygen, K_oxygen=0.03):
        # K_oxygen is in mmol/liter, I have no idea what it should be so making it up
        return Coxygen / (Coxygen + K_oxygen)

    @staticmethod
    def OxygenConsumptionRate():
        # From https://www.pnas.org/doi/pdf/10.1073/pnas.79.4.1166
        oxygen_molecules_needed = 3.8e7  # O2_molecules/cell/second
        return (
            oxygen_molecules_needed / chemistry.Constants.Avogadro * 3600
        )  # millimolesO2/hr/cell

    @staticmethod
    def GrowthRateTempSensitivity(T, Toptimal=36.4, deltaT_sensitivity=1.0):
        # Temperature model for growth/death follows a gaussian around Topt
        # The idea and shape of curve is from Carcano thesis
        return np.exp(-((T - Toptimal) ** 2) / deltaT_sensitivity**2)

    @staticmethod
    def DeathRateTempSensitivity(T):
        """
        Death rate can have different sensitivity
        """
        return 1


def state_vars(
    t,
    state,
    params: parameters.InputParameters,
    feed_fn=None,
    temp_fn=None,
):

    # state
    Xv, Xt, Cglc, Cgln, Clac, Camm, Cmab, Coxygen, V, pH = state

    if not feed_fn or not temp_fn:
        raise IOError("feed/temp model missing")

    F = feed_fn(t, V)
    T = temp_fn(t)

    Osmolarity = (
        Cglc * Species.Glc.phi
        + Cgln * Species.Gln.phi
        + Clac * Species.Lac.phi
        + Camm * Species.NH3.phi
        + Cmab * 1
    )  # TODO: This is probably not right

    # mus
    mu = (
        params.mu_max
        * Cglc
        / (Cglc + params.Ks_glc)
        * Cgln
        / (Cgln + params.Ks_gln)
        * params.Ki_amm
        / (Camm + params.Ki_amm)
        * (
            ProcessDependencies.GrowthRateTempSensitivity(T)
            * ProcessDependencies.GrowthRateOxygenDependence(Coxygen)
            * ProcessDependencies.GrowthRatepHSensitivity(pH)
        )
    )

    mu_d = (
        params.mu_d_max
        * params.Ks_glc
        / (Cglc + params.Ks_glc)
        * params.Ks_gln
        / (Cgln + params.Ks_gln)
        * Camm
        / (Camm + params.Ki_amm)
        * ProcessDependencies.DeathRateTempSensitivity(T)
        + params.mu_d_min
    )

    # qs
    q_glc = (
        params.q_glc_max
        * Cglc
        / (Cglc + params.k_glc)
        * (mu / (mu + params.mu_max) + 0.5)
    )  # TODO: this probably needs to depend somehow on O2?
    q_gln = (
        params.q_gln_max * Cgln / (Cgln + params.k_gln)
    )  # TODO: this probably needs to depend somehow on O2?

    q_lac_uptake = 0.0
    if Cglc < 0.5:
        q_lac_uptake = float(
            params.q_lac_max
        )  # TODO: Error in this equation also causes lactate to go below 0 when starting glucose conc is low (~10 mmol/L)!
    q_lac = params.Y_lac_glc * Cglc / Clac * q_glc - q_lac_uptake
    q_amm = params.Y_amm_gln * q_gln  # nothing consumes ammonia once formed?
    q_mab = params.Y_mab_cell * ProcessDependencies.OsmolarityDependence(Osmolarity)
    if Cglc < 1 or Camm > params.Ki_amm or Cgln < 0.05:
        q_mab = q_mab / 100

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
    feed_fn,
    temp_fn,
):
    # params repacking
    params = parameters.InputParameters(*args)
    # state
    Xv, Xt, Cglc, Cgln, Clac, Camm, Cmab, Coxygen, V, pH = state
    F, T, mu, mu_d, q_glc, q_gln, q_lac, q_amm, q_mab, Osmolarity = state_vars(
        t, state, params, feed_fn, temp_fn
    )

    # diffeqs
    dXv = (mu - mu_d) * Xv
    dXt = mu * Xv - params.K_lys * (
        Xt - Xv
    )  # umm what does this do really? consider combining lysis with viable cell population
    dCglc = -q_glc * Xv + F * (params.Cglc_feed - Cglc) / V
    dCgln = -q_gln * Xv + F * (params.Cgln_feed - Cgln) / V
    dClac = q_lac * Xv - F * Clac / V
    dCamm = q_amm * Xv - F * Camm / V
    dCmab = q_mab * Xv - params.mab_time_decay * Cmab**1.2
    # TODO: Cmab at different temp goes to same number but at different times, need to recheck.

    # Process
    dCoxygen = 0
    dV = F
    dpH = 0  # TODO

    return [
        dXv,
        dXt,
        dCglc,
        dCgln,
        dClac * Clac / (Clac + parameters.EPSILON),
        dCamm * Camm / (Camm + parameters.EPSILON),
        dCmab,
        dCoxygen,
        dV,
        dpH,
    ]
