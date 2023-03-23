import dataclasses
import typing

import numpy as np

from insilicho import units
from insilicho.chemistry import Thermodynamics

# simulation epsilon, a small number takes whatever units we want
EPSILON = float(np.finfo(float).eps)
# A SMALL CONC value to avoid zeros
SMALL_CONC = 0.01  # mM


class UnitValidationMixin:
    @staticmethod
    def units_map():
        raise NotImplementedError

    def __post_init__(self):
        for field in dataclasses.fields(self):
            setattr(self, field.name, getattr(self, field.name))

    def __setattr__(self, name, val):
        if isinstance(val, str):
            unitd_val = units.UNIT(val)
            to_units = self.units_map()[name]
            try:
                self.__dict__[name] = unitd_val.to(to_units).magnitude
            except units.pint.errors.DimensionalityError:
                raise ValueError(
                    f"Dimensionality error in setting {name}, cannot convert from:"
                    f"{unitd_val.units} to: {to_units}"
                )
        else:
            super().__setattr__(name, val)

    def tolist(self):
        return [getattr(self, field.name) for field in dataclasses.fields(self)]


@dataclasses.dataclass(order=True)
class InputParameters(UnitValidationMixin):
    # From SI, Table 2 of "Model uncertainty-based evaluation of process strategies
    # during scale-up of biopharmaceutical process"
    mu_max: typing.Union[float, str] = 0.043  # 1/hour
    mu_d_max: typing.Union[float, str] = 0.06  # 1/hour
    mu_d_min: typing.Union[float, str] = 0.001  # 1/hour

    k_glc: typing.Union[float, str] = 0.20  # mmol/liter
    k_gln: typing.Union[float, str] = 2.5  # mmol/liter
    K_lys: typing.Union[float, str] = 0.001  # 1/hour
    Ks_amm: typing.Union[float, str] = 10.0  # mmol/liter
    Ki_amm: typing.Union[float, str] = 10.0  # mmol/liter
    Ks_glc: typing.Union[float, str] = 0.02  # mmol/liter
    Ks_gln: typing.Union[float, str] = 0.03  # mmol/liter

    q_mab: typing.Union[float, str] = 3.12e-10  # mg/cell/hour
    q_glc_max: typing.Union[float, str] = 0.05e-9  # mmol/cell/hour
    q_gln_max: typing.Union[float, str] = 0.054e-9  # mmol/cell/hour
    q_lac_max: typing.Union[float, str] = 0.2e-9  # mmol/cell/hour

    # Yield parameters
    Y_amm_gln: typing.Union[float, str] = 0.9  # dimensionless (mol amm/ mol gln)
    Y_lac_glc: typing.Union[float, str] = 0.25  # dimensionless (mol lac/ mol glc)

    # Feed conditions
    Cglc_feed: typing.Union[float, str] = 150.0  # mmol/liter
    Cgln_feed: typing.Union[float, str] = 10.0  # mmol/liter

    # Optimal
    T_optimal: typing.Union[float, str] = 36.40
    T_optimal_decay_spread: typing.Union[float, str] = 3.12
    pH_optimal: typing.Union[float, str] = 6.99
    pH_optimal_decay_spread: typing.Union[float, str] = 1.00

    # Ndays to sim
    Ndays: int = 12  # days
    Nsamples: int = 2  # per day

    @staticmethod
    def units_map():
        return {
            "mu_max": "1/hr",
            "mu_d_max": "1/hr",
            "mu_d_min": "1/hr",
            "k_glc": "mmol/L",
            "k_gln": "mmol/L",
            "K_lys": "1/hr",
            "Ks_amm": "mmol/L",
            "Ki_amm": "mmol/L",
            "Ks_glc": "mmol/L",
            "Ks_gln": "mmol/L",
            "q_mab": "mg/hr",
            "q_glc_max": "mmol/hr",
            "q_gln_max": "mmol/hr",
            "q_lac_max": "mmol/hr",
            "Y_amm_gln": "dimensionless",
            "Y_lac_glc": "dimensionless",
            "Cglc_feed": "mmol/L",
            "Cgln_feed": "mmol/L",
            "T_optimal": "degC",
            "T_optimal_decay_spread": "degC",
            "pH_optimal": "dimensionless",
            "pH_optimal_decay_spread": "dimensionless",
        }


@dataclasses.dataclass
class InitialConditions(UnitValidationMixin):
    Xv: float = 3e9
    Xt: float = Xv
    Cglc: float = 150
    Cgln: float = 10
    Clac: float = EPSILON
    Camm: float = EPSILON
    Cmab: float = EPSILON
    Coxygen: float = Thermodynamics.Csat_oxygen(35)
    V: float = 40 / 1000  # liter
    pH: float = 7.0

    @staticmethod
    def units_map():
        return {
            "V": "L",
            "Cglc": "mmol/L",
            "Cgln": "mmol/L",
            "Clac": "mmol/L",
            "Camm": "mmol/L",
            "Cmab": "mmol/L",
            "Coxygen": "mmol/L",
            "pH": "dimensionless",
            "Xv": "1/L",
            "Xt": "1/L",
        }
