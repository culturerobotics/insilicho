import dataclasses
import typing

from insilicho import units
from insilicho.chemistry import Thermodynamics

# simulation epsilon, a small number takes whatever units we want
EPSILON = 1e-12


@dataclasses.dataclass
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
                    f"Dimensionality error in setting {name}, cannot convert from: {unitd_val.units} to: {to_units}"
                )
        else:
            super().__setattr__(name, val)


@dataclasses.dataclass(order=True)
class InputParameters(UnitValidationMixin):
    # defaults from Table 3, Pg 241 of book "Animal Cell Biotechnology"
    # parameter values
    mu_max: typing.Union[float, str] = 0.20  # 1/hour
    mu_d_max: typing.Union[float, str] = 0.03  # 1/hour
    mu_d_min: typing.Union[float, str] = 0.003  # 1/hour

    k_glc: typing.Union[float, str] = 0.19  # mmol/liter
    k_gln: typing.Union[float, str] = 1  # mmol/liter
    K_lys: typing.Union[float, str] = 0.005  # 1/hour
    Ks_amm: typing.Union[float, str] = 10  # mmol/liter
    Ki_amm: typing.Union[float, str] = 10  # mmol/liter
    Ks_glc: typing.Union[float, str] = 0.03  # mmol/liter
    Ks_gln: typing.Union[float, str] = 0.03  # mmol/liter

    q_mab: typing.Union[float, str] = 0.3e-9  # mmol/liter/hour
    q_glc_max: typing.Union[float, str] = 0.25e-9  # mmol/liter/hour
    q_gln_max: typing.Union[float, str] = 0.085e-9  # mmol/liter/hour
    q_lac_max: typing.Union[float, str] = 0.1e-9  # mmol/liter/hour
    mab_time_decay: typing.Union[float, str] = 1e-3  # 1/hour/[C^0.2 units]

    Cglc_feed: typing.Union[float, str] = 150  # mmol/liter
    Cgln_feed: typing.Union[float, str] = 150  # mmol/liter

    Y_amm_gln: float = 2  # -
    Y_lac_glc: float = 2  # -
    Y_mab_cell: float = 1e-9  # - value is arbitrary, doesnt work properly as this number is not defined in paper

    perfect_control: bool = True
    # # Bolus feeding
    # bolus_size: float = 0.03  # L/h
    # num_bolus: int = 10  # -
    # bolus_frequency: float = 24  # hrs

    # Ndays to sim
    Ndays: int = 12  # days
    Nsamples: int = 1  # per day

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
            "q_mab": "mmol/L/hr",
            "q_glc_max": "mmol/L/hr",
            "q_gln_max": "mmol/L/hr",
            "q_lac_max": "mmol/L/hr",
            "Cglc_feed": "mmol/L",
            "Cgln_feed": "mmol/L",
            "mab_time_decay": "1/(hr*(mmol/L)^.2)",
        }


@dataclasses.dataclass
class InitialConditions(UnitValidationMixin):
    V: float = 50 / 1000  # liter
    Xv: float = 8e9
    Xt: float = Xv
    Cglc: float = 100
    Cgln: float = 100
    Clac: float = EPSILON
    Camm: float = EPSILON
    Cmab: float = EPSILON
    Coxygen: float = Thermodynamics.Csat_oxygen(35)
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
