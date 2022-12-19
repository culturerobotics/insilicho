import dataclasses
import typing

import numpy as np

from insilicho import units


@dataclasses.dataclass
class Specie:
    formula: str
    molar_mass: typing.Union[float, units.UnitType]
    phi: typing.Union[int, units.UnitType] = 1
    charge: int = 0
    disassoc_const: float = 0  # M units TODO: get right values for these


class Species:
    # Some convenience species objects
    Glc = Specie("C6H12O6", 180)
    Gln = Specie("C5H10N2O3", 146.14)
    Asp = Specie("C4H8N2O3", 132.12)
    Lac = Specie("C3H6O3", 90.08, phi=2, disassoc_const=1.3e-1)
    NH3 = Specie("NH3", 17.031, disassoc_const=1.9e-2)
    O2 = Specie("O2", 31.999)
    CO2 = Specie("CO2", 44.01)


class Constants:
    Avogadro = 6.02214179e20
    Rg = 0.0821 * 1e3  # L*atm/mmol/K


class Thermodynamics:
    # TODO: resolve units in this class
    @staticmethod
    def pH(C, use_OH=False):
        """returns pH using free proton concentration in mM"""
        if use_OH:
            return 14 + np.log10(C / 1000)
        return -np.log10(C / 1000)

    @staticmethod
    def HenrysCoeff(T, gas=Species.O2.formula):
        """Henry's coefficient (atm) at temp T (in degC)"""
        C = {"O2": 1700, "CO2": 2400}  # atm
        HTs = {"O2": 4.259e4, "CO2": 0.163e4}  # K
        Ts = 289  # K
        return HTs[gas] * np.exp(-C[gas] * ((T + 273.15) ** -1 - Ts**-1))

    @staticmethod
    def Csat_oxygen(T):
        # in mmol/L; from https://www.jstor.org/stable/2835754
        T = T + 273.15  # Kelvin
        lnC = (
            -139.34410
            + (1.575701e5 / T)
            - (6.642308e7 / T**2)
            + (1.243800e10 / T**3)
            - (8.621949e11 / T**4)
        )
        return np.exp(lnC) / Species.O2.molar_mass

    @staticmethod
    def p_to_c(p, T):
        return p / Constants.Rg / (T + 273.15)  # Kelvin

    @staticmethod
    def c_to_p(c, T):
        return c * Constants.Rg * (T + 273.15)  # Kelvin

    @staticmethod
    def molar_to_molal(M, rho_sol=1):
        # TBD with units, right now it just returns M (dilute sol approx)
        return M / rho_sol
