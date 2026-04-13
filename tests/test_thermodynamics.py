import pytest

from insilicho.chemistry import Constants, Thermodynamics


class TestPH:
    def test_neutral_ph(self):
        # 1e-4 mM = 1e-7 M = pH 7
        assert Thermodynamics.pH(1e-4) == pytest.approx(7.0)

    def test_acidic_ph(self):
        # 1 mM = 1e-3 M = pH 3
        assert Thermodynamics.pH(1) == pytest.approx(3.0)

    def test_basic_via_oh(self):
        # 1 mM OH- = 1e-3 M => pOH = 3 => pH = 11
        assert Thermodynamics.pH(1, use_OH=True) == pytest.approx(11.0)


class TestHenrysCoeff:
    def test_o2_at_25c(self):
        h = Thermodynamics.HenrysCoeff(25, gas="O2")
        # Henry's coeff for O2 at 25C is ~roughly 4.3e4 atm
        assert 3e4 < h < 6e4

    def test_co2_at_25c(self):
        h = Thermodynamics.HenrysCoeff(25, gas="CO2")
        # CO2 is more soluble, lower Henry's constant
        assert h < Thermodynamics.HenrysCoeff(25, gas="O2")

    def test_increases_with_temperature(self):
        # Gas solubility decreases with temp => Henry's coeff increases
        assert Thermodynamics.HenrysCoeff(37) > Thermodynamics.HenrysCoeff(25)


class TestCsatOxygen:
    def test_physiological_temp(self):
        csat = Thermodynamics.Csat_oxygen(37)
        # Dissolved O2 saturation at 37C is ~0.2 mM
        assert 0.1 < csat < 0.4

    def test_decreases_with_temperature(self):
        assert Thermodynamics.Csat_oxygen(25) > Thermodynamics.Csat_oxygen(37)


class TestPressureConcentration:
    def test_roundtrip(self):
        T = 37
        c_original = 0.21  # mmol/L
        p = Thermodynamics.c_to_p(c_original, T)
        c_back = Thermodynamics.p_to_c(p, T)
        assert c_back == pytest.approx(c_original)

    def test_ideal_gas_law(self):
        # p = c * R * T(K), verify against manual calc
        T = 25  # degC
        c = 1.0  # mmol/L
        p = Thermodynamics.c_to_p(c, T)
        expected = c * Constants.Rg * (T + 273.15)
        assert p == pytest.approx(expected)
