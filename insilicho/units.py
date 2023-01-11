import typing

import pint

# Setup the pint unit-tracking system.
# All units must come from same registry.
UNIT = pint.UnitRegistry()

UnitType = typing.Any  # not sure
# Some convenience types
temperature_t = UNIT("degC")
concentration_t = UNIT("mmol / L")
volume_t = UNIT("mL")
density_t = UNIT("mg / L")
time_t = UNIT("seconds")
dimensionless_t = UNIT("")
current_t = UNIT("ampere")
power_t = UNIT("watt")
mass_t = UNIT("grams")
pump_rate_t = UNIT("mL / min")
gas_rate_t = UNIT("mL / min")
molar_mass_t = UNIT("mg / mmol")
agitation_t = UNIT("revolutions_per_minute")
# This one is a special case
str_t = object()
mode_t = object()  # can be bools or ints
list_t = object()
dict_t = object()
ID_t = object()
same_as_input_t = object()
NO_UNITS = [str_t, mode_t, list_t, dict_t, ID_t, same_as_input_t]
