import typing

import pint

# Setup the pint unit-tracking system.
# All units must come from same registry.
UNIT = pint.UnitRegistry()

UnitType = typing.Any
