import dataclasses


@dataclasses.dataclass
class Vessel:
    volume: float = 250 / 1000


def kLa(rpm=1000, flow=25, vessel_volume=Vessel.volume):
    # flow in sccm (mL/min), vessel_volume in liter (convert to mL), kLA is in 1/h
    vvm = flow / (vessel_volume * 1000)

    # using Culture's 250mL vessel analysis
    return (
        -143.3
        + 0.24 * rpm
        + 84.3 * vvm
        + 0.035 * vvm * rpm
        + -4.72e-5 * rpm**2
        + -47.3 * vvm**2
    )
