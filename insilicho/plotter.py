import numpy as np
from matplotlib import pyplot as plt


def plot(tspan: np.ndarray, state: np.ndarray, state_vars: np.ndarray) -> plt.figure:
    """Default plot for solver

    Returns:
        plt.figure: A figure object.
    """
    plt.rcParams["figure.figsize"] = [16, 12]
    fig = plt.figure()
    ax = fig.add_subplot(111)  # The big subplot
    ax1 = fig.add_subplot(811)
    ax2 = fig.add_subplot(812)
    ax3 = fig.add_subplot(813)
    ax4 = fig.add_subplot(814)
    ax5 = fig.add_subplot(815)
    ax6 = fig.add_subplot(816)
    ax7 = fig.add_subplot(817)
    ax8 = fig.add_subplot(818)

    # Turn off axis lines and ticks of the big subplot
    ax.spines["top"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["right"].set_color("none")
    ax.tick_params(labelcolor="w", top=False, bottom=False, left=False, right=False)
    ax.set_xlabel("Time [hrs]")
    ax1.set_ylabel("Viable cells [millions/mL]")
    ax2.set_ylabel("Oxygen [mM]")
    ax3.set_ylabel("Osmolarity [mM]")
    ax4.set_ylabel("mAbs [mg/L]")
    ax5.set_ylabel("pH [-]")
    ax6.set_ylabel("Volume [L]")
    ax7.set_ylabel("Species [mM]")
    ax8.set_ylabel("Temp[degC]")

    (Xv, Xt, Cglc, Cgln, Clac, Camm, Cmab, Coxygen, V, pH) = state.transpose()
    Temp = state_vars[:, 1]
    Osmolarity = state_vars[:, 9]

    ax1.plot(tspan, Xv * 1e-9)  # converted to millions/mL by *1e-9
    ax2.plot(tspan, Coxygen)
    ax3.plot(tspan, Osmolarity)
    ax4.plot(tspan, Cmab)
    ax5.plot(tspan, pH)  # pH
    ax6.plot(tspan, V)
    ax7.plot(tspan, Cglc, tspan, Cgln, tspan, Clac, tspan, Camm)
    ax7.legend(["Glc", "Gln", "Lac", "Amm"])
    ax8.plot(tspan, Temp)
    plt.show()

    return fig
