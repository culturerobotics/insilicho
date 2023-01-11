import numpy as np
from matplotlib import pyplot as plt

from insilicho import growth_model, parameters
from insilicho.solver import solve


def solve_and_plot(
    params=None,
    ic=None,
    model=growth_model.model,
    tspan=np.linspace(0, 24 * 10, 1000),
    feed_fn=None,
    temp_fn=None,
    solver_hmax=np.inf,
):
    if params is None:
        params = parameters.InputParameters()
    if ic is None:
        ic = parameters.InitialConditions()

    def default_plot_setup(tspan, state, state_vars):
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

    state, state_vars, info = solve(
        params,
        ic,
        model,
        tspan=tspan,
        feed_fn=feed_fn,
        temp_fn=temp_fn,
        solver_hmax=solver_hmax,
    )
    default_plot_setup(tspan, state, state_vars)
    return state, state_vars, info
