import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint


def sir_model(y, t, N, beta, gamma):
    # The SIR model differential equations.
    S, I, R = y
    dSdt = -beta * S * I / N
    dIdt = (beta * S * I / N) - (gamma * I)
    dRdt = gamma * I
    return dSdt, dIdt, dRdt


def main():
    # Total population, N.
    N = 10000
    # Initial number of infected and recovered individuals, I0 and R0.
    I0, R0 = 1, 0
    # Everyone else, S0, is susceptible to infection initially.
    S0 = N - I0 - R0
    # Contact rate, b
    # eta, and mean recovery rate, gamma, (in 1/days).
    beta, gamma = 0.2, 1 / 10
    # A grid of time points (in days)
    weeks = 26
    t = np.linspace(0, weeks * 7, weeks * 7 * 10)

    # Initial conditions vector
    y0 = S0, I0, R0
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(sir_model, y0, t, args=(N, beta, gamma))
    S, I, R = ret.T

    # Plot the data on three separate curves for S(t), I(t) and R(t)
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)  # , axis_bgcolor='#dddddd', axisbelow=True)
    ax.set_title("$\\beta = {beta}$ / $\\gamma = {gamma}$".format(beta=beta, gamma=beta))
    ax.plot(t, S / N, 'b', alpha=0.5, lw=2, label='Susceptible')
    ax.plot(t, I / N, 'r', alpha=0.5, lw=2, label='Infected')
    ax.plot(t, R / N, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Population [%]')
    ax.set_ylim(0, 1.2)
    # ax.yaxis.set_tick_params(length=2)
    # ax.xaxis.set_tick_params(length=2)
    # ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    ax.legend()
    # legend.get_frame().set_alpha(0.5)
    # for spine in ('top', 'right', 'bottom', 'left'):
    #     ax.spines[spine].set_visible(False)
    plt.show()


if __name__ == "__main__":
    main()
