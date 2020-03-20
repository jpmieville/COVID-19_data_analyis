import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint


def deriv(y, t, N, beta, gamma, sigma):
    # The SEIR model differential equations.
    #
    S, E, I, R = y
    dSdt = -beta * S * I / N
    dEdt = (beta * S * I / N) - (sigma * E)
    dIdt = (sigma * E) - (gamma * I)
    dRdt = gamma * I
    return dSdt, dEdt, dIdt, dRdt


def main():
    # Total population, N.
    N = 10000
    # Initial number of infected and recovered individuals, I0 and R0.
    E0 = 1
    I0 = 0
    R0 = 0
    # Everyone else, S0, is susceptible to infection initially.
    S0 = N - -E0 - I0 - R0
    # Contact rate, b
    # eta, and mean recovery rate, gamma, (in 1/days).
    beta = 3
    gamma = 1 / 10
    sigma = 1 / 5
    print(beta / gamma)
    # A grid of time points (in days)
    weeks = 52
    t = np.linspace(0, weeks * 7, weeks * 7 * 10)

    # Initial conditions vector
    y0 = S0, E0, I0, R0
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(deriv, y0, t, args=(N, beta, gamma, sigma))
    S, E, I, R = ret.T

    # Plot the data on three separate curves for S(t), I(t) and R(t)
    fig, ax = plt.subplots()
    #
    # ax = fig.add_subplot(111)  # , axis_bgcolor='#dddddd', axisbelow=True)
    # infected =
    ax.plot(t, S / N, 'b', alpha=0.5, lw=2, label='Susceptible')
    ax.plot(t, E / N, 'c', alpha=0.5, lw=2, label='Exposed')
    ax.plot(t, I / N, 'r', alpha=0.5, lw=2, label='Infected')
    # ax.plot(t, I.cumsum() / N, 'r', alpha=0.5, lw=2, label='Infected cumulated')
    ax.plot(t, R / N, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
    # ax.plot(t, (S + E + I + R) / N, 'y', alpha=0.5, lw=2, label='Total')
    ax.set_title("$\\beta = {beta}$ / $\\gamma = {gamma}$ / $\\sigma = {sigma}$".format(beta=beta, gamma=gamma, sigma=sigma))
    ax.set_xlabel('Time in days')
    ax.set_ylabel('Relative population')
    ax.set_ylim(0, 1.05)
    # ax.yaxis.set_tick_params(length=2)
    # ax.xaxis.set_tick_params(length=2)
    # ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    # for spine in ('top', 'right', 'bottom', 'left'):
    #     ax.spines[spine].set_visible(False)
    print(f"{int(t[np.argmax(E)])} max Susceptible")
    print(f"{int(t[np.argmax(I)])} max Infected")
    plt.show()


if __name__ == "__main__":
    main()
