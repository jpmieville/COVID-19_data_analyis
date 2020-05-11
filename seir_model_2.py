import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint, cumtrapz


def beta_function(t, start, beta):
    if t <= start:
        return beta
    elif t > start:
        return ((beta - 0.01)  * np.exp(-0.05*(t - start))) + (0.01)

def deriv(y, t, N, beta, gamma, sigma, start):
    # The SEIR model differential equations.
    #
    S, E, I, R = y
    dSdt = -beta_function(t, start, beta) * S * (I + E) / N
    dEdt = (beta_function(t, start, beta) * S * (I + E) / N) - (sigma * E)
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
    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
    beta = 0.2
    gamma = 1 / 20
    sigma = 1 / 10
    print(beta / gamma)
    start = 25
    # A grid of time points (in days)
    weeks = 52
    t = np.linspace(0, weeks * 7, weeks * 7 * 10)

    # Initial conditions vector
    y0 = S0, E0, I0, R0
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(deriv, y0, t, args=(N, beta, gamma, sigma, start))
    S, E, I, R = ret.T

    # Plot the data on three separate curves for S(t), I(t) and R(t)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7))
    #
    # ax = fig.add_subplot(111)  # , axis_bgcolor='#dddddd', axisbelow=True)
    # infected =
    ax1.plot(t, S / N, 'b', label='Susceptible')
    ax1.plot(t, E / N, 'c', label='Exposed')
    ax1.plot(t, I / N, 'r', label='Infected')
    # ax.plot(t, I.cumsum() / N, 'r', alpha=0.5, lw=2, label='Infected cumulated')
    ax1.plot(t, R / N, 'g', label='Recovered with immunity')
    # ax.plot(t, (S + E + I + R) / N, 'y', alpha=0.5, lw=2, label='Total')
    ax1.set_title("$\\beta = {beta}$ / $\\gamma = {gamma}$ / $\\sigma = {sigma}$".format(beta=beta, gamma=gamma, sigma=sigma))
    ax1.set_xlabel('Time in days')
    ax1.set_ylabel('Relative population')
    ax1.set_ylim(0, 1.05)
    # ax.yaxis.set_tick_params(length=2)
    # ax.xaxis.set_tick_params(length=2)
    # ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    legend = ax1.legend()
    legend.get_frame().set_alpha(0.5)
    # for spine in ('top', 'right', 'bottom', 'left'):
    #     ax.spines[spine].set_visible(False)
    print(f"{int(t[np.argmax(E)])} max Susceptible")
    print(f"{int(t[np.argmax(I)])} max Infected")
    # plt.show()
    # I_cum = cumtrapz(I, t)
    # print(max(I_cum))
    # ax2.plot(t[:-1], I_cum / max(I_cum), 'b')
    ax2.plot(t, I / N, 'r')
    ax2.plot(t, np.array([beta_function(i, start, beta) for i in t]), "m")
    ax2.plot(t, R / N, 'g')
    # ax2.set_yscale('log')
    # ax2.set_xscale('log')
    # ax3.plot(R / N, S / N, )
    # ax3.set_yscale('log')
    plt.show()


if __name__ == "__main__":
    main()
