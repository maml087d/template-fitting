import numpy as np
import matplotlib.pyplot as plt
import json
from ptbinnedfittinghelpers import loadcoeffs

def A0(cos_theta):
    P0 = np.mean((1 - 3 * cos_theta**2)/2)
    return P0 * 20/3 + 2/3


def A1(theta, phi):
    P1 = np.mean(np.sin(2 * theta) * np.cos(phi))
    return P1 * 5


def A2(theta, phi):
    P2 = np.mean(np.sin(theta)**2 * np.cos(2 * phi))
    return 10 * P2


def A3(theta, phi):
    P3 = np.mean(np.sin(theta) * np.cos(phi))
    return 4 * P3


def A4(cos_theta):
    return 4 * np.mean(cos_theta)


def A5(theta, phi):
    P5 = np.mean(np.sin(theta)**2 * np.sin(2*phi))
    return 5 * P5


def A6(theta, phi):
    P6 = np.mean(np.sin(2 * theta) * np.sin(phi))
    return 5 * P6


def A7(theta, phi):
    P7 = np.mean(np.sin(theta) * np.sin(phi))
    return 4 * P7


def calc_coef(theta, phi):
    cos_theta = np.cos(theta)
    a0 = A0(cos_theta)
    a1 = A1(theta, phi)
    a2 = A2(theta, phi)
    a3 = A3(theta, phi)
    a4 = A4(cos_theta)
    a5 = A5(theta, phi)
    a6 = A6(theta, phi)
    a7 = A7(theta, phi)

    return a0, a1, a2, a3, a4, a5, a6, a7


def diffcross_thetaphi(norm, cos_theta, phi, A0, A1, A2, A3, A4, A5, A6, A7):
    theta = np.arccos(cos_theta)
    P8 = 1 + cos_theta**2
    P0 = (A0 * (1 - 3 * cos_theta**2))/2
    P1 = A1 * np.sin(2 * theta) * np.cos(phi)
    P2 = A2 * np.sin(theta)**2 * np.cos(2 * phi) / 2
    P3 = A3 * np.sin(theta) * np.cos(phi)
    P4 = A4 * cos_theta
    P5 = A5 * np.sin(theta)**2 * np.sin(2 * phi)
    P6 = A6 * np.sin(2 * theta) * np.sin(phi)
    P7 = A7 * np.sin(theta) * np.sin(phi)
    return norm * (P0 + P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8)


def diffcross_costheta(norm, cos_theta, A):
    p0 = A[0]/2 * (1 - 3 * cos_theta**2)
    p4 = A[4] * cos_theta
    p8 = 1 + cos_theta**2
    p = norm * np.array([p0, p4, p8])
    return 3 * p / 8


def diffcross_phi(norm, phi, A):
    p8 = np.ones_like(phi)
    p2 = A[2]/4 * np.cos(2 * phi)
    p3 = 3 * np.pi * A[3] / 16 * np.cos(phi)
    p5 = A[5] / 2 * np.sin(2 * phi)
    p7 = 3 * np.pi * A[7] / 16 * np.sin(phi)
    p = norm/(2 * np.pi) * np.array((p2, p3, p5, p7, p8))
    return p


def main():
    # getting data from json
    prefix="datbin0"
    cuts='no' # '' or 'no'
    # with open('histodata/data.json', 'r') as f:
    #     data = json.load(f)

    with open('data_inclusive/histo_' + cuts + 'cuts/' + prefix + '_cos.json', 'r') as f:
        hist_costheta = json.load(f)

    with open('data_inclusive/histo_' + cuts + 'cuts/' + prefix + '_phi.json', 'r') as f:
        hist_phi = json.load(f)

    with open('data_inclusive/histo_' + cuts + 'cuts/' + prefix + '_xsec.json', 'r') as f:
        xsec_bin = json.load(f)

    for k in hist_costheta.keys():
        hist_costheta[k] = np.array(hist_costheta[k])

    for k in hist_phi.keys():
        hist_phi[k] = np.array(hist_phi[k])

    for k in xsec_bin.keys():
        xsec_bin[k] = np.array(xsec_bin[k])

    # data = np.array(data[1:][:])
    #
    # theta = data[:, 0]
    # phi = data[:, 1]

    sigma_unpol = xsec_bin["yvals"][0]

    # calc and print angular coeffs
    # coeff = calc_coef(theta, phi)
    # coeff_out = ""
    # for i in range(len(coeff)):
    #     print(f"A{i} = {coeff[i]}")
    #     coeff_out += f"{coeff[i]}\n"
    # coeff_out += f"{sigma_unpol}"
    #
    # with open("coeffs.txt", "w") as f:
    #     f.write(coeff_out)

    coeff = loadcoeffs("data_inclusive/histo_nocuts", 1)[0][:-1]
    print(coeff)

    plot_cos = np.linspace(-1, 1, 1000)
    plot_phi = np.linspace(0, 2 * np.pi, 1000)

    cos_dist = diffcross_costheta(sigma_unpol, plot_cos, coeff)
    phi_dist = diffcross_phi(sigma_unpol, plot_phi, coeff)

    sum_cos = np.cumsum(cos_dist, axis=0)
    sum_phi = np.cumsum(phi_dist, axis=0)

    ind_cos = [0, 4, 8]
    ind_phi = [2, 3, 5, 7, 8]

    p8_phi_scale = 0.1

    # plotting
    with plt.rc_context({
        'font.size': 12,  # Base font size
        'axes.labelsize': 12,  # Axis label font size
        'axes.titlesize': 14,  # Title font size
        'xtick.labelsize': 10,  # X-tick label font size
        'ytick.labelsize': 10,  # Y-tick label font size
        'legend.fontsize': 10  # Legend font size
    }):
        fig = plt.figure(figsize=(12, 8))
        axs = []
        axs.append(fig.add_subplot(3, 2, 1))
        axs.append(fig.add_subplot(3, 2, 2))
        axs.append(fig.add_subplot(3, 2, 3))
        axs.append(fig.add_subplot(3, 2, 4))
        axs.append(fig.add_subplot(3, 2, 5))
        axs.append(fig.add_subplot(3, 2, 6))

        axs[0].set_title('cos($\\theta$)')
        axs[1].set_title('$\phi$')

        # axs[0].set_xlabel('cos($\\theta$)')
        # axs[2].set_xlabel('cos($\\theta$)')
        axs[4].set_xlabel('cos($\\theta$)')
        # axs[1].set_xlabel('$\phi$')
        # axs[3].set_xlabel('$\phi$')
        axs[5].set_xlabel('$\phi$')

        for ax in axs[0:2]:
            ax.set_ylabel('Polynimial value', fontsize=12)

        axs[2].set_ylabel('d$\sigma$ / dcos($\\theta$)  [fb]')
        axs[3].set_ylabel('d$\sigma$ / d$\phi$  [fb]')

        axs[4].set_ylabel('d$\sigma$ / dcos($\\theta$)  [fb]')
        axs[5].set_ylabel('d$\sigma$ / d$\phi$  [fb]')

        for i in range(len(cos_dist[:, 0])):
            axs[0].plot(plot_cos, cos_dist[i, :], label=f"$P_{ind_cos[i]}$")

        for i in range(len(phi_dist[:, 0])-1):
            axs[1].plot(plot_phi, phi_dist[i, :], label=f"$P_{ind_phi[i]}$")
        axs[1].plot(plot_phi, phi_dist[-1, :] * p8_phi_scale, label=f"$P_8$ * {p8_phi_scale:.1f}")

        for i in range(len(sum_cos[:, 0])):
            label = f'$P_{ind_cos[0]}$'
            for j in range(i):
                label += f"+$P_{ind_cos[j+1]}$"

            axs[2].plot(plot_cos, sum_cos[i, :], label=label)

        for i in range(len(sum_phi[:, 0])):
            label = f'$P_{ind_phi[0]}$'
            for j in range(i):
                label += f"+$P_{ind_phi[j + 1]}$"
            axs[3].plot(plot_phi, sum_phi[i, :], label=label)

        axs[4].stairs(hist_costheta["yvals"], hist_costheta["bins"], label="Histogram from MC sample")
        axs[5].stairs(hist_phi["yvals"], hist_phi["bins"], label="Histogram from MC sample")

        axs[4].plot(plot_cos, sum_cos[-1, :], label="analytical distribution from angular coefficients")
        axs[5].plot(plot_phi, sum_phi[-1, :], label="analytical distribution from angular coefficients")

        for ax in axs:
            ax.legend()
        fig.savefig("/home/maltem/Documents/bachlorthesis/polynoms.png", dpi=300, transparent=True , bbox_inches='tight')
        plt.show()


if __name__ == "__main__":
    main()