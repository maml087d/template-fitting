import numpy as np
import json
import matplotlib.pyplot as plt
from numba.cuda import shared

import ang_coeff_calc as ang
from ptbinnedfittinghelpers import *


def main():
    # plotting templates
    ptbinned = False
    cuts = True
    bin = 0         #bins \in [0,nptbins-1] or 'of'

    if cuts:
        _cuts = ''
    else:
        _cuts = 'no'

    if ptbinned:
        pathtotemps =  '/home/maltem/Nextcloud/TU_Dresden/Bachelorthesis/minuitfitting/templates/ptbinned/Zjj_ptbin_temps_' + _cuts + 'cuts0'
        pathtocoeffs = '/home/maltem/Nextcloud/TU_Dresden/Bachelorthesis/minuitfitting/data_ptbinned/histo_nocuts'
        pathtoptbins = '/home/maltem/Nextcloud/TU_Dresden/Bachelorthesis/minuitfitting/data_ptbinned/ptbins.txt'
        pathtodata = '/home/maltem/Nextcloud/TU_Dresden/Bachelorthesis/minuitfitting/data_ptbinned/histo_' +_cuts + 'cuts'
        ptbins = loadptbins(pathtoptbins)
        nptbins = len(ptbins) - 1


    else:
        pathtotemps = '/home/maltem/Nextcloud/TU_Dresden/Bachelorthesis/minuitfitting/templates/inclusive/Zjj_inclusive_temps_' + _cuts + 'cuts1'
        pathtocoeffs = '/home/maltem/Nextcloud/TU_Dresden/Bachelorthesis/minuitfitting/data_inclusive/histo_nocuts'
        pathtodata = '/home/maltem/Nextcloud/TU_Dresden/Bachelorthesis/minuitfitting/data_inclusive/histo_' +_cuts + 'cuts'
        nptbins = 1
    angbins = np.array(loadbins(pathtodata))
    if bin == 'of':
        print(nptbins)
        _, temps = loadtemplates(pathtotemps, nptbins, _load='1D', overflow=True)
        _, coeff = loadcoeffs(pathtocoeffs, nptbins, overflow=True)
        _, data_of = loaddata(pathtodata, nptbins, _load='1D', overflow=True)
        sigma_unpol = data_of['xsec'][0]
        temps_cos = np.array(temps['cos'])
        temps_phi = np.array(temps['phi'])
    elif isinstance(bin, int):
        temps = loadtemplates(pathtotemps, nptbins, _load='1D', overflow=False)
        print(np.array(temps['cos']).shape)
        coeff = loadcoeffs(pathtocoeffs, nptbins, overflow=False)[bin, :]
        data = loaddata(pathtodata, nptbins, _load='1D', overflow=False)
        sigma_unpol = data['xsec'][bin][0]
        temps_cos = np.array(temps['cos'])[bin, :]
        temps_phi = np.array(temps['phi'])[bin, :]
    elif is_iterable(bin):
        raise NotImplementedError

    if not ptbinned:
        cos_errs = []
        phi_errs = []
        for i in range(9):
            with open(pathtotemps+'/ptbin0_cos'+str(i)+'.json', 'r') as f:
                cos_errs.append(json.load(f)["yerrs"])
            with open(pathtotemps+'/ptbin0_phi'+str(i)+'.json', 'r') as f:
                phi_errs.append(json.load(f)["yerrs"])
        cos_errs = np.array(cos_errs)
        phi_errs = np.array(phi_errs)
        print(cos_errs.shape, phi_errs.shape)
        print("_--------------------------------------_________")

    # filenames = []
    # for i in range(9):
    #     filenames.append('temps_cos' + str(i) + '.json')
    #     filenames.append('temps_phi' + str(i) + '.json')
    #
    # temps = {}
    # for f in filenames:
    #     with open(path_todata + f, 'r') as g:
    #         temps[f] = json.load(g)
    #
    # coeff = []
    # with open('coeffs.txt', 'r') as g:
    #     for i in range(8):
    #         coeff.append(float(g.readline()))
    #
    # with open('data/data_all/xsec_bin.json', 'r') as f:
    #     xsec_bin = json.load(f)

    # sigma_unpol = xsec_bin["yvals"][0]

    plot_phi = np.linspace(0, 2*np.pi, 1000)
    plot_cos = np.linspace(-1, 1, 1000)

    cos_dist = ang.diffcross_costheta(sigma_unpol, plot_cos, coeff)
    phi_dist = ang.diffcross_phi(sigma_unpol, plot_phi, coeff)
    with plt.rc_context({
        'font.size': 12,  # Base font size
        'axes.labelsize': 12,  # Axis label font size
        'axes.titlesize': 14,  # Title font size
        'xtick.labelsize': 10,  # X-tick label font size
        'ytick.labelsize': 10,  # Y-tick label font size
        'legend.fontsize': 10  # Legend font size
    }):
    # for fsdfsf in [0]:
        fig1, axs1 = plt.subplots(9, 2, figsize=(18,24))
        # fig1.suptitle('Templates without cuts', fontsize=25, fontweight='bold')
        # axs1 = [fig1.add_subplot(9, 2, 1), fig1.add_subplot(9, 2, 2),
        #         fig1.add_subplot(9, 2, 3), fig1.add_subplot(9, 2, 4),
        #         fig1.add_subplot(9, 2, 5), fig1.add_subplot(9, 2, 6),
        #         fig1.add_subplot(9, 2, 7), fig1.add_subplot(9, 2, 8),
        #         fig1.add_subplot(9, 2, 9), fig1.add_subplot(9, 2, 10),
        #         fig1.add_subplot(9, 2, 11), fig1.add_subplot(9, 2, 12),
        #         fig1.add_subplot(9, 2, 13), fig1.add_subplot(9, 2, 14),
        #         fig1.add_subplot(9, 2, 15), fig1.add_subplot(9, 2, 16),
        #         fig1.add_subplot(9, 2, 17), fig1.add_subplot(9, 2, 18)]

        axs1[0, 0].set_title('$\cos \\theta$')
        axs1[0, 1].set_title('$\phi$')
        for i in range(9):
            axs1[i, 0].set_ylabel(f'template $P_{i}$')
            axs1[i, 1].set_ylabel(f'template $P_{i}$')
            print(angbins.shape)
            dx = angbins[:,1] - angbins[:,0]
            x = (angbins[:,:-1].T + dx / 2).T
            print(x.shape)
            print(x)
            # axs1[i, 0].stairs(temps_cos[i, :], angbins[0, :], label="template")
            # axs1[i, 1].stairs(temps_phi[i, :], angbins[1, :], label="template")
            axs1[i, 0].errorbar(x[0,:], temps_cos[i, :], xerr=dx[0]/2, yerr=cos_errs[i, :], fmt='o', label="template")
            axs1[i, 1].errorbar(x[1,:], temps_phi[i, :], xerr=dx[1]/2, yerr=phi_errs[i, :],fmt='o', label="template")


        axs1[0, 0].plot(plot_cos, cos_dist[0], ls='--',label="theoretical Polynomial")
        axs1[1, 0].axhline(y=0, color='grey', linestyle='--')
        axs1[2, 0].axhline(y=0, color='grey', linestyle='--')
        axs1[3, 0].axhline(y=0, color='grey', linestyle='--')
        axs1[4, 0].plot(plot_cos, cos_dist[1], ls='--',label="theoretical Polynomial")
        axs1[5, 0].axhline(y=0, color='grey', linestyle='--')
        axs1[6, 0].axhline(y=0, color='grey', linestyle='--')
        axs1[7, 0].axhline(y=0, color='grey', linestyle='--')
        axs1[8, 0].plot(plot_cos, cos_dist[2], ls='--',label="theoretical Polynomial")

        axs1[0, 1].axhline(y=0, color='grey', linestyle='--')
        axs1[1, 1].axhline(y=0, color='grey', linestyle='--')
        axs1[2, 1].plot(plot_phi, phi_dist[0], ls='--',label="theoretical Polynomial")
        axs1[3, 1].plot(plot_phi, phi_dist[1], ls='--',label="theoretical Polynomial")
        axs1[4, 1].axhline(y=0, color='grey', linestyle='--')
        axs1[5, 1].plot(plot_phi, phi_dist[2], ls='--',label="theoretical Polynomial")
        axs1[6, 1].axhline(y=0, color='grey', linestyle='--')
        axs1[7, 1].plot(plot_phi, phi_dist[3], ls='--',label="theoretical Polynomial")
        axs1[8, 1].plot(plot_phi, phi_dist[4], ls='--',label="theoretical Polynomial")

        axs1[8, 0].set_xlabel('$\cos \\theta$')
        axs1[8, 1].set_xlabel('$\phi$')

        for axs in axs1:
            for ax in axs:
                ax.legend()
        print(f'save fig as: templates_inclusive_{"cuts" if cuts else "nocuts"}.png')
        fig1.savefig(f'templates_inclusive_{"cuts" if cuts else "nocuts"}.png', dpi=300, transparent=True, bbox_inches='tight')
        plt.show()

if __name__ == "__main__":
    main()