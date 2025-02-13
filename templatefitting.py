# %config InlineBackend.figure_formats = ['svg']
from iminuit import Minuit
import numpy as np
import matplotlib.pyplot as plt
# import json

class TemplateFits_1set():

    def __init__(self, xe, templates_cos, data_cos, coeffs=np.ones(9), lumi=140, _2d=False):
        if not _2d:
            self.xe = np.array(xe)
            self.dx = xe[1] - xe[0]
            self.data_cos = np.array(data_cos)
            self.templates_cos = np.array(templates_cos)
            self.lumi = lumi
            self.errordef = Minuit.LIKELIHOOD
            self.coeffs = coeffs
        else:
            raise Exception("Not implemented yet")

    def __call__(self, xi):
        xi = np.array(xi) / self.coeffs
        lam_cos = np.sum(self.templates_cos * xi * self.lumi * self.dx, axis=1)
        for i in range(len(lam_cos)):
            if lam_cos[i] < 0:
                print(lam_cos[i])
        l = -1 * np.sum((np.log(lam_cos).T * self.data_cos).T - lam_cos)
        return l


class TemplateFits():

    def __init__(self, xe_cos, xe_phi, templates_cos, data_cos, templates_phi, data_phi, coeffs=np.ones(9), lumi=140,
                 _2d=False):
        if not _2d:
            self.xe_cos = np.array(xe_cos)
            self.xe_phi = np.array(xe_phi)
            self.dx_cos = xe_cos[1] - xe_cos[0]
            self.dx_phi = xe_phi[1] - xe_phi[0]
            self.data_cos = np.array(data_cos)
            self.templates_cos = np.array(templates_cos)
            self.data_phi = np.array(data_phi)
            self.templates_phi = np.array(templates_phi)
            self.lumi = lumi
            self.errordef = Minuit.LIKELIHOOD
            self.coeffs = coeffs
        else:
            raise Exception("Not implemented yet")

    def __call__(self, xi):
        xi = np.array(xi) / self.coeffs
        lam_cos = np.sum(self.templates_cos * xi * self.lumi * self.dx_cos, axis=1)
        lam_phi = np.sum(self.templates_phi * xi * self.lumi * self.dx_phi, axis=1)
        for i in range(len(lam_cos)):
            if lam_cos[i] < 0:
                print(lam_cos[i])
        l = -1 * (np.sum((np.log(lam_cos).T * self.data_cos).T - lam_cos) + np.sum(
            (np.log(lam_phi).T * self.data_phi).T - lam_phi))
        return l


class TemplateFits_2d():

    def __init__(self, dx, templates, data, coeffs=np.ones(9), lumi=140, _2d=False):
        if not _2d:
            self.dx = dx
            self.data = np.array(data)
            self.templates = np.array(templates)
            self.lumi = lumi
            self.errordef = Minuit.LIKELIHOOD
            self.coeffs = coeffs

        else:
            raise Exception("Not implemented yet")

    def __call__(self, xi):
        xi = np.array(xi) / self.coeffs
        lam = xi[-1] * np.sum(self.templates * np.append(xi[:-1], 1) * self.lumi * self.dx, axis=1)
        for i in range(len(lam)):
            null = []
            if lam[i] == 0:
                null.append(i)
                # print(i)
            # if lam[i] < 0:
            #     print(f"lam[{i}]: ", lam[i])
        if null != []:
            pass
            # print("lam[i] = 0: ", null)
        # try:
        l = -1 * np.sum((np.log(lam).T * self.data).T - lam)
        # except RuntimeWarning as e:
        #     print(e)
        #     print("runtime warning")
        #     lam_new = []
        #     data_new = []
        #     for j in range(len(lam)):
        #         if not j in null:
        #             lam_new.append(lam[j])
        #             data_new.append(self.data[j])
        #     lam_new = np.array(lam_new)
        #     data_new = np.array(data_new)
        #
        #     l = -1 * np.sum((np.log(lam_new).T * data_new).T - lam_new)
        return l


def rebin_histogram(counts, n):
    if len(counts) % n != 0:
        counts = counts[:len(counts) - (len(counts) % n)]  # Trim to make it divisible by n
    return counts.reshape(-1, n).sum(axis=1)


def rebin_2d_histogram(hist2d, m, n, edges=None):
    # Check if the histogram dimensions are divisible by m and n; if not, trim them
    new_shape = (hist2d.shape[0] // m * m, hist2d.shape[1] // n * n)
    hist2d_trimmed = hist2d[:new_shape[0], :new_shape[1]]

    # Reshape and sum to rebin
    rebinned_hist = hist2d_trimmed.reshape(new_shape[0] // m, m, new_shape[1] // n, n).sum(axis=(1, 3))

    # if edges rebin edges
    if edges:
        rebinned_x_edges = edges[0, :][::n]
        rebinned_y_edges = edges[1, :][::n]
        return rebinned_hist, np.array((rebinned_x_edges, rebinned_y_edges))

    return rebinned_hist


def rebin_helper(hist, m, n):
    if hist.ndim == 1:
        newdim = int(np.sqrt(len(hist)))
        temp = hist.reshape(newdim, newdim)
        new = rebin_2d_histogram(temp, m, n)
        new = new.flatten()
        return new

    elif hist.ndim == 2:
        t = []
        for i in range(hist.shape[0]):
            newdim = int(np.sqrt(len(hist[i, :])))
            temp = hist[i, :].reshape(newdim, newdim)
            new = rebin_2d_histogram(temp, m, n)
            new = new.flatten()
            t.append(new)
        t = np.array(t)
        return t


def rebin(data, temps, m, n):
    return rebin_helper(data, m, n), rebin_helper(temps, m, n)