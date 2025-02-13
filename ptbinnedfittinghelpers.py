import numpy as np
import json

# loading data
def loadbins(path='/home/maltem/Nextcloud/TU_Dresden/Bachelorthesis/minuitfitting/data_ptbinned/histo_nocuts'):
    with open(path+'/datbin0_2D.json', 'r') as f:
        bins = json.load(f)['bins']
    return bins

def loaddata(path='/home/maltem/Nextcloud/TU_Dresden/Bachelorthesis/minuitfitting/data_ptbinned/histo_nocuts', numberbins=36, _load='2D', overflow=False):
    print(f'loading {_load} data from {path} with {numberbins} bins')
    if _load == 'full':
        data = {
            'xsec': [],
            '2D': [],
            'cos': [],
            'phi': []
        }
    elif _load == '2D':
        data = {
            'xsec': [],
            '2D': []
        }
    elif _load == 'cos':
        data = {
            'xsec': [],
            'cos': []
        }
    elif _load == 'phi':
        data = {
            'xsec': [],
            'phi': []
        }
    elif _load == '1D':
        data = {
            'xsec': [],
            'cos': [],
            'phi': []
        }
    if overflow:
        data_of = {}
    for key in data.keys():
        for j in range(numberbins):

            with open(path+'/datbin'+str(j)+'_'+key+'.json', 'r') as f:
                dat = json.load(f)
                for n in dat.keys():
                    if n.find("vals") != -1:
                        # print(n)
                        k = n
                        break
                else: print("no key found", dat.keys()); raise Exception()
                print(key, k)
                data[key].append(np.array(dat[k]))
        if overflow:

            with open(path+'/overflow_'+key+'.json', 'r') as f:
                dat = json.load(f)
                for n in dat.keys():
                    if n.find("vals") != -1:
                        print(n)
                        k = n
                        break
                else: print("no key found", dat.keys()); raise Exception()
                print(key, k)
                data_of[key] = np.array(dat[k])
    if overflow:
        return data, data_of
    return data

def countlines(file):
    with open(file, 'r') as f:
        s = sum(1 for line in f)
    return s

def loadptbins(file):
    ptbins = []
    with open(file, 'r') as f:
        for line in f:
            ptbins.append(float(line))
    return ptbins

def loadtemplates(path, numberbins=36, _load='2D', overflow=False):
    print(f'loading {_load} templates from {path} with {numberbins} bins')

    if _load == 'full':
        data = {
            '2D': [[] for i in range(numberbins)],
            'cos': [[] for i in range(numberbins)],
            'phi': [[] for i in range(numberbins)]
        }
    elif _load == '2D':
        data = {
            '2D': [[] for i in range(numberbins)]
        }
    elif _load == 'cos':
        data = {
            'cos': [[] for i in range(numberbins)]
        }
    elif _load == 'phi':
        data = {
            'phi': [[] for i in range(numberbins)]
        }
    elif _load == '1D':
        data = {
            'cos': [[] for i in range(numberbins)],
            'phi': [[] for i in range(numberbins)]
        }

    else:
        print(f'_load argument {_load} not recognized\n'
              f'Please choose from full, 2D, cos, phi, 1D')
        raise Exception("_load arg not recognized")
    if overflow:
        data_of = {}
        for key in data.keys():
            data_of[key] = [[] for i in range(9)]
    for key in data.keys():
        for j in range(9):
            for i in range(numberbins):
                print("LOADING: ", path+'/ptbin'+str(i)+'_'+key+str(j)+'.json')
                with open(path+'/ptbin'+str(i)+'_'+key+str(j)+'.json', 'r') as f:
                    dat = json.load(f)
                    for n in dat.keys():
                        if n.find("vals") != -1:
                            # print(n)
                            k = n
                            break
                    else:
                        print("no key found", dat.keys())
                        raise Exception()
                    print(key, k)
                    data[key][i].append(np.array(dat[k]))
            if overflow:
                print("LOADING: ", path+'/ptoverflow_' + key +str(j)+'.json')

                with open(path+'/ptoverflow_' + key +str(j)+'.json', 'r') as f:
                    dat = json.load(f)
                    for n in dat.keys():
                        if n.find("vals") != -1:
                            print(n)
                            k = n
                            break
                    else: print("no key found", dat.keys()); raise Exception()
                    data_of[key][j] = np.array(dat[k])
    if overflow:
        return data, data_of
    return data

def loadcoeffs(pathtocoeffs, nptbins, overflow=False):
    #TODO
    coeffs = np.zeros((nptbins, 9))
    with open(pathtocoeffs+'coeffs.txt', 'r') as f:
        for i in range(9):
            for j in range(nptbins):
                coeffs[j][i] = float(f.readline())
    if overflow:
        coeffs_of = np.zeros(9)
        with open(pathtocoeffs+'coeffs_of.txt', 'r') as f:
            for i in range(9):
                coeffs_of[i] = float(f.readline())
        return coeffs, coeffs_of
    return coeffs

def load(pathtodata, pathtotemplates, pathtoptbins, pathtocoeffs, _load='2D', overflow=False):
    print(_load)
    ptbins = loadptbins(pathtoptbins)
    nptbins = len(ptbins) - 1
    ang_bins = loadbins(pathtodata)
    if overflow:
        data, data_of = loaddata(pathtodata, nptbins, _load, overflow)
        temps, temps_of = loadtemplates(pathtotemplates, nptbins, _load, overflow)
        coeffs, coeffs_of = loadcoeffs(pathtocoeffs, nptbins, overflow)
        return data, data_of, temps, temps_of, ptbins, ang_bins, coeffs, coeffs_of
    else:
        data = loaddata(pathtodata, nptbins, _load, overflow)
        temps = loadtemplates(pathtotemplates, nptbins, _load, overflow)
        coeffs = loadcoeffs(pathtocoeffs, nptbins, overflow)
        return data, temps, ptbins, ang_bins, coeffs

def is_iterable(obj):
    try:
        e_ = (e for e in obj)
        return True
    except TypeError:
        return False

