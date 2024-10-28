#!/usr/bin/env python3

'''
This file contains the style settings for the paper plots
'''

import matplotlib.pyplot as plt
import numpy as np

marker_sizes ={
    'default' : 6,
    'truth' : 0
} 

linestyles = {
    'default' : "",
    'truth' : "--"
}

lwidth = 1.0
alph = 1.0

colors = {
    "truth": "black",
    "meas" : "red",
    "unfold": "blue",
    }

markers = {
    'diff' : 'o',
    'ref' : 's',
    'truth' : '',
}


superscripts = {
    'truth' : 'MC',
    'meas' : 'Uncorr',
    'unfold' : 'Corr',
    'pythia' : 'Truth',
}


base_sty = dict(alpha = alph)

diff_marker_sty = dict(marker= markers['diff'] , markersize = marker_sizes['default'] , linestyle = linestyles['default'], linewidth = lwidth)
ref_marker_sty = dict(marker = markers['ref'] , markersize = marker_sizes['default'] , linestyle = linestyles['default'], linewidth = lwidth)
truth_marker_sty = dict(marker = markers['truth'] , markersize = marker_sizes['truth'] , markerfacecolor='none', markeredgecolor='none' , linestyle = linestyles['truth'], linewidth = 2.0)
two_part_marker_sty = dict(markerfacecolor='white', markeredgewidth=1.0)

meas_color_sty = dict(color=colors["meas"])
truth_color_sty = dict(color=colors["truth"])
unfold_color_sty = dict(color=colors["unfold"])

diff_pythia_two_sty = {**base_sty, **diff_marker_sty, **truth_color_sty, **two_part_marker_sty}
diff_pythia_four_sty = {**base_sty, **diff_marker_sty, **truth_color_sty}
diff_meas_two_sty = {**base_sty, **diff_marker_sty, **meas_color_sty, **two_part_marker_sty}
diff_meas_four_sty = {**base_sty, **diff_marker_sty, **meas_color_sty}
diff_unfold_two_sty = {**base_sty, **diff_marker_sty, **unfold_color_sty, **two_part_marker_sty}
diff_unfold_four_sty = {**base_sty, **diff_marker_sty, **unfold_color_sty}
diff_truth_sty = {**base_sty, **truth_marker_sty, **truth_color_sty}

ref_truth_sty = {**base_sty, **truth_marker_sty, **truth_color_sty}
ref_meas_two = {**base_sty, **ref_marker_sty, **meas_color_sty, **two_part_marker_sty}
ref_meas_four = {**base_sty, **ref_marker_sty, **meas_color_sty}


panel_labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)', '(q)', '(r)', '(s)', '(t)', '(u)', '(v)', '(w)', '(x)', '(y)', '(z)']



def GetLabel(flow_type, n , sup = None , k = None):
    if flow_type not in ['d', 'r']:
        print("Invalid flow type:", flow_type)
        return None
    label=r"$\mathcal{v}"
    if flow_type == 'd':
        label+=r"\mathbf{'}"
    label+=r"_{\mathcal{" +f"{n}" + r"}"
    if k is not None:
        label+=r"\mathcal{\{" + f"{k}" + r" \}}"
    label+=r"}"
    if sup is not None:
        label+=r"^{\mathrm{" + f"{sup}" + r"}}"
    label+=r"$"
    return label

def GetYaxisLabel(flow_type, n):
    if flow_type not in ['d', 'r']:
        print("Invalid flow type:", flow_type)
        return None

    label=r"$\mathcal{v}"
    if flow_type == 'd':
        label+=r"\mathbf{'}^{\mathrm{Jet}}"
    label+=r"_{\mathcal{" +f"{n}" + r"}}$"
    
    return label

def GetpTaxisLabel():
    return r"$\mathcal{p}_{T,\mathrm{Jet}}$ [GeV/c]"

myStyle = {
    'axes.facecolor':'white',
    'axes.edgecolor': 'black',
    'figure.facecolor' : 'white',
    'figure.edgecolor' : 'white',
    'figure.frameon' : False,
    'figure.subplot.left' : 0.16,
    'figure.subplot.right' : 0.95,
    'figure.subplot.bottom' : 0.16,
    'figure.subplot.top' : 0.95,

    'axes.labelpad' : 1.4,
    'axes.labelcolor' : 'black',
    'axes.labelsize': 12,
    'axes.linewidth' : 1.0,
    'axes.titlesize' : 'medium',
    'xaxis.labellocation': 'center',
    'yaxis.labellocation': 'top',
    'axes.titlepad': 0.0,
    
    # Helvetica
    'text.color': 'black',
    'font.family': 'Times New Roman',
    'mathtext.fontset': 'stix',
    'font.size': 12,
    'legend.fontsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'text.usetex': False,
    'mathtext.it' : 'sans:italic',
    'mathtext.default' : 'regular',

    # ticks
    'xtick.major.size' : 7,
    'xtick.minor.size' : 5,
    'xtick.major.width' : 1,
    'xtick.minor.width' : 1,
    'xtick.direction' : 'in',
    'ytick.major.size' : 7,
    'ytick.minor.size' : 5,
    'ytick.major.width' : 1,
    'ytick.minor.width' : 1,
    'ytick.direction' : 'in',
    'ytick.right' : True,
    'xtick.top' : True,
    'xtick.minor.visible' : True,
    'ytick.minor.visible' : True,
    'xtick.major.pad': 2,
    'ytick.major.pad': 2,

    'lines.linewidth': 1.0,
    'lines.linestyle': '-',
    'lines.markersize' : 5,
    'lines.marker' : 'o',
    'lines.markeredgewidth' : 1,

    'savefig.transparent': False,
    'errorbar.capsize': 2,

    # legend
    'legend.frameon': False,
    'legend.loc': 'best',
    'legend.numpoints': 1,
    'legend.scatterpoints': 1,
    'legend.fontsize': 12,
    'legend.handlelength': 1.5,
    'legend.handletextpad': 0.5,
    'legend.labelspacing': 0.5,
    'legend.borderpad': 0.4,
    'legend.borderaxespad': 0.5,
    'legend.columnspacing': 0.7,
    'legend.markerscale': 1.0,
    'legend.shadow': False
}

def GetHistosInFile(filename):
    import ROOT
    import numpy as np
    
    f = ROOT.TFile(filename)
    histo_1dnames = [key.GetName() for key in f.GetListOfKeys() if key.GetClassName() == 'TH1D' or key.GetClassName() == 'TH1F']
    histo_2dnames = [key.GetName() for key in f.GetListOfKeys() if key.GetClassName() == 'TH2D' or key.GetClassName() == 'TH2F']
    T1D_hists = {}
    T2D_hists = {}
    for name in histo_1dnames:
        th1d = f.Get(name)
        x = np.array([th1d.GetBinCenter(i) for i in range(1,th1d.GetNbinsX()+1)])
        x_err = np.array([th1d.GetBinWidth(i)/2 for i in range(1,th1d.GetNbinsX()+1)])
        y = np.array([th1d.GetBinContent(i) for i in range(1,th1d.GetNbinsX()+1)])
        y_err = np.array([th1d.GetBinError(i) for i in range(1,th1d.GetNbinsX()+1)])
        X = [x,x_err]
        Y = [y,y_err]
        T1D_hists[name] = {'x': X[0], 'x_err': X[1], 'y': Y[0], 'y_err': Y[1]}
    for name in histo_2dnames:
            th2d = f.Get(name)
            xbins = np.array([th2d.GetXaxis().GetBinLowEdge(i) for i in range(1,th2d.GetNbinsX()+1)])
            ybins = np.array([th2d.GetYaxis().GetBinLowEdge(i) for i in range(1,th2d.GetNbinsY()+1)])
            values = []
            value_errors = []
            for i in range(1,th2d.GetNbinsX()+1):
                tmpvalues = []
                tmpvalue_errors = []
                for j in range(1,th2d.GetNbinsY()+1):
                    tmpvalues.append(th2d.GetBinContent(i,j))
                    tmpvalue_errors.append(th2d.GetBinError(i,j))
                values.append(tmpvalues)
                value_errors.append(tmpvalue_errors)
            values = np.array(values)
            value_errors = np.array(value_errors)  
            # values = np.array([th2d.GetBinContent(i,j) for i in range(1,th2d.GetNbinsX()+1) for j in range(1,th2d.GetNbinsY()+1)])
            # value_errors = np.array([th2d.GetBinError(i,j) for i in range(1,th2d.GetNbinsX()+1) for j in range(1,th2d.GetNbinsY()+1)])
            T2D_hists[name] = {'xbins': xbins, 'ybins': ybins, 'values': values, 'value_errors': value_errors}
     
    f.Close()
    return T1D_hists, T2D_hists

def GetFigSize(fig_width_pt=510.0, denominator=3.0):
    inches_per_pt = 1.0/72.27
    golden_mean = (np.sqrt(5)-1.0)/denominator
    fig_width = fig_width_pt*inches_per_pt
    fig_height = fig_width*golden_mean
    return [fig_width, fig_height]

def ReadFlowResults(filename):
    th1d, _ = GetHistosInFile(filename)
    ref_flow = ReadTTree(filename, 'tree')
    return th1d, ref_flow
    
def ReadTTree(filename, treename, vars=None):
    import uproot
    import numpy as np 
    import awkward as ak

    f = uproot.open(filename)
    tree = f[treename]
    data = ak.to_list(tree.arrays(vars))[0]
    return data
      
def GetUnfoldedIndices(th1d):
    names = [f'{n}' for n in th1d.keys() if 'unfolded' in n]
    indices = [int(n.split('_')[-1]) for n in names]
    indices = np.unique(indices).astype(int)
    return indices

def GetPlot(data, xmin = None, xmax = None, ymin = None, ymax = None):
    x = data['x']
    y = data['y']
    xerr = data['x_err']
    yerr = data['y_err']
    if xmin is not None:
        cut = np.where(x > xmin)
        x = x[cut]
        y = y[cut]
        xerr=xerr[cut]
        yerr=yerr[cut]

    if xmax is not None:
        cut = np.where(x < xmax)
        x = x[cut]
        y = y[cut]
        xerr=xerr[cut]
        yerr=yerr[cut]

    if ymin is not None:
        cut = np.where(y > ymin)
        x = x[cut]
        y = y[cut]
        xerr=xerr[cut]
        yerr=yerr[cut]

    if ymax is not None:
        cut = np.where(y < ymax)
        x = x[cut]
        y = y[cut]
        xerr=xerr[cut]
        yerr=yerr[cut]

    return x, y ,xerr, yerr
