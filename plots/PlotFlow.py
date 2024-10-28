#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os 
import sys
import datetime
import argparse

import paperstyle as mystl

plt.rcParams.update(mystl.myStyle)
plt.rc('axes', unicode_minus=False)

def PlotRefFow(fin, fout=None):
    # print("Plotting Reference Flow")
    _, ref_data = mystl.ReadFlowResults(fin)
    fig_size = mystl.GetFigSize(510, 4.5)
    fig , axs = plt.subplots(1,3, figsize=fig_size, dpi=200, constrained_layout=True)
    
    for iharm in range(3):
        ax = axs[iharm]
        
        two_name = f'v{iharm+2}_two_ref'
        four_name = f'v{iharm+2}_four_ref'
        truth_name = f'v{iharm+2}_ref_truth'
        two_label = mystl.GetLabel('r',iharm +2 , k = 2)
        four_label = mystl.GetLabel('r',iharm +2 , k = 4)
        t_label = mystl.GetLabel('r', iharm+ 2, sup=mystl.superscripts['truth'])


        y_m2, y_m2_err = ref_data[two_name], ref_data[f'{two_name}_err'] 
        y_m4 , y_m4_err = ref_data[four_name], ref_data[f'{four_name}_err']
        y_t = ref_data[truth_name]
        ax.errorbar(1.5, y_m4, yerr=y_m4_err, **mystl.ref_meas_four, label = four_label)
        ax.errorbar(0.5, y_m2, yerr=y_m2_err, **mystl.ref_meas_two, label = two_label)
        
        ax.plot([0,2], [y_t, y_t], **mystl.ref_truth_sty, label = t_label) 
        
        # turn off minor ticks
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.tick_params(axis='x', which='both', top=False)
        ax.set_ylim(ax.get_ylim()[0]*0.99, ax.get_ylim()[1]*1.01)
        ax.set_xlim(0,2)
        x_ticklabels = [ two_label, four_label]
        x_ticks = [0.5,1.5]
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize=10)
        yticks = [f'{y:.4f}' for y in ax.get_yticks()]
        ax.set_yticks(ax.get_yticks())
        ax.set_yticklabels(yticks, fontsize=8)

        px = 0.85
        py = 0.92
        ax.text(px,py,mystl.panel_labels[iharm], transform=ax.transAxes, weight='bold', fontsize=8)
        if iharm == 0:
            tx = 0.06
            ty = 0.92
            ax.text(tx,ty,'Reference Flow', transform=ax.transAxes, weight='bold', fontsize=8)
            ax.text(tx,ty-0.08, "PYTHIA+TennGen", transform=ax.transAxes, fontsize=8)
            ax.text(tx,ty-0.16,r'$Au+Au$ $\sqrt{s_{NN}}$ = 200 GeV', transform=ax.transAxes, fontsize=8)
            ax.text(tx,ty-0.24,r'20-30% Central', transform=ax.transAxes, fontsize=8)

    if fout is not None:
        plt.savefig(fout)
        plt.close()
    else:
        plt.show()
    return

def PlotDiffFlow(fin, uiter, fout = None , cuts = {'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None}):
    # print("Plotting Differential Flow")
    th1ds, _ = mystl.ReadFlowResults(fin)

    fig_size = mystl.GetFigSize(510, 4.5)
    fig , axs = plt.subplots(1,3, figsize=fig_size, dpi=200, constrained_layout=True)

    xmin = cuts['xmin']
    xmax = cuts['xmax']
    ymin = cuts['ymin']
    ymax = cuts['ymax']


    for iharm in range(3):

        ax = axs[iharm]

        two_meas = f'h1_v{iharm+2}_two_diff_meas'
        four_meas = f'h1_v{iharm+2}_four_diff_meas'
        two_unfolded =  f'h1_v{iharm+2}_two_diff_unfolded_{uiter}'
        four_unfolded =  f'h1_v{iharm+2}_four_diff_unfolded_{uiter}'
        two_pythia = f'h1_v{iharm+2}_two_diff_truth'
        four_pythia = f'h1_v{iharm+2}_four_diff_truth'
        truth_name = f'h1_jet_v{iharm+2}_func'

        x_2m, y_2m, _, y_2m_err = mystl.GetPlot(th1ds[two_meas], xmin, xmax, ymin, ymax)
        x_4m, y_4m, _, y_4m_err = mystl.GetPlot(th1ds[four_meas],  xmin, xmax, ymin, ymax)
        x_2u, y_2u, _, y_2u_err = mystl.GetPlot(th1ds[two_unfolded] , xmin, xmax, ymin, ymax)
        x_4u, y_4u, _, y_4u_err = mystl.GetPlot(th1ds[four_unfolded] , xmin, xmax, ymin, ymax)
        x_2p, y_2p, _, y_2p_err = mystl.GetPlot(th1ds[two_pythia], xmin, xmax, ymin, ymax)
        x_4p, y_4p, _, y_4p_err = mystl.GetPlot(th1ds[four_pythia] , xmin, xmax, ymin, ymax)
        x_t, y_t, _ , _ = mystl.GetPlot(th1ds[truth_name], 0,50)


        ax.errorbar(x_2m, y_2m, yerr=y_2m_err, **mystl.diff_meas_two_sty, label = mystl.GetLabel("d", iharm+2, k=2, sup = mystl.superscripts['meas']))
        ax.errorbar(x_4m, y_4m, yerr=y_4m_err, **mystl.diff_meas_four_sty, label =  mystl.GetLabel("d", iharm+2, k=4, sup = mystl.superscripts['meas']))
        ax.errorbar(x_2u, y_2u, yerr=y_2u_err, **mystl.diff_unfold_two_sty, label =  mystl.GetLabel("d", iharm+2, k=2, sup = mystl.superscripts['unfold']))
        ax.errorbar(x_4u, y_4u, yerr=y_4u_err, **mystl.diff_unfold_four_sty, label = mystl.GetLabel("d", iharm+2, k=4, sup = mystl.superscripts['unfold']))
        ax.errorbar(x_2p, y_2p, yerr=y_2p_err, **mystl.diff_pythia_two_sty, label = mystl.GetLabel("d", iharm+2, k=2, sup = mystl.superscripts['pythia']))
        ax.errorbar(x_4p, y_4p, yerr=y_4p_err, **mystl.diff_pythia_four_sty, label = mystl.GetLabel("d", iharm+2, k=4, sup = mystl.superscripts['pythia']))
        ax.plot(x_t, y_t, **mystl.diff_truth_sty, label = mystl.GetLabel("d", iharm+2, sup = mystl.superscripts['truth']))

        ax.set_xlabel(mystl.GetpTaxisLabel(), loc='center')
        px = 0.85
        py = 0.92
        ax.text(px,py,mystl.panel_labels[iharm], transform=ax.transAxes, weight='bold')
        if iharm == 0:
            ax.set_ylabel(mystl.GetYaxisLabel('d', iharm+2), loc='center')
            tx = 0.06
            ty = 0.92
            # axs.set_xlim(0, 0.2)
            ax.text(tx,ty,'Differential Flow', transform=ax.transAxes, weight='bold', fontsize=8)
            ax.text(tx,ty-0.08, "PYTHIA+TennGen", transform=ax.transAxes, fontsize=8)
            ax.text(tx,ty-0.16,r'$Au+Au$ $\sqrt{s_{NN}}$ = 200 GeV', transform=ax.transAxes, fontsize=8)
            ax.text(tx,ty-0.24,r'20-30% Central', transform=ax.transAxes, fontsize=8)
            ax.legend(loc='lower left', fontsize=8, frameon=False, ncol=4)

    if fout is not None:
        plt.savefig(fout)
        plt.close()

    else:
        plt.show()
    return

def PlotDiffMeasFlow(fin, iharm, fout = None,  cuts = {'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None}):

    # print("Plotting Differential Measured Flow")
    th1ds, _ = mystl.ReadFlowResults(fin)

    fig_size = mystl.GetFigSize(247, 3.5)
    fig = plt.figure(figsize=fig_size, dpi=200, constrained_layout=True)
    ax = fig.add_subplot(111)

    two_meas = f'h1_v{iharm}_two_diff_meas'
    four_meas = f'h1_v{iharm}_four_diff_meas'
    two_pythia = f'h1_v{iharm}_two_diff_truth'
    four_pythia = f'h1_v{iharm}_four_diff_truth'
    truth_name = f'h1_jet_v{iharm}_func'
    
    
    xmin = cuts['xmin']
    xmax = cuts['xmax']
    ymin = cuts['ymin']
    ymax = cuts['ymax']

    x_2m, y_2m, _, y_2m_err = mystl.GetPlot(th1ds[two_meas], xmin, xmax, ymin, ymax)
    x_4m, y_4m, _, y_4m_err = mystl.GetPlot(th1ds[four_meas], xmin, xmax, ymin, ymax)
    x_2p, y_2p, _, y_2p_err = mystl.GetPlot(th1ds[two_pythia], xmin, xmax, ymin, ymax)
    x_4p, y_4p, _, y_4p_err = mystl.GetPlot(th1ds[four_pythia], xmin, xmax, ymin, ymax)
    x_t, y_t, _ , _ = mystl.GetPlot(th1ds[truth_name], xmin, xmax, ymin, ymax)


    ax.errorbar(x_2m, y_2m, yerr=y_2m_err, **mystl.diff_meas_two_sty, label = mystl.GetLabel("d", iharm, k=2, sup = mystl.superscripts['meas']))
    ax.errorbar(x_4m, y_4m, yerr=y_4m_err, **mystl.diff_meas_four_sty, label = mystl.GetLabel("d", iharm, k=4, sup = mystl.superscripts['meas']))
    ax.errorbar(x_2p, y_2p, yerr=y_2p_err, **mystl.diff_pythia_two_sty, label = mystl.GetLabel("d", iharm, k=2, sup = mystl.superscripts['pythia']))
    ax.errorbar(x_4p, y_4p, yerr=y_4p_err, **mystl.diff_pythia_four_sty, label = mystl.GetLabel("d", iharm, k=4, sup = mystl.superscripts['pythia']))
    ax.plot(x_t, y_t, **mystl.diff_truth_sty, label = mystl.GetLabel("d", iharm, sup = mystl.superscripts['truth']))
    

    px = 0.05
    py = 0.05
    ax.text(px,py,f"Measured $v_{iharm}$", transform=ax.transAxes, weight='bold', fontsize=8)

    ax.legend(loc='lower left', fontsize=8, frameon=False, ncol=2)

    ax.set_ylabel(mystl.GetYaxisLabel('d', iharm), loc='center')
    ax.set_xlabel(mystl.GetpTaxisLabel(), loc='center')

    if fout is not None:
        plt.savefig(fout)
        plt.close()
    else:
        plt.show()

    return

def PlotDiffUnfoldFlow(fin, iharm, fout = None,  cuts = {'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None}):

    # print("Plotting Differential Unfolded Flow")
    th1ds, _ = mystl.ReadFlowResults(fin)

    fig_size = mystl.GetFigSize(247, 3.5)
    fig = plt.figure(figsize=fig_size, dpi=200, constrained_layout=True)
    ax = fig.add_subplot(111)

    uindices = mystl.GetUnfoldedIndices(th1ds)    
    truth_name = f'h1_jet_v{iharm}_func'
    
    
    xmin = cuts['xmin']
    xmax = cuts['xmax']
    ymin = cuts['ymin']
    ymax = cuts['ymax']
    
    for u in uindices:
        u_two = f'h1_v{iharm}_two_diff_unfolded_{u}'
        x_2m, y_2m, _, y_2m_err = mystl.GetPlot(th1ds[u_two], xmin, xmax, ymin, ymax)
        u_four = f'h1_v{iharm}_four_diff_unfolded_{u}'
        x_4m, y_4m, _, y_4m_err = mystl.GetPlot(th1ds[u_four], xmin, xmax, ymin, ymax)
        ax.errorbar(x_2m, y_2m, yerr=y_2m_err, label = mystl.GetLabel("d", iharm, k=2, sup = f"iter {u}"))
        ax.errorbar(x_4m, y_4m, yerr=y_4m_err, label = mystl.GetLabel("d", iharm, k=4, sup = f"iter {u}"))


    x_t, y_t, _ , _ = mystl.GetPlot(th1ds[truth_name], xmin, xmax, ymin, ymax)
    ax.plot(x_t, y_t, **mystl.diff_truth_sty, label = mystl.GetLabel("d", iharm, sup = mystl.superscripts['truth']))


    ax.set_ylabel(mystl.GetYaxisLabel('d', iharm), loc='center')
    ax.set_xlabel(mystl.GetpTaxisLabel(), loc='center')
    px = 0.05
    py = 0.05
    ax.text(px,py,f"Unfolded $v_{iharm}$", transform=ax.transAxes, weight='bold', fontsize=8)

    ax.legend(loc='lower left', fontsize=8, frameon=False, ncol=4)

    

    if fout is not None:
        plt.savefig(fout)
        plt.close()
    
    else:
        plt.show()
    return

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Plot Flow Results')
    parser.add_argument('-o', '--output', type=str, help='Output directory', default=None)
    parser.add_argument('-f', '--files', nargs='+', help='Input files', default=None)
    parser.add_argument('-r', '--rotations', nargs='+', help='Rotations', default=None)
    args = parser.parse_args()

    OUTPUT_PATH='/lustre/isaac/scratch/tmengel/jet-vn-cumulant/plots/'
    output_dir=OUTPUT_PATH
    if args.output is None:
        now = datetime.datetime.now().strftime("%h_%d_%Y_%H%M")
        output_dir = f'{OUTPUT_PATH}/{now}'
    else:
        output_dir = args.output
    
    os.makedirs(output_dir, exist_ok=True)

    
    ROTATIONS = ['ConstantA', 'ConstantB','ConstantC', 'Cosine','Linear','Logarithmic']
    rotations = []
    if args.rotations is not None:
        rotations = args.rotations
    else:
        rotations = ROTATIONS

    FILES_DIR ='/lustre/isaac/scratch/tmengel/jet-vn-cumulant/ana/rootfiles/'
    RES_FILE='flow_results.root'
    files = []
    if args.files is not None:
        files = args.files
    else:
        files = [f'{FILES_DIR}/{rot}/{RES_FILE}' for rot in rotations]


    print (f'Output directory: {output_dir}')
    print (f'Files: {files}')
    print (f'Rotations: {rotations}')
    for f, rot in zip(files, rotations):
        if not os.path.exists(f):
            print(f'File {f} does not exist. Skipping...')
            continue

        PlotRefFow(f, f'{output_dir}/{rot}_ref_flow.png')

        no_cuts = {'xmin': 0, 'xmax': 60, 'ymin': None, 'ymax': None}

        th1ds, _ = mystl.ReadFlowResults(f)
        uindices = mystl.GetUnfoldedIndices(th1ds)
        for uiter in uindices:
            PlotDiffFlow(f, uiter, f'{output_dir}/{rot}_diff_flow_iter{uiter}.png', no_cuts)
        
        for iharm in range(2,5):
            PlotDiffMeasFlow(f, iharm, f'{output_dir}/{rot}_diff_meas_v{iharm}.png', no_cuts)
            PlotDiffUnfoldFlow(f, iharm, f'{output_dir}/{rot}_diff_unfold_v{iharm}.png', no_cuts)
        
        ymin = 0
        ymax = 0
        for iharm in range(2,5):
            jet_func = f'h1_jet_v{iharm}_func'
            cuts = {'xmin': 0,
                    'xmax': 60,
                    'ymin': 0.5*np.min(th1ds[jet_func]['y']),
                    'ymax': 1.5*np.max(th1ds[jet_func]['y'])}
            
            if cuts['ymin'] < ymin:
                ymin = cuts['ymin']
            if cuts['ymax'] > ymax:
                ymax = cuts['ymax']
            PlotDiffMeasFlow(f, iharm, f'{output_dir}/{rot}_diff_meas_v{iharm}_zoomed.png', cuts)
            PlotDiffUnfoldFlow(f, iharm, f'{output_dir}/{rot}_diff_unfold_v{iharm}_zoomed.png', cuts)

        all_cuts = {'xmin': 0, 'xmax': 60, 'ymin': ymin, 'ymax': ymax}
        for uiter in uindices:
            PlotDiffFlow(f, uiter, f'{output_dir}/{rot}_diff_flow_iter{uiter}_zoomed.png', all_cuts)

        # Copy files to output directory
        str_cmd = f'cp {f} {output_dir}/{rot}_flow_results.root'
        os.system(str_cmd)

        print(f'Copied {f} to {output_dir}/{rot}_flow_results.root')

    print('Done')
    sys.exit(0)