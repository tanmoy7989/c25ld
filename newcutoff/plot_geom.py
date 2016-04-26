import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator

# GLOBALS
analysis_dir = os.path.expanduser('./geom_prop/c25_longtraj')
Measures = ['Rg', 'R_EE', 'kappa', 'SASA_atom', 'pmf_Rg_R_EE']
cgtypes = ['AA', 'SP', 'SPLD', 'LD']
fftypes = ['lj', 'wca']
linestyles = {'AA': ('red', 'None', 'o'), 
              'SP': ('blue','-', 'None'), 
              'SPLD': ('green', '--', 'None' ), 
              'LD': ('black', ':', 'None')} # (color, linestyle, markerstyle)

axlbls = {'Rg': r'$Rg (\AA)$', 
          'R_EE': r'$R_{EE} (\AA)$',
          'kappa': r'$\kappa$', 
          'SASA_atom': r'$SASA (\AA^2/atom)$',
		  'pmf_Rg_R_EE': r'$W(Rg)$' + ' (kcal/mol)', 
	      'Dist': 'Distribution'}                        

plot_fmt = 'svg'
pad_left = 0.12
pad_bottom = 0.25
hspace = 0.2; wspace = 0.3

kB = 0.001987 #kcal/mol
Temp = 298.0
conv_factor = kB * Temp

# Extract data
chain_length = int(sys.argv[1])
savedir = os.path.abspath(sys.argv[2])

# Plot 1D histograms
for measure in Measures[:-1]:
    fig = plt.figure(figsize = (8,4), facecolor = 'w', edgecolor = 'w')
    axlist = {'lj': plt.subplot('121'), 'wca': plt.subplot('122')}
    for fftype in fftypes:
        ax = axlist[fftype]
        for cgtype in cgtypes:
            filename = os.path.join(savedir, 'hist', '%s_%s_%s.pickle' % (measure, fftype, cgtype))
            try:
                x,y = pickle.load(open(filename, 'r'))
            except IOError:
                print 'c%d (%s, %s) data not there for %s' % (chain_length, fftype, cgtype, measure)
                continue
            ax.plot(x, y, linewidth = 3, markersize = 8, label = cgtype, color = linestyles[cgtype][0], 
                    linestyle = linestyles[cgtype][1], marker = linestyles[cgtype][2])
            
        # design
        ax.set_xlabel(axlbls[measure], fontsize = 'large'); ax.set_ylabel(axlbls['Dist'], fontsize = 'large')           
        ax.set_title('c%d %s' % (chain_length, fftype), fontsize = 'large')
        ax.xaxis.set_major_locator(MaxNLocator(nbins = 6, prune = 'both'))
        ax.yaxis.set_major_locator(MaxNLocator(nbins = 6, prune = 'both'))
        
        if fftype == 'lj': 
            ax.legend(loc = 'best', prop = {'size': 'large'})
        plt.subplots_adjust(left = pad_left, bottom = pad_bottom, wspace = wspace, hspace = hspace)
        figname = os.path.join(savedir, '%s.%s' % (measure, plot_fmt))
        plt.savefig(figname, dpi = 300)


# Plot 2D histogram
del measure, ax
measure = 'pmf_Rg_R_EE'
for fftype in fftypes:
    fig = plt.figure(figsize = (8,8), facecolor = 'w', edgecolor = 'w')
    index = '22'
    for i, cgtype in enumerate(cgtypes):
        filename = os.path.join(savedir, 'hist', '%s_%s_%s.pickle' % (measure, fftype, cgtype))
        try:
            (x,y), pmf  = pickle.load(open(filename, 'r'))
        except IOError:
            print 'c%d (%s, %s) data not there for %s' % (chain_length, fftype, cgtype, measure)
            continue
        ax = plt.subplot(index+str(i+1))
        X = (x - x.min())/(x.max() - x.min())
        Y = (y - y.min())/(y.max() - y.min())
        Z = -np.log(pmf) * conv_factor
        im = ax.imshow(Z, cmap = cm.spectral_r, extent = [X.min(),X.max(),Y.min(),Y.max()])
        
        # design
        ax.annotate(cgtype, xy = (0.2, 0.2), fontsize = 'large', weight = 'bold')
        ax.xaxis.set_major_locator(MaxNLocator(nbins = 5, prune = 'both'))
        ax.yaxis.set_major_locator(MaxNLocator(nbins = 5, prune = 'both'))
        ax.set_xlabel(''); ax.set_ylabel('')
        if i == 0 or i == 1: 
            ax.set_xticklabels([])

        # colorbar
        cax = fig.add_axes([0.95, 0.1, 0.02, 0.8])
        cbar = fig.colorbar(im, cax = cax, orientation = 'vertical')
        cbar.set_label(r'$pmf(Rg, R_{EE}) (kcal/mol)$', size = 'large', rotation = 90)
    
    # central x and y labels 	
    plt.figtext(0.35, 0.03, 'Dimensionless ' + r'$Rg$', size = 'x-large')
    plt.figtext(0.03, 0.35, 'Dimensionless ' + r'$R_{EE}$', size = 'x-large', rotation = 90)

    plt.subplots_adjust(wspace = wspace, hspace = 0., left = pad_left, bottom = pad_bottom)
    figname = os.path.join(savedir, '%s_%s.%s' % (measure, fftype, plot_fmt))
    plt.savefig(figname, dpi = 300)

plt.show()
