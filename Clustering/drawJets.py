import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

eta_min = -3.
eta_max = 3.
bins_eta = 71
bins_phi = 72

extent = eta_min, eta_max, -np.pi, np.pi

eta_edges = np.linspace(eta_min, eta_max, bins_eta + 1)
phi_edges = np.linspace(-np.pi, np.pi, bins_phi + 1)

eta = np.linspace(eta_min, eta_max, bins_eta + 1)[:-1] + (eta_max - eta_min) / (2 * bins_eta)
phi = np.linspace(-np.pi, np.pi, bins_phi + 1)[:-1] + (np.pi / bins_phi)


def drawJets(jets, name, label=None, eta=None, phi=None, pt=None, genjet_pt=None,genjet_eta=None, genjet_phi=None, genjet_pid=None, genjet_matched=None, nJets = -1):
    n = 1
    p = 1

    fig = plt.figure(figsize=(12, 8))

    ax = None

    colors = cm.rainbow(np.linspace(0, 1, max(2,len(jets))))
    cmap = LinearSegmentedColormap.from_list('cmap', colors, len(colors))

    if n > 3:
        ax = fig.add_subplot(2, int((n+1)/2), p, sharey=ax)
    else:
        ax = fig.add_subplot(1, n, p, sharey=ax)
    area = np.zeros((eta_edges.shape[0] - 1, phi_edges.shape[0] - 1),
                    dtype=np.float64)

    if nJets>-1:
        jets = jets[:nJets]
    for ijet, jet in enumerate(jets):
        constit = jet.constituents_array()
        jetarea, _, _ = np.histogram2d(constit['eta'], constit['phi'],
                                       bins=(eta_edges, phi_edges))
        area += (jetarea > 0) * (ijet + 1)

    ax.imshow(np.ma.masked_where(area == 0, area).T, cmap=cmap,
              extent=extent, aspect=(eta_max - eta_min) / (2*np.pi),
              interpolation='none', origin='lower')

    if not eta is None:
        ax.scatter(eta, phi,
                   s=30* np.log10(pt)/np.log10(pt.max())  )
    
    if not genjet_eta is None:
        ax.scatter(genjet_eta[genjet_matched], genjet_phi[genjet_matched],
                   s=50, marker="x",color='black')
        ax.scatter(genjet_eta[np.logical_not(genjet_matched)], genjet_phi[np.logical_not(genjet_matched)],
                   s=50, marker="o",color='black')

        if genjet_pid is None:
            genjetList = "Gen Jets\n"
            for i in range(len(genjet_eta)):
                genjetList += "pt=%.2f, eta=%.3f, phi=%.3f\n"%(genjet_pt[i], genjet_eta[i],genjet_phi[i])
            if len(genjet_eta)==0:
                genjetList = "No Gen Jets"
        else:
            genjetList = "Gen Partons\n"
            for i in range(len(genjet_eta)):
                genjetList += "pdg %i, pt=%.2f, eta=%.3f, phi=%.3f\n"%(genjet_pid[i], genjet_pt[i], genjet_eta[i],genjet_phi[i])
            if len(genjet_eta)==0:
                genjetList = "No Gen Partons"

        ax.text(0.5, 0.85, genjetList,
                verticalalignment='top', horizontalalignment='center',
                transform=ax.transAxes,
                fontsize=12)


    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.set_ylabel(r'$\phi$')
    ax.set_xlabel(r'$\eta$')

    if not label==None:
        ax.text(0.5, 0.05, label,
                verticalalignment='bottom', horizontalalignment='center',
                transform=ax.transAxes,
                fontsize=12)

    for ijet, jet in enumerate(jets):
        ax.text(jet.eta/6.+0.5, jet.phi/(2*np.pi)+0.5, "%.2f"%jet.pt,
                verticalalignment='center', horizontalalignment='center',
                transform=ax.transAxes,
                fontsize=10)
    ax.legend(["matched gen", "unmatched gen"],loc=9)
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_yticklabels() for a in fig.axes[1:]], visible=False)
    fig.tight_layout()
    fig.savefig(name)

#drawJets(1, 6, fig, label, jets, df.tower_eta, df.tower_phi, df.tower_pt, dfGen.genjet_eta, dfGen.genjet_phi)
