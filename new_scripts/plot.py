import matplotlib.pyplot as plt


def plot(hillslope_pars, mesh_pars, fig=None, axs=None, color='k'):
    """plotting routine"""
    if fig is None:
        fig = plt.figure()
    if axs is None:
        axs = fig.subplots(2,1, sharex=True)
        fig.subplots_adjust(hspace=0)

    axs[0].plot(mesh_pars['x'], mesh_pars['z'])
    axs[0].set_xticklabels([])
    axs[0].set_ylabel('elev [m]')
    axs[0].yaxis.set_label_position("right")

    axs[1].plot(mesh_pars['x'], mesh_pars['y'])
    axs[1].set_xlabel('hillslope length [m]')
    axs[1].set_ylabel('width')
    axs[1].yaxis.set_label_position("right")