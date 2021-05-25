#!/usr/bin/env python3
"""Plots observation files written in this comparison."""

import numpy as np
import pandas
from matplotlib import pyplot as plt

def load(filename):
    """Loads an observation file"""
    try:
        df = pandas.read_csv(filename, comment='#', delimiter=',')
    except pandas.errors.ParserError:
        df = pandas.read_csv(filename, comment='#', delimiter=' ')

    if len(df.keys()) == 1:
        # try spaces instead?
        df = pandas.read_csv(filename, comment='#', delimiter=' ')

    if len(df.keys()) == 1:
        raise RuntimeError(f"Unable to read file {filename}, no observations present?")

    return df


def fig_setup(nplots, figsize=None, fig=None):
    ncols = int(np.ceil(np.sqrt(nplots)))
    nrows = nplots // ncols + 1

    if (figsize is None):
        figsize = (8,5)

    if fig is None:
        fig = plt.figure(figsize=figsize)
    axs = fig.subplots(nrows,ncols).flatten()
    return fig, axs

def plot(df, color, name, fig=None, axs=None):
    if axs is None:
        fig, axs = fig_setup(len(df.keys())-1, fig)

    times_label = df.keys()[0]
    times = df[times_label]
    if (times_label == 'time [s]'):
        times = times / 86400
        times_label = 'time [d]'

    for i, obs in enumerate(df.keys()[1:]):
        axs[i].plot(times, df[obs], '-', color=color, label=name)
        axs[i].set_xlabel(times_label)
        axs[i].set_ylabel(obs)

    return i, fig, axs

if __name__ == '__main__':
    import sys, os
    import colors

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("DIRECTORIES", nargs="+", type=str,
                        help="List of directories containing observation files.")

    args = parser.parse_args()
    colorlist = colors.enumerated_colors(len(args.DIRECTORIES))

    fig = None
    axs = None
    for obs_file, color in zip(obs_files, colorlist):
        if not obs_file.endswith('.dat'):
            obs_file = os.path.join(obs_file, 'observation_column_all.dat')
        df = load(obs_file)
        last_plot, fig, axs = plot(df, color, obs_file, fig=fig, axs=axs)

    axs[last_plot].legend()
    plt.show()
