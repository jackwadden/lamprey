import h5py, argparse
import matplotlib.pyplot as plt
import numpy as np


def main(args):

    fast5_file = h5py.File(args.hdf5_file, 'r')
    plot_count = 0
    fig, axs = plt.subplots(args.max_plots, 1, figsize=(20,30))
    for readname in fast5_file:

        signal = fast5_file[readname]['Raw']['Signal'][:]
        axs[plot_count].plot(signal)

        plot_count += 1
        if plot_count >= args.max_plots: break
    plt.savefig('signals.png')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('hdf5_file')
    parser.add_argument('--max_plots', default=5, type=int)
    args = parser.parse_args()
    main(args)
