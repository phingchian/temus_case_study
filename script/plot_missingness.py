import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def plot_hist(ax, filepath, col, title='', ylabel='Frequency', xlabel='Proportion', xlim=None, apply_log=False):
    """
    Plot a histogram of the selected column from a file

    Args:
        ax (matplotlib.axes._subplots.AxesSubplot): The axis to plot on.
        filepath (str): Path to the input csv file.
        col (str): Column name to plot
        title (str, optional): Title of the plot. Defaults to ''.
        ylabel (str, optional): Label for y-axis. Defaults to 'Frequency'.
        xlabel (str, optional): Label for x-axis. Defaults to 'Proportion'.
        xlim (_type_, optional): Limits for the x-axis. Defaults to None.
        apply_log (bool, optional): Whether to apply log transformation to the data. Defaults to False.
    """
    data = pd.read_csv(filepath, sep='\\s+')[col]
    if apply_log:
        data = np.log10(data)
    ax.hist(data, bins=20)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_xlim(xlim)

def main(smiss_path, vmiss_path, outfile='missingness_plot.png'):
    """Main function to plot missingness histograms for individual and SNP missingness.

    Args:
        smiss_path (str): Path to the .smiss file.
        vmiss_path (str): Path to the .vmiss file.
        outfile (str, optional): Output filename. Default to 'missingness_plot.png'.
    """
    labels = ['Individual', 'SNP']
    _, ax = plt.subplots(1, 2, figsize=(12, 6))
    for i, filepath in enumerate([smiss_path, vmiss_path]):
        plot_hist(ax[i], filepath, 'F_MISS', title=f'{labels[i]} Missingness', xlim=(0, 1))
    plt.tight_layout()
    plt.savefig(outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot missingness histograms from PLINK2 --missing output files.")
    parser.add_argument('smiss_path', type=str, help="Path to the .smiss file")
    parser.add_argument('vmiss_path', type=str, help="Path to the .vmiss file")
    parser.add_argument('outfile', type=str, help="Output filename")
    args = parser.parse_args()
    
    main(args.smiss_path, args.vmiss_path, args.outfile)