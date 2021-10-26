# for data manipulation:
import pandas as pd

# importing for scientific and numeric manipulations
from scipy.stats import kde
import numpy as np

# importing components for visualization
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import PatchCollection
import seaborn as sns

# for suppressing warnings
import warnings

from diaid_pasef.method_evaluation import boxes

# importing components for visualization
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams.update({'font.size': 16.5})

# for suppressing warnings
warnings.filterwarnings('ignore')


def kernel_density_calculation(
    library_subset: pd.DataFrame,
    nbins: int,
) -> np.ndarray:  # todo: how to describe multiple values?
    """Calculates the kernel density estimation of a data frame representing a
        filtered proteomics library or
        single-shot measurement.

    Parameters:
    library_subset (pd.DataFrame): pre-filtered data frame with unified column
        names.
    nbins (int): number of bins for the kernel density estimation.

    Returns:
    xi, yi, zi (numpy.ndarray): coordinates of the kernel density estimation
        where zi indicates the density.
    """
    x = library_subset['mz']
    y = library_subset['IM']

    k = kde.gaussian_kde([x,y])
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    return xi, yi, zi


def plot_density(
    xi: np.ndarray,
    yi: np.ndarray,
    zi: np.ndarray,
    plot_parameters: dict,
    file_name: str,
    gui: bool = False
) -> None:
    """Creates a density plot from a data frame representing a filtered
        proteomics library or single-shot measurement.

    Parameters:
    xi, yi, zi (numpy.ndarray): coordinates of the kernel density estimation
        where zi indicates the density
    plot_parameters (dict): dictionary, which contains all input parameters for
        creating plots (e.g., displayed m/z-range, ion mobility-range)
    file_name: file path and file name where the plot should be saved
    gui (bool): whether to use in the GUI or not. Defaults is False.
    """
    fig, ax = plt.subplots()
    plt.pcolormesh(
        xi, yi, zi.reshape(xi.shape),
        vmax = 0.02,
        cmap=plt.cm.viridis
    )
    ax.set_xlim(plot_parameters["plot_mz"][0], plot_parameters["plot_mz"][1])
    ax.set_ylim(plot_parameters["plot_IM"][0], plot_parameters["plot_IM"][1])
    ax.set_facecolor('#440256')
    plt.xlabel('$\mathregular{\it{m/z}}$')
    plt.ylabel('$\mathregular{1/K_0}$ [Vs $\mathregular{cm^{-2}}$]')
    plt.colorbar().set_label('Density', labelpad=-28, y=1.14, rotation=0)
    plt.savefig(file_name, bbox_inches='tight', pad_inches=0, dpi=300)
    if gui:
        return fig
    else:
        plt.clf()


def plot_precursor_distribution_as_histogram(
    library_subset: pd.DataFrame,
    plot_parameters: dict,
    file_name: str,
    gui: bool = False
) -> None:
    """Plot histogram with the precursor distribution in the m/z dimension
        sorted by charge state.

    Parameters:
    library_subset (pd.DataFrame): pre-filtered data frame with unified column
        names.
    plot_parameters (dict): dictionary, which contains all input parameters for
        creating plots (e.g., displayed m/z-range, ion mobility-range)
    file_name: file path and file name where the plot should be saved
    gui (bool): whether to use in the GUI or not. Defaults is False.
    """
    fig = plt.figure()
    colors = ['#7AC7C9', '#4EA7BB', '#267FA5']
    sns.distplot(
        library_subset['mz'],
        hist=True,
        kde=False,
        rug=False,
        color=colors[0],
        label='all precursors'
    )
    sns.distplot(
        library_subset['mz'][library_subset['Charge'] == 2],
        hist=True,
        kde=False,
        rug=False,
        color=colors[1],
        label='doubly charged precursors'
    )
    sns.distplot(
        library_subset['mz'][library_subset['Charge'] == 3],
        hist=True,
        kde=False,
        rug=False,
        color=colors[2],
        label='triply charged precursors'
    )
    plt.xlim([plot_parameters["plot_mz"][0], plot_parameters["plot_mz"][1]])
    plt.ylabel('No. of precursors')
    plt.legend(bbox_to_anchor=(0.8, 1.42))
    plt.xlabel('$\mathregular{\it{m/z}}$') # comment
    fig.savefig(file_name, bbox_inches='tight', pad_inches=0, dpi=300)
    if gui:
        plt.legend(
            loc='upper right',
            fontsize='x-small'
        )
        return fig
    else:
        plt.clf()


def plot_density_and_method(
    df: pd.DataFrame,
    xi: np.ndarray,
    yi: np.ndarray,
    zi: np.ndarray,
    plot_parameters: dict,
    file_name: str,
) -> None:
    """Plot histogram with the precursor distribution in the m/z dimension
        sorted by charge state.

    Parameters:
    df (pd.DataFrame): data frame that contains the scan type (PASEF), scan
        number and the corresponding diaPASEF window coordinates for each
        window per scan.
    xi, yi, zi (numpy.ndarray): coordinates of the kernel density estimation
        where zi indicates the density.
    plot_parameters (dict): dictionary, which contains all input parameters for
        creating plots (e.g., displayed m/z-range, ion mobility-range)
    file_name: file path and file name where the plot should be saved.
    """

    rect, _, _, _, _ = boxes(df)

    fig, ax = plt.subplots()
    plt.pcolormesh(
        xi, yi, zi.reshape(xi.shape),
        vmax=0.02,
        cmap=plt.cm.viridis
    )
    pc = PatchCollection(
        rect,
        facecolor=plot_parameters["window_color"],
        alpha=plot_parameters["window_transparency"],
        edgecolor=plot_parameters["window_frame_color"]
    )
    ax.add_collection(pc)

    ax.set_xlim(plot_parameters["plot_mz"][0], plot_parameters["plot_mz"][1])
    ax.set_ylim(plot_parameters["plot_IM"][0], plot_parameters["plot_IM"][1])
    ax.set_facecolor('#440256')
    plt.xlabel('$\mathregular{\it{m/z}}$')
    plt.ylabel('$\mathregular{1/K_0}$ [Vs $\mathregular{cm^{-2}}$]')
    plt.colorbar().set_label('Density', labelpad=-28, y=1.14, rotation=0)
    plt.savefig(file_name, bbox_inches='tight', pad_inches=0, dpi=300)
    plt.clf()
