# for data manipulation
import pandas as pd

# importing for scientific and numeric manipulations
from scipy.stats import kde
import numpy as np

# importing components for visualization
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import PatchCollection

# for suppressing warnings
import warnings

from pydiaid.diapasef.method_evaluation import boxes

# importing components for visualization
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams.update({'font.size': 16.5})
mpl.rcParams['pdf.fonttype'] = 42

# for suppressing warnings
warnings.filterwarnings('ignore')


def kernel_density_calculation(
    library_subset: pd.DataFrame,
    nbins: int,
) -> np.ndarray:
    """Calculate the kernel density estimation of a data frame representing a
        filtered proteomics library or single-shot measurement.

    Parameters:
    library_subset (pd.DataFrame): pre-filtered data frame with unified column
        names.
    nbins (int): number of bins for the kernel density estimation.

    Returns:
    xi, yi, zi (numpy.ndarray): coordinates of the kernel density estimation,
        where zi indicates the density.
    """
    x = library_subset['mz']
    y = library_subset['IM']

    k = kde.gaussian_kde([x, y])
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
    """Create a density plot from a data frame representing a filtered
        proteomics library or single-shot measurement.

    Parameters:
    xi, yi, zi (numpy.ndarray): coordinates of the kernel density estimation,
        where zi indicates the density.
    plot_parameters (dict): dictionary, which contains all input parameters for
        creating plots (e.g., displayed m/z-range, ion mobility-range).
    file_name: file path and file name where the plot should be saved
    gui (bool): whether to use in the GUI or not. Defaults is False.
    """
    fig, ax = plt.subplots()

    plt.pcolormesh(
        xi, yi, zi.reshape(xi.shape),
        vmax=0.02,
        cmap=plt.cm.viridis
    )
    ax.set_xlim(plot_parameters["plot_mz"][0], plot_parameters["plot_mz"][1])
    ax.set_ylim(plot_parameters["plot_IM"][0], plot_parameters["plot_IM"][1])
    ax.set_facecolor('#440256')
    plt.xlabel('$\mathregular{\it{m/z}}$')
    plt.ylabel('$\mathregular{1/K_0}$ [Vs $\mathregular{cm^{-2}}$]')
    plt.colorbar().set_label('Density', labelpad=15, rotation=270)
    plt.savefig(file_name, bbox_inches='tight', pad_inches=0.1, dpi=300)
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
        creating plots (e.g., displayed m/z-range, ion mobility-range).
    file_name: file path and file name where the plot should be saved
    gui (bool): whether to use in the GUI or not. Defaults is False.
    """
    fig = plt.figure()

    plt.hist(
        [
            library_subset['mz'][library_subset['Charge'] == 3],
            library_subset['mz'][library_subset['Charge'] == 2],
            library_subset['mz']
        ],
        bins=100,
        histtype='step',
        fill=True,
        alpha=0.5,
        color=['#267FA5', '#4EA7BB', '#7AC7C9'],
        label=[
            'triply charged precursors',
            'doubly charged precursors',
            'all precursors'
        ]
    )

    plt.xlim([plot_parameters["plot_mz"][0], plot_parameters["plot_mz"][1]])
    plt.ylabel('No. of precursors')
    plt.legend(bbox_to_anchor=(0.8, 1.42))
    plt.xlabel('$m/z$')
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
    """Plot the dia-PASEF acquisition scheme on top of a density plot
    representing a filtered proteomics library or single-shot measurement.

    Parameters:
    df (pd.DataFrame): data frame that contains the scan type (PASEF), scan
        number and the corresponding di-aPASEF window coordinates for each
        window per scan.
    xi, yi, zi (numpy.ndarray): coordinates of the kernel density estimation,
        where zi indicates the density.
    plot_parameters (dict): dictionary, which contains all input parameters for
        creating plots (e.g., displayed m/z-range, ion mobility-range).
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
    plt.colorbar().set_label('Density', labelpad=15, rotation=270)
    plt.savefig(file_name, bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.clf()


def kernel_density_calculation_multiple_charge(
    library_subset: pd.DataFrame,
    nbins: int,
    charge_setting: int,
) -> np.ndarray:
    """Calculate the kernel density estimation for specific multiple charged
        precursors.

    Parameters:
    library_subset (pd.DataFrame): pre-filtered data frame with unified column
        names.
    nbins (int): number of bins for the kernel density estimation.

    Returns:
    xi, yi, zi (numpy.ndarray): coordinates of the kernel density estimation,
        where zi indicates the density.
    """
    x = library_subset['mz'][library_subset['Charge'] == charge_setting]
    y = library_subset['IM'][library_subset['Charge'] == charge_setting]

    k = kde.gaussian_kde([x, y])
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    return xi, yi, zi


def generate_dia_boxes(
    df_parameters_final, 
    im_steps
    ):
    """
    Generates boxes for dia-PASEF scans. Each box represents the coordinates of the quadrupole 
    per TOF trigger. This function adapts the synchro-PASEF box generation for dia-PASEF's 
    rectangular windows.
    
    Parameters:
    df_parameters_final (pd.DataFrame): Method parameters from dia-PASEF method file. Should include
        columns: 'MS Type', 'Cycle Id', 'Start IM', 'End IM', 'Start Mass', 'End Mass'
    im_steps (int): The number of steps (boxes) within each ion mobility range
    
    Returns:
    list: A list of boxes, where each box is a list containing:
        [mz_start, im_position, mz_width, im_step_size, cycle_id]
    pd.DataFrame: Filtered dataframe containing only PASEF scan entries
    """
    # Filter out MS1 scans and reset index
    df_temp = df_parameters_final[df_parameters_final["MS Type"].str.contains("PASEF", case=False)].copy()
    df_temp = df_temp.reset_index(drop=True)
    
    boxes = []
    
    # Calculate ion mobility step size
    im_range = df_temp["End IM"].max() - df_temp["Start IM"].min()
    general_im_step_size = im_range / im_steps

    for index, row in df_temp.iterrows():
        # Calculate m/z width for this window
        mz_width = row["End Mass"] - row["Start Mass"]
        im_height = row["End IM"] - row["Start IM"]
        num_im_steps = round(im_height/general_im_step_size)
        individual_im_step_size = im_height/num_im_steps
        
        # Generate boxes for each step in the ion mobility range
        for step in range(1, num_im_steps+1):
            # Calculate current ion mobility position
            im_position = row["End IM"] - (step * individual_im_step_size)
            
            # In dia-PASEF, mz_start is constant within each window
            mz_start = row["Start Mass"]
    
            boxes.append([
                mz_start,           # m/z start position
                im_position,        # current IM position
                mz_width,          # m/z width of window
                individual_im_step_size,      # IM step size
                row["Cycle Id"]    # cycle ID
            ])

    
    return boxes, df_temp