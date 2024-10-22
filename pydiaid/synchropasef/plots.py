# for data manipulation:
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import kde
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches

from matplotlib import cm
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon


import holoviews as hv
hv.extension("bokeh")

# importing components for visualization
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams.update({'font.size': 16.5})
mpl.rcParams['pdf.fonttype'] = 42


def kernel_density_calculation(
    library_subset: pd.DataFrame,
    nbins: int,
) -> np.ndarray:  # todo: how to describe multiple values?
    """Calculates the kernel density estimation of a data frame representing a
        filtered proteomics library or single-shot measurement.

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

    # Creating density plot
    nbins = nbins
    k = kde.gaussian_kde([x, y])
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    positions = np.vstack([xi.ravel(), yi.ravel()])
    zi = np.reshape(k(positions).T, xi.shape)

    return xi, yi, zi


def density_plot(
    xi: np.ndarray,
    yi: np.ndarray,
    zi: np.ndarray,
    plot_parameters: dict,
    file_name: str,
    color_scheme_name="viridis_scheme",
    gui: bool = False,
) -> plt.figure:
    """
    Creates and saves a density plot from a kernel density estimation, for either a filtered proteomics 
    library or a single-shot measurement. The grid values (xi, yi, zi) of the density estimation form 
    the core matrix of the plot.

    It takes parameters for the plot, like m/z-range and ion mobility-range, from a given dictionary. 
    An optional parameter allows the user to select the color scheme of the plot from predefined schemes. 
    The function also allows the plot to be created within a graphical user interface (GUI) by being 
    returned as a matplotlib figure.

    Parameters:
    xi, yi, zi (numpy.ndarray): Coordinates of the kernel density estimation where zi indicates the density.
    plot_parameters (dict): Dictionary which contains all input parameters for creating plots (e.g., 
        displayed m/z-range, ion mobility-range).
    file_name (str): File path and file name where the plot should be saved.
    color_scheme_name (str, optional): The name of the color scheme to be used. Defaults to ‘viridis_scheme’.
    gui (bool, optional): A flag to determine if the plot is created within a GUI. If True, the function 
        will return the plot as a matplotlib figure. Defaults to False.

    Returns:
    plt.figure: If gui is True, the function will return the created plot as a matplotlib figure. 
    Otherwise, it will return None.
    """

    color_scheme = {
        "viridis_scheme": [plt.cm.viridis, '#440256'],
        "white_scheme": [plt.cm.binary, "white"],
    }

    fig, ax = plt.subplots()

    plt.pcolormesh(
        xi, yi, zi.reshape(xi.shape),
        vmax=0.02,
        cmap=color_scheme[color_scheme_name][0],
        shading='auto'
    )
    ax.set_xlim(plot_parameters["plot_mz"][0], plot_parameters["plot_mz"][1])
    ax.set_ylim(plot_parameters["plot_IM"][0], plot_parameters["plot_IM"][1])
    ax.set_facecolor(color_scheme[color_scheme_name][1])
    plt.xlabel('m/z')
    plt.ylabel('$\mathregular{1/K_0}$ [Vs $\mathregular{cm^{-2}}$]')
    plt.colorbar().set_label('Density', labelpad=15, rotation=270)
    plt.savefig(file_name, bbox_inches='tight', pad_inches=0.1, dpi=300)
    if gui:
        return fig
    else:
        ax


def plot_precursor_distribution_as_histogram(
    library_subset: pd.DataFrame,
    plot_parameters: dict,
    file_name: str,
    gui: bool = False
) -> None:
    """
    Generates and saves a histogram representing the distribution of precursors in the m/z dimension, 
    sorted by charge state.

    This function creates a histogram from a dataframe subset that includes m/z values for different 
    charge states of precursors. It uses a dictionary to fetch parameters for the histogram like the 
    range of the m/z axis.

    This function also allows the histogram to be created within a graphical user interface (GUI) by 
    being returned as a matplotlib figure.

    Parameters:
    library_subset (pd.DataFrame): A pre-filtered dataframe with unified column names that holds the m/z 
        values and corresponding charge states.
    plot_parameters (dict): A dictionary which contains all input parameters for creating the histogram 
        (e.g., displayed m/z-range, ion mobility-range).
    file_name (str): File path and file name where the histogram should be saved.
    gui (bool, optional): A flag to determine if the histogram is created within a GUI. If True, the 
        function will return the histogram as a matplotlib figure. Defaults to False.

    Returns:
    None or plt.figure: If gui is True, the function will return the created histogram as a matplotlib 
        figure. Otherwise, it will return None.
    """
    fig, ax = plt.subplots()

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
    plt.xlabel('m/z')
    fig.savefig(file_name, bbox_inches='tight', pad_inches=0, dpi=300)
    if gui:
        plt.legend(
            loc='upper right',
            fontsize='x-small'
        )
        return fig
    else:
        plt.clf()


def density_plot_plus_scan_area(
    xi: np.ndarray,
    yi: np.ndarray,
    zi: np.ndarray,
    plot_parameters: dict,
    file_name: str,
    scan_area_definition,
    color_scheme_name="viridis_scheme",
    gui: bool = False,
) -> plt.figure:
    """
    Creates and saves a density plot from a kernel density estimation, for either a filtered 
    proteomics library or a single-shot measurement and superimposes the plot with the scan area diagonals.

    This function creates a density plot from grid values of density estimation while considering the 
    scan area. Furthermore, it allows the user to select the color scheme of the plot from predefined 
    schemes. The function also allows for the plot to be created within a graphical user interface (GUI) 
    by being returned as a matplotlib figure.

    Parameters:
    xi, yi, zi (numpy.ndarray): Coordinates of the kernel density estimation where zi indicates the density.
    plot_parameters (dict): A dictionary which contains all input parameters for creating plots 
        (e.g., displayed m/z-range, ion mobility-range).
    file_name (str): File path and file name where the plot should be saved.
    scan_area_definition (numpy.ndarray): The values representing the scan area.
    color_scheme_name (str, optional): The name of the color scheme to be used. Defaults to ‘viridis_scheme’.
    gui (bool, optional): A flag to determine if the plot is created within a GUI. If True, the function 
        will return the plot as a matplotlib figure. Defaults to False.

    Returns:
    None or plt.figure: If gui is True, the function will return the created plot as a matplotlib figure. 
    Otherwise, returns None.
    """

    color_scheme = dict()
    color_scheme = {
        "viridis_scheme": [plt.cm.viridis, '#440256'],
        "white_scheme": [plt.cm.binary, "white"],
    }

    fig, ax = plt.subplots()

    plt.pcolormesh(
        xi, yi, zi.reshape(xi.shape),
        vmax=0.02,
        cmap=color_scheme[color_scheme_name][0],
        shading='auto'
    )

    df_top = pd.DataFrame(scan_area_definition[:2], columns=['x_top', 'y_top'])
    df_bottom = pd.DataFrame(scan_area_definition[2:], columns=['x_bottom', 'y_bottom'])
    ax.plot('x_top', 'y_top', data=df_top, color='white')
    ax.plot('x_bottom', 'y_bottom', data=df_bottom, color='white')

    ax.set_xlim(plot_parameters["plot_mz"][0], plot_parameters["plot_mz"][1])
    ax.set_ylim(plot_parameters["plot_IM"][0], plot_parameters["plot_IM"][1])
    ax.set_facecolor(color_scheme[color_scheme_name][1])
    plt.xlabel('m/z')
    plt.ylabel('$\mathregular{1/K_0}$ [Vs $\mathregular{cm^{-2}}$]')
    plt.colorbar().set_label('Density', labelpad=15, rotation=270)
    plt.savefig(file_name, bbox_inches='tight', pad_inches=0.1, dpi=300)
    if gui:
        return fig
    else:
        plt.clf()


def plot_acquisition_scheme(
    df,
    path,
    colorscheme="white_only",
    gui: bool = False,
):
    """
    Generates and saves an acquisition scheme plot from a dataframe representing a mass spectrometry 
    (MS) acquisition sequence.

    This function creates an acquisition scheme plot that allows to visualize the sequence and m/z range 
    of MS1 and MS2 scans. Each bar represents the m/z range scanned in an MS2 event.

    Parameters:
    df (pd.DataFrame): A dataframe including sequencing details.
    path (str): File path and file name where the plot should be saved.
    colorscheme (str, optional): The name of the color scheme to be used. Defaults to 'white_only'.
    gui (bool, optional): A flag to determine if the plot is created within a GUI. If True, the function 
        will return the plot as a matplotlib figure. Defaults to False.

    Returns:
    None or plt.figure: If gui is True, the function will return the created plot as a matplotlib 
    figure. Otherwise, returns None.
    """
    color_scheme = dict()
    color_scheme = {
        "white_only": "white",
        "white_grey": ["white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white"]
    }

    ms1_list = list()
    ms2_list = list()

    df_temp = df[["mass pos.1 start [m/z]", "mass pos.1 end [m/z]"]][df["type"]!="ms"].astype(float)
    for index, row in df.iterrows():
        if row["type"] == "ms":
            ms1_temp = Rectangle((min(df_temp["mass pos.1 start [m/z]"]), -(index)*10), max(df_temp["mass pos.1 end [m/z]"])-min(df_temp["mass pos.1 start [m/z]"]), 10)
            ms1_list.append(ms1_temp)
        else:
            ms2_temp = Rectangle((float(row["mass pos.1 start [m/z]"]), -(index)*10), float(row["mass pos.1 end [m/z]"])-float(row["mass pos.1 start [m/z]"]), 10)
            ms2_list.append(ms2_temp)

    plt.clf()
    ax = plt.axes([0, 0, 1, 0.12*(len(ms1_list)+len(ms2_list))])
    ms1 = PatchCollection(ms1_list, facecolors="lightyellow", alpha=0.9)
    ax.add_collection(ms1)
    ms2 = PatchCollection(ms2_list, facecolors=color_scheme[colorscheme], alpha=0.5)
    ax.add_collection(ms2)

    ax.set_facecolor('#440256')
    plt.grid(alpha=0.5)
    plt.axis('equal')
    ax.get_yaxis().set_visible(False)

    tick_distance = 30
    ticks = np.arange(min(df_temp["mass pos.1 start [m/z]"]), max(df_temp["mass pos.1 end [m/z]"])+1, tick_distance)
    ax.set_xticks(ticks)

    labels = [item.get_text() for item in ax.get_xticklabels()]
    for num in range(0,len(ticks)):
        labels[num] = str(num*tick_distance)
    ax.set_xticklabels(labels)

    ax.set_xlabel("window width [Th]")

    ms1_patch = mpatches.Patch(color="lightyellow", label='MS1 scan', alpha=0.9)
    ms2_patch = mpatches.Patch(color="white", label='MS2 scans', alpha=0.5)
    ax.legend(handles=[ms1_patch, ms2_patch], bbox_to_anchor=(1, 0, 0, 1.6), facecolor='#440256', labelcolor="white", framealpha=1)

    ax.set_title(str(len(ms1_list))+ " MS1 scans & " + str(len(ms2_list)) + " MS2 scans")
    ax.figure.savefig(path, bbox_inches='tight', pad_inches=0, dpi=300)

    if gui:
        return ax.figure
    else:
        plt.clf()


def plot_method_and_precursors(
    xi: np.ndarray,
    yi: np.ndarray,
    zi: np.ndarray,
    plot_parameters: dict,
    df,
    path,
    alpha=0.2,
    window_color = "white",
    color_scheme_name="viridis_scheme",
    gui: bool = False,
):
    """
    Generates and saves a density plot from a kernel density estimation, for either a filtered 
    proteomics library or a single-shot measurement. Additionally, it overlays the scan windows 
    on the density plot.

    This function creates a density plot with overlaid Polygones of scan windows, which provides a 
    visual representation precursor coverage. The scan windows are color-coded based on the parameters.

    Parameters:
    xi, yi, zi (numpy.ndarray): Coordinates of the kernel density estimation where zi indicates the density.
    plot_parameters (dict): A dictionary which contains all input parameters for creating the plot 
        (e.g., displayed m/z-range, ion mobility-range).
    df (pd.DataFrame): A dataframe including scan window details.
    path (str): File path and file name where the plot should be saved.
    alpha (float, optional): The opacity of the scan window Polygon. Defaults to 0.2.
    window_color (str, optional): The color of the scan window Polygon. Defaults to ‘white’.
    color_scheme_name (str, optional): The name of the color scheme to be used. Defaults to ‘viridis_scheme’.
    gui (bool, optional): A flag to determine if the plot is created within a GUI. If True, the function 
        will return the plot as a matplotlib figure. Defaults to False.

    Returns:
    None or plt.figure: If gui is True, the function will return the created plot as a matplotlib 
        figure. Otherwise, returns None.
    """
    df_temp = df[
        ["mobility pos.1 [1/K0]",
        "mass pos.1 start [m/z]",
        "mass pos.1 end [m/z]",
        "mobility pos.2 [1/K0]",
        "mass pos.2 start [m/z]"]
    ][df["type"]!="ms"].astype(float)
    polygon_coordinates = list()
    for index, row in df_temp.iterrows():
        polygon_coordinates.append(
            Polygon(
            [
            (row["mass pos.1 start [m/z]"], row["mobility pos.1 [1/K0]"]),
            (row["mass pos.1 end [m/z]"], row["mobility pos.1 [1/K0]"]),
            (row["mass pos.2 start [m/z]"]+(row["mass pos.1 end [m/z]"]-row["mass pos.1 start [m/z]"]), row["mobility pos.2 [1/K0]"]),
            (row["mass pos.2 start [m/z]"], row["mobility pos.2 [1/K0]"])
        ]
        )
        )

    color_scheme = {
        "viridis_scheme": [plt.cm.viridis, '#440256'],
        "white_scheme": [plt.cm.binary, "white"],
    }

    fig, ax = plt.subplots()

    plt.pcolormesh(
        xi, yi, zi.reshape(xi.shape),
        vmax=0.02,
        cmap=color_scheme[color_scheme_name][0],
        shading='auto'
    )
    ax.set_xlim(plot_parameters["plot_mz"][0], plot_parameters["plot_mz"][1])
    ax.set_ylim(plot_parameters["plot_IM"][0], plot_parameters["plot_IM"][1])
    ax.set_facecolor(color_scheme[color_scheme_name][1])
    plt.xlabel('m/z')
    plt.ylabel('$\mathregular{1/K_0}$ [Vs $\mathregular{cm^{-2}}$]')
    plt.colorbar().set_label('Density', labelpad=15, rotation=270)

    pc = PatchCollection(polygon_coordinates, edgecolor = window_color, facecolor = "None")
    ax.add_collection(pc)
    pc = PatchCollection(polygon_coordinates, facecolor = window_color, alpha = alpha)

    ax.add_collection(pc)

    ax.figure.savefig(path, bbox_inches='tight', pad_inches=0.1, dpi=500)
    if gui:
        return ax.figure # return fig, ax
    else:
        plt.clf()


def histogram_precursor_slicing(
    mz_position, 
    library,
    acquisition_parameters,
    df,
    path,
    gui: bool = False,
):
    """
    Generates and saves a histogram plot showing the distribution of Ion Mobility (IM) peak width for a 
    proteomics library. Overlays scan lines indicating widths for each scan line and therefore indicates
    the proportion of possibly sliced precursors.

    This function translates the IM peak widths to ms using the acquisition parameters, calculates the 
    synchro scan widths using a specified m/z posiiton and then overlays the histogram with lines 
    indicating the scan widths for each scan. It provides a direct comparison between the distribution 
    of peptide precursor peak widths and the synchro scan widths at a specific m/z, thus visualizing the 
    slicing behavior.

    Parameters:
    mz_position (float): The m/z position where synchro scan width is calculated.
    library (pd.DataFrame): The proteomic library to calculate peak width in ms for.
    acquisition_parameters (dict): Dictionary containing all the input parameters for acquisition method
    df (pd.DataFrame): Data frame, carrying the acquition method parameters such as scan cooridnates and 
        scan number.
    path (str): File path and filename where the plot should be saved.
    gui (bool, optional): A flag to determine if the plot is created within a GUI. If True, the function 
        will return the plot as a matplotlib figure. Defaults to False.

    Returns:
    None or plt.figure: If gui is True, the function will return the created plot as a matplotlib 
        figure. Otherwise, returns None.
    """
    
    df_temp = df[
            ["mobility pos.1 [1/K0]",
            "mass pos.1 start [m/z]",
            "mass pos.1 end [m/z]",
            "mobility pos.2 [1/K0]",
            "mass pos.2 start [m/z]"]
        ][df["type"]=="vista"].astype(float).reset_index(drop=True)

    number_of_scans = len(df_temp)

    mobility_equation_dict = get_scan_lines(df_temp)
    
    library_temp = library.copy()
    IM_step = (acquisition_parameters["IM"][1]-acquisition_parameters["IM"][0])/acquisition_parameters["ramp_steps"]
    library_temp['IMlength in ms'] = (library_temp['IMlength'] / IM_step) * 0.11
    median_IM_peak = np.median(library_temp['IMlength in ms'])

    fig, ax = plt.subplots()
    plt.hist(
        library_temp['IMlength in ms'],
        bins=75,
        histtype='step',
        fill=True,
        alpha=0.5,
        color='#267FA5',
        label="IM_length"
    )

    plt.ylabel('count')
    plt.xlim(0,25)
    plt.xlabel('ion mobility peak width [ms]')
    ax.axvline(x=median_IM_peak, color='black', linestyle='--', label="median IM peak width")
    colors = list(mcolors.TABLEAU_COLORS.values())*5

    scans = range(1, number_of_scans+1)
    num = 0
    for scan in scans:
        top_border = mobility_equation_dict['scan'+str(scan)+'top'][0] * mz_position + mobility_equation_dict['scan'+str(scan)+'top'][1]
        bottom_border = mobility_equation_dict['scan'+str(scan)+'bottom'][0] * mz_position + mobility_equation_dict['scan'+str(scan)+'bottom'][1]
        scan_width_in_ms = round(((top_border-bottom_border)/ IM_step) * 0.11, 3)
        ax.axvline(x=scan_width_in_ms, color=colors[num], label = "scan"+str(scan)+"= "+str(scan_width_in_ms))
        num+=1
    plt.legend()
    fig.savefig(path, bbox_inches='tight', pad_inches=0, dpi=300)

    if gui:
        return fig
    else:
        plt.clf()


def generate_boxes(
    df_parameters_final,
    im_steps
):
    """
    Generates boxes. Each box represents the coordinates of the quadrupole per TOF trigger. Collectively
    the boxes can be used to represent synchro scan diagonals.

    This function takes in the df_parameters_final (which contains the method parameters) and the number 
    of steps for ion mobility and derives the coordinates for each box.

    Each entry in the 'boxes' list defins a box, containing the mz start position, the current ion 
    mobility position, the mz width (mass range within a box), the distance between each step in ion
    mobility, and the current index (which corresponds to the synchro scan).

    Parameters:
    df_parameters_final (pd.DataFrame): Method parameters that include synchro scan information. Columns 
        should include 'mobility pos.1 [1/K0]', 'mass pos.1 start [m/z]', 'mass pos.1 end [m/z]', 
        'mobility pos.2 [1/K0]', and 'mass pos.2 start [m/z]' of MS2 scans (not MS1).
    im_steps (int): The number of steps (boxes) within the ion mobility range.

    Returns:
    list: A list of boxes, where each box is a list containing the mz start position, the current ion 
        mobility position, the mz width, the distance between each step in ion mobility, and the current 
        box/scan index.
    pd.DataFrame: A dataframe holding the method parameters including the scan area definitions, with
        MS2 scan entries only.
    """
    df_temp = df_parameters_final[["mobility pos.1 [1/K0]", "mass pos.1 start [m/z]", "mass pos.1 end [m/z]", "mobility pos.2 [1/K0]", "mass pos.2 start [m/z]"]][df_parameters_final["type"]!="ms"].astype(float)
    df_temp_reset_index = df_temp.reset_index(drop=True)
    mobility_equation_dict = get_scan_lines(df_temp_reset_index)
    
    boxes = []
    im_distance_between_steps = (df_temp_reset_index["mobility pos.2 [1/K0]"].iloc[0]-df_temp_reset_index["mobility pos.1 [1/K0]"].iloc[0])/im_steps
    for index, row in df_temp_reset_index.iterrows():
        im_position_temp = row["mobility pos.2 [1/K0]"]
        mz_width = row["mass pos.1 end [m/z]"] - row["mass pos.1 start [m/z]"]
        for step_num in range(im_steps):
            im_position_temp = im_position_temp - im_distance_between_steps
            mz_position_start = (im_position_temp - mobility_equation_dict['scan'+str(index+1)+'top'][1])/mobility_equation_dict['scan'+str(index+1)+'top'][0]
            boxes.append([mz_position_start, im_position_temp, mz_width, im_distance_between_steps, index+1])
    
    return boxes, df_temp_reset_index


def generate_gif_single_windows(
    xi,
    yi,
    zi,
    plot_parameters,
    boxes,
    scans,
    save_at,
    facecolor="#FF0098",
    window_color = None,
    scans_plotted_at_once=5,
):
    """
    Generates a grid of plots, where each plot highlights a specific area in a 2D scatter plot matrix. 
    Each time a synchro scan and a specific area of the synchro scan are highlighted. The underlying 
    scatter plot matrix shows the relation between peptide abundance and peptide position in m/z and 
    ion mobility.

    The function iterates through the synchro scan numbers and corresponding quadrupole coordinates at
    given ion mobility positions and highlights each of them in a separate plot. The color of the 
    highlighted box and the boxes of all synchro scans can be defined in the parameters.

    Note that it is possible to have several scans highlighted in one plot, as defined by the 
    'scans_plotted_at_once' parameter. Hence, the function creates a time series of plots.

    Parameters:
    xi (numpy.array): The x values of the 2D scatter plot, m/z value of peptide.
    yi (numpy.array): The y values of 2D scatter plot, ion mobility value of peptide.
    zi (numpy.array): The z values of the 2D scatter plot indicating the density.
    plot_parameters (dict): Contains the limits for the plot in m/z and ion mobility.
    boxes (list): List of coordinates of rectangles (boxes) to be plotted, earlier generated with 
        generate_boxes function.
    scans (list): Each integer value corresponds for one scan.
    save_at (str): Path string of the directory where to save the generated plots.
    facecolor (str): Color of the box that is currently being focused on.
    window_color (str): Color scheme for other scans. None defaults to white, 'colorful' and 
        'white_grey' can be used.
    scans_plotted_at_once (int): The number of scans that are being highlighted in one plot.

    Returns:
    None. Produces and saves plots to the defined path.
    """
    cnt = 1
    color_list = {
        "colorful": ["blue", "red", "green", "blue", "yellow", "orange", "red", "green", "blue", "yellow", "orange", "blue", "red", "green", "blue", "yellow", "orange", "red", "green", "blue", "yellow", "orange"],
        "white_grey": ["white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white", "silver", "white"]
    }
    plt.clf()
    for scan in scans:
        for index, rectangle in enumerate(boxes):
            if rectangle[4]==scan:
                ax = plt.axes()

                plt.pcolormesh(
                    xi, yi, zi.reshape(xi.shape),
                    vmax=0.02,
                    cmap=plt.cm.viridis,
                    shading='auto'
                )
                ax.set_xlim(plot_parameters["plot_mz"][0], plot_parameters["plot_mz"][1])
                ax.set_ylim(plot_parameters["plot_IM"][0], plot_parameters["plot_IM"][1])
                ax.set_facecolor('#440256')
                plt.xlabel('$\mathregular{\it{m/z}}$')
                plt.ylabel('$\mathregular{1/K_0}$ [Vs $\mathregular{cm^{-2}}$]')
                plt.colorbar().set_label('Density', labelpad=15, rotation=270)

                num = 0
                for item in scans:
                    box_temp = [
                            Rectangle(
                                (rect[0],
                                rect[1]),
                                rect[2], rect[3]) for rect in boxes if rect[4]==item
                        ]
                    if window_color == None:
                        box = PatchCollection(box_temp, facecolors='white', alpha=0.1, edgecolor="None")
                    else:
                        box = PatchCollection(box_temp, facecolors=color_list[window_color][num], alpha=0.5, edgecolor="None")
                        num+=1

                    ax.add_collection(box)

                filtered_boxes = [item for item in boxes[index: index + scans_plotted_at_once] if item[4]==rectangle[4]]

                box_temp = [
                        Rectangle(
                            (rect[0],
                            rect[1]),
                            rect[2], rect[3]) for rect in filtered_boxes
                ]

                pc = PatchCollection(box_temp, facecolor=facecolor)
                ax.add_collection(pc)
                ax.figure.savefig(save_at+'{0:0>5}'.format(cnt)+".png", bbox_inches='tight', pad_inches=0.1, dpi=100)
                plt.clf()
                cnt += 1


def get_scan_lines(df):
    """
    Generates linear equations for the scan lines. The equation for a scan line represents the mobility 
    and m/z position in the quadrupole per TOF trigger for a single MS2 scan. The equations are derived 
    from the scan area coordinates (mobility and m/z position) for each synchro scan.

    This function extracts the mass and mobility positions for each scan from the method parameter
    dataframe (df), calculates the linear equations of the scan lines, and stores the equations in a 
    dictionary where the keys are the scan line name (eg… 'scan1top') and the values are the equations of 
    the scan lines in the form of (m/z per mobility slope, intercept).

    Parameters:
    df (pd.DataFrame): The input dataframe that contains the ma/z and ion mobility positions for each 
        scan. Columns should include 'mobility pos.1 [1/K0]', 'mass pos.1 start [m/z]', 
        'mass pos.1 end [m/z]', 'mobility pos.2 [1/K0]', and 'mass pos.2 start [m/z]' of MS2 scans (not MS1).

    Returns:
    dict: A dictionary containing 'scanXtop' and 'scanXbottom' as keys, where X is the scan index, and 
    their corresponding values are the linear equations of the scan lines derived from the mass positions 
    and the mobility positions for each scan. The equation is represented in the form of a tuple which is
    (mass per mobility Slope, intercept).
    """
    dict_equation = dict()

    for index, row in df.iterrows():
        dict_equation[f"scan{index+1}top"] = get_slope_and_intercept(
            mz1 = row["mass pos.1 start [m/z]"],
            im1 = row["mobility pos.1 [1/K0]"],
            mz2 = row["mass pos.2 start [m/z]"],
            im2 = row["mobility pos.2 [1/K0]"],
        )
        dict_equation[f"scan{index+1}bottom"] = get_slope_and_intercept(
            mz1 = row["mass pos.1 end [m/z]"],
            im1 = row["mobility pos.1 [1/K0]"],
            mz2 = row["mass pos.2 start [m/z]"]+(row["mass pos.1 end [m/z]"]-row["mass pos.1 start [m/z]"]),
            im2 = row["mobility pos.2 [1/K0]"],
        )
    return dict_equation


def get_slope_and_intercept(mz1, im1, mz2, im2):
    """
    Calculates the slope and intercept of linear equations based on two points with ion mobility and 
    m/z coordinates.

    Parameters:
    mz1 (float): The m/z value of the first coordinate.
    im1 (float): The ion mobility of the first coordinate.
    mz2 (float): The m/z value of the second coordinate.
    im2 (float): The ion mobility of the second coordinate.

    Returns:
    tuple: A tuple containing the slope and the intercept of the linear equation.
    """
    slope = (im2 - im1) / (mz2 - mz1)
    intercept = im1 - slope * mz1
    return slope, intercept



