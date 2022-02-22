# importing components for visualization
import holoviews as hv
import plotly.express as px

# for data manipulation:
import pandas as pd

# importing for scientific and numeric manipulations
from scipy.stats import kde
import numpy as np

hv.extension("bokeh")

# ---------------------------------------------------------------------------
# example code for you, Jane :)
xi, yi, zi = kernel_density_calculation_GUI(
    library,
    method_conf["graphs"]["numbins"]
    )


qmesh = density_plot_GUI(xi, yi, zi)


plot_window_scheme_on_top_of_density(df_parameters_final, qmesh, plot_parameters)


histogram_plots(library, "mz")


histogram_plots(library, "IM")

# --------------------------------------------------------------------------


def kernel_density_calculation_GUI(
    library_subset: pd.DataFrame,
    nbins: int,
) -> np.ndarray:  # todo: how to describe multiple values?
    """This function calculates the kernel density estimation of a data frame
    representing a filtered proteomics library or single-shot measurement.

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

    # Creating density plot
    nbins=150
    k = kde.gaussian_kde([x,y])
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    positions = np.vstack([xi.ravel(), yi.ravel()])
    zi = np.reshape(k(positions).T, xi.shape)

    return xi, yi, zi


def boxes_GUI(
    df: pd.DataFrame,
) -> list:
    """This function calculates dia-PASEF window specific information. It
        creates a list, containing information for plotting the dia-PASEF
        window scheme on top of a kernel density estimation plot.

    Parameters:
    df (pd.DataFrame): data frame that contains the scan type (PASEF), scan
        number, and the corresponding diaPASEF window coordinates for each
        window per scan.

    Returns:
    rect (list): This list contains all information to print the diaPASEF
        windows on top of a kernel density estimation plot.
    """

    rect = list()

    count = dict()
    for item in df["Cycle Id"].values:
        count[item] = count.get(item, 0) + 1

    for i in count.keys():
        df_2 = df[df["Cycle Id"] == i]
        df_2["xw"] = df_2.apply(lambda x: x["End Mass"]-x["Start Mass"], axis=1)
        df_2["yw"] = df_2.apply(lambda x: x["End IM"]-x["Start IM"], axis=1)
        df_2["rect"] = df_2.apply(lambda x: [
            x["Start Mass"],
            x["Start IM"],
            x["xw"],
            x["yw"]
            ], axis=1)
        rect += list(df_2["rect"].values)

    return rect


def density_plot_GUI(
    xi: np.ndarray,
    yi: np.ndarray,
    zi: np.ndarray
):
    """Creates an interactive density plot from a data frame representing a
    filtered proteomics library or single-shot measurement.

    Parameters:
    xi, yi, zi (numpy.ndarray): coordinates of the kernel density estimation,
        where zi indicates the density
    """

    qmesh = hv.QuadMesh((xi, yi, zi)).options(
        cmap=plt.get_cmap('viridis'),
        colorbar = True,
        width = 500,
        line_width = 0.00001,
        bgcolor= '#440256',
        xlabel = "m/z, Th",
        ylabel= "Inversed IM, V·s·cm\u207B\u00B2",
        colorbar_opts = {'title': 'Density'}
    )

    return qmesh


def plot_window_scheme_on_top_of_density(
    df: pd.DataFrame,
    qmesh: hv.element.raster.QuadMesh,
    plot_parameters: dict,
):
    """Plot the dia-PASEF acquisition scheme on top of a density plot
    representing a filtered proteomics library or single-shot measurement.

    Parameters:
    df (pd.DataFrame): data frame that contains the scan type (PASEF), scan
        number and the corresponding diaPASEF window coordinates for each
        window per scan.
    qmesh (holoviews.element.raster.QuadMesh): interactive density plot, where
        the dia-PASEF windows are plotted on top.
    plot_parameters (dict): dictionary, which contains all input parameters for
        creating plots (e.g., displayed m/z-range, ion mobility-range)
    """

    rect = boxes_GUI(df)

    polys = hv.Polygons(
        [
        hv.Box(
            rectangle[0],
            rectangle[1],
            (rectangle[2], rectangle[3])
            ) for rectangle in rect
        ]
    )
    polys.opts(
        color=plot_parameters["window_color"],
        line_width=1,
        alpha=plot_parameters["window_transparency"]
        )

    qmesh.opts(colorbar_opts = {'title': 'Density'})

    plot = qmesh * polys

    plot.opts(bgcolor='#440256',
              xlabel = "m/z, Th",
              ylabel= "Inversed IM, V·s·cm\u207B\u00B2"
             )

    return plot


def histogram_plots(
    library_subset: pd.DataFrame,
    x_axis: str,
) -> None:
    """Plot histograms with the precursor distribution in the m/z or ion
        mobility dimension sorted by charge state.

    Parameters:
    library_subset (pd.DataFrame): pre-filtered data frame with unified column
        names.
    x_axis (str): the histogram can be build for the m/z- or ion Mobility
        dimesion. This value selects the correct column in the data frame.
    """

    axis_dict = {
        "mz": "m/z, Th",
        "IM": "Inversed IM, V·s·cm\u207B\u00B2",
    }

    x_axis_label = axis_dict[x_axis]


    fig = px.histogram(library_subset, x = x_axis, color = 'Charge')
    # Overlay both histograms
    fig.update_layout(barmode = 'overlay',
                      xaxis_title_text = x_axis_label, # xaxis label
                      yaxis_title_text = 'No. of precursors', # yaxis label
                     )
    # Reduce opacity to see both histograms
    fig.update_traces(opacity = 0.75)
    fig.show()
