# for data manipulation
import pandas as pd

# importing for scientific and numeric manipulations
import numpy as np

# importing components for visualization
from matplotlib.patches import Rectangle


def calculate_precursor_within_scan_area(
    library: pd.DataFrame,
    mz: tuple,
    im: tuple,
) -> dict:
    """This function calculates the number of precursors within the scan area
    in percent.

    Parameters:
    library (pd.DataFrame): a pre-filtered data frame with unified column names
        containing all required precursor information.
    mz (tuple): lower and upper value of the m/z range.
    im (tuple): lower and upper value of the ion mobility range.

    Returns:
    dict: this dictionary has only one item. The key is the value description,
        and the value is the ratio of precursors [%] within the set scan area.
    """
    filtered_selection = library['mz'] > mz[0]
    filtered_selection &= library['mz'] < mz[1]
    filtered_selection &= library['IM'] > im[0]
    filtered_selection &= library['IM'] < im[1]

    prec_scan_area = (len(library[filtered_selection])/len(library['mz']))*100

    return {'precursors within m/z-range [%]': round(prec_scan_area, 2)}


def calculate_percentage_multiple_charged_precursors(
    library_subset: pd.DataFrame,
) -> dict:
    """This function calculates the ratio of precursors with different charge
        states in percent.

    Parameters:
    library_subset (pd.DataFrame): a pre-filtered data frame with unified
        column names containing all required precursor information.

    Returns:
    dict: The key is the value description, and the value is the ratio of
        precursors [%] with a specific charge state.
    """
    dict_charge_prec = dict()
    for i in range(1, 5):
        dict_charge_prec[str(i)] = round(
            (
                len(
                    library_subset[library_subset['Charge'] == i]
                ) / len(library_subset)
            ) * 100, 2
        )
    return dict_charge_prec


def coverage(
    df: pd.DataFrame,
    library: pd.DataFrame,
) -> dict:
    """This function calculates different evaluation values to estimate the
        optimization potential of the tested dia-PASEF method based on the
        theoretically covered (doubly/ triply charged) precursors.

    Parameters:
    df (pd.DataFrame): data frame that contains the scan type (PASEF), scan
        number, and the corresponding dia-PASEF window coordinates for each
        window per scan.
    library_subset (pd.DataFrame): a pre-filtered data frame with unified
        column names containing all required precursor information.

    Returns:
    dict: This dictionary describes multiple values to evaluate diaPASEF
        methods in respect to the proteomics library.
    """
    _, borders, small_window, big_window, mean_window = boxes(df)

    coverage_list = list()
    for i in borders:
        # subset for precursors within boundaries and append to dataset:
        # mz_min:i[0], mz_max:i[2], im_min:i[1], im_max:i[3]
        filtered_selection = library['mz'] >= i[0]
        filtered_selection &= library['mz'] <= i[2]
        filtered_selection &= library['IM'] >= i[1]
        filtered_selection &= library['IM'] <= i[3]
        coverage_list.append(library[filtered_selection])
    coverage_unique_prec = pd.concat(coverage_list, ignore_index=True)
    coverage_removed_prec_duplicates = coverage_unique_prec.drop_duplicates(['Peptide', 'Charge'])

    dict_prec_coverage = {
        "unique proteins in the library": len(library.drop_duplicates(['Proteins'])),
        "unique precursors in the library": len(library.drop_duplicates(['Peptide', 'Charge'])),
        "smallest diaPASEF window": round(small_window, 2),
        "biggest diaPASEF window": round(big_window, 2),
        "average diaPASEF window size": round(mean_window, 2),
        "No. of covered proteins": len(coverage_removed_prec_duplicates.drop_duplicates(['Proteins'])),
        "No. of covered precursors": len(coverage_removed_prec_duplicates),
    }
    dict_prec_coverage["all proteins covered"] = format(np.round(dict_prec_coverage["No. of covered proteins"]/len(library.drop_duplicates(['Proteins']))*100, 1)) + "%"
    dict_prec_coverage["all precursors covered"] = format(np.round(dict_prec_coverage["No. of covered precursors"]/len(library)*100, 1)) + "%"
    try:
        dict_prec_coverage["No. of covered, doubly charged precursors"] = len(coverage_removed_prec_duplicates[coverage_removed_prec_duplicates['Charge'] == 2])
        dict_prec_coverage["all doubly charged precursors covered"] = format(np.round(dict_prec_coverage["No. of covered, doubly charged precursors"] / len(library[library['Charge'] == 2]) * 100, 1)) + "%"
    except Exception as e:
            dict_prec_coverage["all doubly charged precursors covered"] = 0
            print(e, "doubly charged precursours are not present") 
    try:
        dict_prec_coverage["No. of covered, triply charged precursors"] = len(coverage_removed_prec_duplicates[coverage_removed_prec_duplicates['Charge'] == 3])
        dict_prec_coverage["all triply charged precursors covered"] = format(np.round(dict_prec_coverage["No. of covered, triply charged precursors"] / len(library[library['Charge'] == 3]) * 100, 1)) + "%"
    except Exception as e:
            dict_prec_coverage["all triply charged precursors covered"] = 0
            print(e, "triply charged precursours are not present") 
    try:
        dict_prec_coverage["No. of covered, quadruply charged precursors"] = len(coverage_removed_prec_duplicates[coverage_removed_prec_duplicates['Charge'] == 4])
        dict_prec_coverage["all quadruply charged precursors covered"] = format(np.round(dict_prec_coverage["No. of covered, quadruply charged precursors"] / len(library[library['Charge'] == 4]) * 100, 1)) + "%"
    except Exception as e:
            dict_prec_coverage["all quadruply charged precursors covered"] = 0
            print(e, "quadruply charged precursours are not present") 
    try:
        dict_prec_coverage["No. of covered, singly charged precursors"] = len(coverage_removed_prec_duplicates[coverage_removed_prec_duplicates['Charge'] == 1])
        dict_prec_coverage["all singly charged precursors covered"] = format(np.round(dict_prec_coverage["No. of covered, singly charged precursors"] / len(library[library['Charge'] == 1]) * 100, 1)) + "%"
    except Exception as e:
            dict_prec_coverage["all singly charged precursors covered"] = 0
            print(e, "singly charged precursours are not present") 
    return dict_prec_coverage


def boxes(
    df: pd.DataFrame,
) -> int:
    """This function calculates dia-PASEF window specific information regarding
        the window size and creates two lists. One is for plotting the dia-PASEF
        window scheme on top of a kernel density estimation plot. The second is
        to calculate the number of proteins and precursors covered by the
        dia-PASEF window scheme.

    Parameters:
    df (pd.DataFrame): data frame that contains the scan type (PASEF), scan
        number, and the corresponding dia-PASEF window coordinates for each
        window per scan.

    Returns:
    rect (list): This list contains all information to print the dia-PASEF
        windows on top of a kernel density estimation plot.
    border (list(tuples)): The list "borders" contains a list of tuples where
        each tuple represents the edge coordinates for each dia-PASEF window.
    small_window (int): size of the smallest dia-PASEF window [Da].
    big_window (int): size of the biggest dia-PASEF window [Da].
    mean_window (int): average dia-PASEF window size [Da].
    """
    rect = list()
    borders = list()
    window_width = list()
    count = dict()

    for item in df["Cycle Id"].values:
        count[item] = count.get(item, 0)+1

    for i in count.keys():
        df_temp = df[df["Cycle Id"] == i]
        df_temp["xw"] = df_temp["End Mass"] - df_temp["Start Mass"]
        df_temp["yw"] = df_temp["End IM"] - df_temp["Start IM"]

        df_temp["borders"] = df_temp.apply(
            lambda x:
                (x["Start Mass"], x["Start IM"], x["End Mass"], x["End IM"]),
                axis=1
            )
        df_temp["rect"] = df_temp.apply(
            lambda x: Rectangle(
                (x["Start Mass"], x["Start IM"]), x["xw"], x["yw"]),
                axis=1
            )
        rect += list(df_temp["rect"].values)
        borders += list(df_temp["borders"].values)
        window_width += list(df_temp["xw"].values)
        small_window = min(window_width)
        big_window = max(window_width)
        mean_window = np.mean(window_width)
    return rect, borders, small_window, big_window, mean_window
