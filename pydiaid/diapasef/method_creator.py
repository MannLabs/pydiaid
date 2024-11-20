# for data manipulation
import pandas as pd

# importing for scientific and numeric manipulations
import numpy as np

import os

# import pydiaid.oadia.method_generator as method_generator


def method_creation(
    library_mz_values: list,
    method_parameters: dict,
    A1: float,
    A2: float,
    B1: float,
    B2: float,
) -> pd.DataFrame:
    """This function calculates the dia-PASEF window coordinates based on the
        input parameters with dynamic isolation widths, adjusted to the
        precursor density in the m/z dimension. All dia-PASEF windows together
        add up to the dia-PASEF window scheme.

    Parameters:
    library_mz_values (list): list of all mz precursor values.
    method_parameters (dict): dictionary, which contains all input parameters
        for creating a dia-PASEF method (e.g., m/z-range, ion mobility-range,
        number of dia-PASEF scans, number of ion mobility steps, the overlap of
        the diaPASEF windows, the shift of the final window scheme in the ion
        mobility axis to account for the blind spot of the changing quadrupole
        at the top of each dia-PASEF window).
    A1 (float): bottom left y-coordinate of the rectangle, which defines the
        scan area and the slope of the dia-PASEF scan line.
    A2 (float): bottom right y-coordinate of the rectangle, which defines the
        scan area and the slope of the diaP-ASEF scan line.
    B1 (float): top left y-coordinate of the rectangle, which defines the scan
        area and the slope of the dia-PASEF scan line.
    B2 (float): top right y-coordinate of the rectangle, which defines the scan
        area and the slope of the dia-PASEF scan line.

    Returns:
    pd.DataFrame: data frame that contains the scan type (PASEF), scan number,
        and the corresponding diaPASEF window coordinates for each window per
        scan.
    """

    # No. of diaPASEF windows
    num_splits = method_parameters["im_steps"] * method_parameters["num_dia_pasef_scans"]

    # calculate the start and end mz range for each diaPASEF window
    # [400,420][420,450], ... :
    x_splits = divide_mz_range_in_equally_filled_parts(
        library_mz_values,
        method_parameters["mz"],
        num_splits
    )
    # x_splits = method_generator.create_variable_bins(
    #     library_mz_values, 
    #     method_parameters["mz"],
    #     num_splits,
    #     # method_parameters["min_width"],
    #     max_width = method_parameters["max_width"]
    # )
    # calculate the position of each diaPASEF window. 1st scan: position 0, 12;
    # 2nd scan: position 1, 13:
    window_position = [[
        scan,
        method_parameters["num_dia_pasef_scans"] *
        (method_parameters["im_steps"]-1) + scan]
        for scan in range(0, method_parameters["num_dia_pasef_scans"], 1)
    ]

    m_top, c_top = line_equation(
        (method_parameters["mz"][0], B1),
        (method_parameters["mz"][1], B2)
    )
    m_bottom, c_bottom = line_equation(
        (method_parameters["mz"][0], A1),
        (method_parameters["mz"][1], A2)
    )

    results_sorted = calculate_diaPASEF_window_coordinates(
        method_parameters["im_steps"],
        method_parameters["num_dia_pasef_scans"],
        window_position,
        m_bottom,
        c_bottom,
        m_top,
        c_top,
        x_splits
    )

    return enhance_diaPASEF_windows_to_edges_of_IM_ranges(
        results_sorted,
        method_parameters["ion_mobility"],
        method_parameters["overlap"]
    )


def calculate_diaPASEF_window_coordinates(
    num_steps: int,
    num_scans: int,
    window_position: list,
    m_bottom: float,
    c_bottom: float,
    m_top: float,
    c_top: float,
    x_splits: list,
) -> pd.DataFrame:
    """This function calculates the dia-PASEF window parameters. The x start and
        end coordinate of each window (= rectangle) is already defined in
        x_splits, and the y coordinates are calculated through the intersection
        point of the x coordinates with the respective line equation.

    Parameters:
    num_steps (int): number of ion mobility windows per dia-PASEF scan.
    num_scans (int): number of diaPASEF scans.
    window_position (list(int)): position of dia-PASEF windows for each dia-PASEF
        scan. E.g., 3 dia-PASEF scans * 2 ion mobility windows per dia-PASEF scan
        = 6 dia-PASEF windows: [[0,3],[1,4],[2,5]], the first list describes the
        window positions of the first dia-PASEF scan.
    m_bottom (float): the slope of the line equation formed by A1 and A2;
        corresponds to the slope of the bottom dia-PASEF window border of the
        last ion mobility windows per dia-PASEF scan.
    c_bottom (float): y-axis intersection of the line equation formed by A1 and
        A2.
    m_top (float): the slope of the line equation formed by B1 and B2;
        corresponds to the slope of the top dia-PASEF window border of the first
        ion mobility windows per dia-PASEF scan.
    c_top (float): y-axis intersection of the line equation formed by B1 and B2.
    x_splits (list(float)): start and end m/z value of each dia-PASEF window
        sorted after dia-PASEF window position.
        E.g., [[408.2423383, 484.3129605], [484.3129605, 506.2626105],
        [506.2626105, 522.7934549], ...]

    Returns:
    pd.DataFrame: data frame containing the scan type (PASEF), scan number, and
        the corresponding dia-PASEF window coordinates for each window per scan.
    """
    y_increments = dict()
    results = list()
    scan_bottom = dict()

    for i in range(0, num_steps, 1):
        for scan in range(0, num_scans, 1):
            temp = list()
            if i == 0:
                x_bottom = x_splits[window_position[scan][0]][0]
                y_bottom = m_bottom*x_bottom + c_bottom
                x_top = x_splits[window_position[scan][1]][-1]
                y_top = m_top*x_top + c_top
                y_increments[scan] = (y_top - y_bottom) / num_steps
                scan_bottom[scan] = y_bottom

            y_bottom = scan_bottom[scan]
            y_upper = y_bottom + y_increments[scan]

            temp.append("PASEF")
            temp.append(scan + 1)
            temp.append(round(y_bottom, 2))
            temp.append(round(y_upper, 2))
            temp.append(round(x_splits[i*num_scans + scan][0], 2))
            temp.append(round(x_splits[i*num_scans + scan][-1], 2))
            scan_bottom[scan] = y_upper
            results.append(temp)

    results_sorted_decreasing_mz_value = sorted(
        results,
        key=lambda x: float(x[4]),
        reverse=True
    )
    results_sorted_increasing_diaPASEF_scan_No = sorted(
        results_sorted_decreasing_mz_value,
        key=lambda x: int(x[1])
    )
    return results_sorted_increasing_diaPASEF_scan_No


def enhance_diaPASEF_windows_to_edges_of_IM_ranges(
    df_without_column_names: pd.DataFrame,
    ion_mobility: tuple,
    overlap: float,
) -> pd.DataFrame:
    """This function recalculates the upper border of the first ion mobility
        windows per dia-PASEF scan and the lower border of the last ion mobility
        windows per dia-PASEF scan to enhance the dia-PASEF windows to the
        borders of the ion mobility range. This enhances the overall precursor
        coverage without loss in sensitivity or acquisition speed.

    Parameters:
    df_without_column_names (pd.DataFrame): data frame containing the scan type
        (PASEF), scan number, and the corresponding dia-PASEF window coordinates
        for each window per scan.
    ion_mobility (tuple): lower and upper value of the ion mobility range.
    overlap (float): this value indicates how much the dia-PASEF windows should
        overlap in the m/z dimension (unit: Da).

    Returns:
    pd.DataFrame: data frame that contains the scan type (PASEF), scan number,
        and the corresponding dia-PASEF window coordinates for each window per
        scan extended to the ion mobility borders.
    """
    df = pd.DataFrame(
        df_without_column_names,
        columns = ["MS Type", "Cycle Id", "Start IM", "End IM", "Start Mass", "End Mass"]
    )

    count = dict()
    for item in df["Cycle Id"].values:
        count[item] = count.get(item, 0)+1

    results = list()
    for i in count.keys():
        df_temp = df[df["Cycle Id"] == i]
        df_temp["End IM"].iloc[0] = ion_mobility[1]
        df_temp["Start IM"].iloc[count[i]-1] = ion_mobility[0]
        df_temp["End Mass"] = df_temp["End Mass"].apply(
            lambda x: float(x) + overlap
        )
        results.append(df_temp)
    results_concatenated = pd.concat(
        results,
        ignore_index=True
    ).values.tolist()

    return pd.DataFrame(
        results_concatenated,
        columns=["MS Type", "Cycle Id", "Start IM", "End IM", "Start Mass", "End Mass"]
    )


def divide_mz_range_in_equally_filled_parts(
    x_series: list,
    mz_range: tuple,
    num_splits: int,
) -> list:
    """This function calculates the lower and upper m/z value of each ion
        mobility window equalizing the number of precursors per diaPASEF window
        across the m/z dimension.

    Parameters:
    x_series (list(float)): list of all m/z precursor values.
    mz_range (tuple): lower and upper value of the m/z range.
    num_splits (int): number of dia-PASEF windows. Calculated by multiplying the
        number of dia-PASEF scans with the number of ion mobility windows per
        diaPASEF scan.

    Returns:
    list(float): start and end m/z value of each dia-PASEF window sorted after
        diaPASEF window position. E.g., [[408.2423383, 484.3129605],
        [484.3129605, 506.2626105], [506.2626105, 522.7934549], ...]
    """
    # filter for precursors within the m/z range
    x_filtered = [ele for ele in sorted(list(x_series.values)) if mz_range[0] <= ele <= mz_range[1]]
    # reminder of division, No. of precursors, which are hard to distributed
    # equally, [0,1,2,3,4,5,6] num_split = 3 -> reminder: 1
    rem = len(x_filtered) % num_splits
    # to distribute precursors equally; the last m/z value is added to the
        # list until the division is without reminder.-> [0, 1, 2, 3, 4, 5, 6, 6, 6]
    x_filtered += [x_filtered[-1] for i in range(num_splits - rem)]

    results = list()
    for arr in np.split(np.array(x_filtered), num_splits):
    # split the list [0, 1, 2, 3, 4, 5, 6, 6, 6] in num_split parts (here: 3):
        # [0 1 2], [3 4 5], [6 6 6]; save upper and lower m/z border in results
        # -> [0,2],[3,5],[6,6]
        results.append(
            [
                arr[0],
                arr[-1]] if len(results) == 0 else [results[-1][1], arr[-1]
            ]
        )
    return results


def line_equation(
    l1: tuple,
    l2: tuple,
) -> float:
    """This function calculates the line equation of two coordinates.
    Parameters:
    l1 (tuple): x and y coordinate of the 1st point. Line encoded as l=(x,y).
    l2 (tuple): x and y coordinate of the 2nd point.

    Returns:
    m (float): the slope of the line equation.
    c (float): y-axis intersection of the line equation formed by l1 and l2.
    """
    m = float((l2[1] - l1[1])) / float(l2[0] - l1[0])
    c = (l2[1] - (m * l2[0]))
    return m, c


def create_parameter_dataframe(
    df_parameters_final: pd.DataFrame,
    file_name: str,
) -> None:
    """This function writes the dia-PASEF parameter file as .txt. This file
        serves as input for the generation of new dia-PASEF methods with
        timsControl (timsTof Pro 2 control software).
    Parameters:
    df_parameters_final (pd.DataFrame): data frame that contains the scan type
        (PASEF), scan number, and the corresponding dia-PASEF window coordinates
        for each window per scan extended to the ion mobility borders.
    file_name (str): file path and file name where the dia-PASEF parameter file
        should be saved.
    """
    results_input = df_parameters_final.astype(str)
    results = results_input.values.tolist()

    with open(os.path.join(os.path.dirname(__file__), "static/DIAParameterspy3TC.txt")) as f:
        content = f.readlines()

    ref_line = [
        len(col.replace("\n",""))
        for col in content[5].split(',')
    ][:-1]
    content_new = list()
    for result in results:
         content_new.append(
            ",".join([result[i].rjust(ref_line[i])
            for i in range(len(result))] + [content[5].split(',')[-1]])
        )
    final_content = content[:4] + content_new

    with open(file_name, 'w') as filehandle:
        filehandle.writelines(final_content)
