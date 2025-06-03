# for data manipulation:
import pandas as pd
import numpy as np

# for creating new folders:
import os

import math
import pydiaid.synchropasef.method_evaluator as method_evaluator


def create_folder(
    paths: list,
) -> None:
    """This function creates a folder structure, where all output files are
        printed.

    Parameters:
    paths (list): a list of all path names that should be created.
    """
    for create_item in paths:
        if os.path.exists(create_item):
            next
        else:
            os.makedirs(create_item)


def caclulate_scan_area_definitions(
    scan_width,
    library_filtered_for,
    library
):
    """
    Calculates and defines the scan area based on the provided scan width, filtered library containing 
    peptide information.

    Parameters:
    scan_width (float) : The width of the scan area.
    library_filtered_for (list): A list of filter values that are used to select entries from the 
        library. The library can be filtered for specific peptide charges.
    library (pd.DataFrame): The original DataFrame containing the library/peptide information 
        which includes columns such as 'Charge'.

    Returns:
    list: A list containing the parameters for a linear function that define the top and bottom limits 
        of the scan area.

    This function first filters the library DataFrame based on the 'Charge' column using values from the 
    library_filtered_for list. After filtering, it calculates the scan area definition based on the 
    filtered library and scan width. It then iterates over a range, shifting the scan area slightly
    left and right each time, calculating the library coverage and recording the best shift value that 
    gives the maximum coverage. Finally, it returns the optimized top and bottom limits of the scan area.

    Note that 'method_evaluator.coverage()' is an external function used to calculate the coverage of the 
    library against the scan area definition. It assumes that this function is defined elsewhere in your 
    code.
    """

    filtered_selection = (library["Charge"] == 0)
    for filter_value in library_filtered_for:
        filtered_selection |= (library["Charge"] == filter_value)
    library_filtered = library[filtered_selection]

    top_limit_IM, bottom_limit_IM, x_lin_reg, y_lin_reg = scan_area_definition(library_filtered, scan_width/2, start_mz=-200, end_mz=1760)
    print("original scan area", top_limit_IM+bottom_limit_IM)
    print("please wait")

    #finding the best scan area by shifting the scan area a bit to the left an right
    list_best_coverage = list()
    for x in range(-50, 55, 5):
        top_limit_IM_temp = [(item[0]+x, item[1]) for item in top_limit_IM]
        bottom_limit_IM_temp = [(item[0]+x, item[1]) for item in bottom_limit_IM]
        lib = method_evaluator.coverage(library, top_limit_IM_temp, bottom_limit_IM_temp)
        list_best_coverage.append(len(lib[lib["insideIM"]==True])/len(library))
    df_best_coverage = pd.DataFrame(data={"shift_in_mz": list(range(-50, 55, 5)), "coverage": list_best_coverage})
    print("coverage_tests", df_best_coverage)
    best_shift = df_best_coverage["shift_in_mz"][df_best_coverage["coverage"]==df_best_coverage["coverage"].max()].iloc[0]
    print("best_shift", best_shift)
    top_limit_IM = [(item[0]+best_shift, item[1]) for item in top_limit_IM]
    bottom_limit_IM = [(item[0]+best_shift, item[1]) for item in bottom_limit_IM]
    lib = method_evaluator.coverage(library, top_limit_IM, bottom_limit_IM)

    return top_limit_IM+bottom_limit_IM


def create_method(
    save_at,
    scan_area,
    method_parameters_scan_area,
    method_parameters_general,
    library
):
    """Generates a synchro-PASEF method file for proteomics data acquisition with a timsTof and saves 
    it at the specified location. This method file can be loaded into the instrument controlling software
    timsControl.

    Parameters:
    save_at (str): path, where the method file will be saved.

    scan_area (list): The parameters for linear functions that define the top and bottom limits of 
        the scan area.

    method_parameters_scan_area (dict): A dictionary containing various parameters of the scan area:
    'scans': Number of synchro scans.
    'window_type': Type of isolation window (equidistant, variable, pre-defined).
    'scan_ratio': Ratio of scan widths. For instance: [1,1,1] for equidistant and [1,2,1] for pre-defined.
    'scan_mode': Mode of scan (classical_synchro-PASEF, highly_accurate_synchro-PASEF or 
        individual_synchro-PASEF).
    'window_modification': Modification of window (None, overlap, staggered).
    'window_overlap': Overlap of window (0.2 means 20% of smallest window).
    'no_of_combined_scans': Number of combined scans for a staggered and classical-synchro-PASEF scheme.
    'window_pattern': Pattern of the window for individual_synchro-PASEF.

    method_parameters_general (dict): A dictionary containing general parameters of the method,
    such as the ion miblity range as 'dict_im_limits' and the 'MS1_positions' as list.

    library (pd.DataFrame): The DataFrame containing the library/peptide information.

    Returns:
    scan_ratio, scans: Returns the modified scan ratio and number of scans which can be used for 
        further processing.

    This function uses various parameters, along with additional functions 
    scan_ratios_with_diagonal_scans(), calculate_scan_area(), generate_isolation_windows(), and 
    create_method_file() to create a method file for proteomics data acquisition. These parameters are
    the scan area limits, method parameters, general method parameters, and the peptide library.
    Note that all these functions are assumed to be defined elsewhere in the code.
    """
    dict_scan_area = {
        'top_limit_IM': scan_area[0:2],
        'bottom_limit_IM': scan_area[2:]
    }
    scans = method_parameters_scan_area['scans']
    window_type = method_parameters_scan_area['window_type']
    scan_ratio = method_parameters_scan_area['scan_ratio']
    scan_mode = method_parameters_scan_area['scan_mode']
    window_modification = method_parameters_scan_area['window_modification']
    window_overlap = method_parameters_scan_area['window_overlap']
    no_of_combined_scans = method_parameters_scan_area['no_of_combined_scans']
    window_pattern = method_parameters_scan_area['window_pattern']
    dict_im_limits = method_parameters_general["dict_im_limits"]

    if scan_mode == "individual_synchro-PASEF":
        scans = max(max(window_pattern))+1
    else:
        next

    if window_type == "equidistant":
        scan_ratio = [1] * scans

    elif window_type == "variable":
        scan_ratio = scan_ratios_with_diagonal_scans(
            library,
            scan_area[0][0],
            scan_area[1][0],
            scan_area[2][0],
            scan_area[0][1],
            scan_area[1][1],
            scans
        )

    elif window_type == "pre-definied":
        next

    df_scan_area = calculate_scan_area(
        dict_scan_area,
        dict_im_limits,
    )
    
    df_method_parameters = generate_isolation_windows(
        df_scan_area,
        scans,
        scan_ratio,
        scan_mode,
        window_modification,
        window_overlap,
        no_of_combined_scans,
        window_pattern,
    )
    
    create_method_file(
        save_at,
        method_parameters_general["MS1_positions"],
        df_method_parameters,
        scan_area,
        method_parameters_scan_area,
        method_parameters_general
    )
    
    return scan_ratio, scans


def scan_area_definition(library_temp, width_mz, start_mz=180, end_mz=1560):
    """Defines the scan area for a given peptide library.

    Parameters:
    library_temp (pd.DataFrame): DataFrame of the peptide library.
    width_mz (float): width in m/z of the scan area.
    start_mz (int, optional): starting value in m/z of linear finction that defines the scan area. 
        Default is 180.
    end_mz (int, optional): ending value in m/z of linear finction that defines the scan area. Default 
        is 1560.

    Returns:
    tuple: Returns top_limit_IM and bottom_limit_IM which are parameters for two linear funcitons defining 
    the limits of the scan area, x_lin_reg and y_lin_reg are the x and y values obtained through linear 
    regression.

    The function  first performs a linear regression on the mz and ion mobility (IM) values of the peptide 
    library. Then it calculates the scan area based on these parameters.
    Specifically, it determines the m/z values for the top and bottom limits of the scan area and the 
    corresponding ion mobility (IM) values.
    """
    X = library_temp["mz"]
    Y = library_temp["IM"]
    model = np.polyfit(X, Y, 1)
    predict = np.poly1d(model)
    x_lin_reg = range(180, 1560)
    y_lin_reg = predict(x_lin_reg)

    #Scan area
    top_start_mz = start_mz - width_mz
    top_end_mz = end_mz - width_mz
    left_IM = predict(start_mz)
    bottom_start_mz = start_mz + width_mz
    bottom_end_mz = end_mz + width_mz
    right_IM = predict(end_mz)
    top_limit_IM = [(top_start_mz, left_IM), (top_end_mz, right_IM)]
    bottom_limit_IM = [(bottom_start_mz, left_IM), (bottom_end_mz, right_IM)]

    return top_limit_IM, bottom_limit_IM, x_lin_reg, y_lin_reg


def scan_ratios_with_diagonal_scans(
    library,
    start_low,
    temp_high,
    end_low,
    im_start,
    im_end,
    scans,
    divisions=100, # 100 is a running time that is okay to wait 4-7min
):
    """Determines the scan ratios given a library, scan area boundaries and number of scans. Aim is to 
    define isolation window widths with a balanced number of precursors. For equidistant isolation 
    windows, the isolation windows have all the same width and therefore a scan ratio of [1,1,1], whereas 
    the scan_ratio can also be [1,2,1] to express a middle isolation window with double of the width.

    Parameters:
    library (pd.DataFrame): DataFrame of the peptide library.
    start_low (float): Start of the bottom m/z interval for calculating the scan ratios.
    temp_high (float): Start of the top m/z interval for calculating the scan ratios.
    end_low (float): End of the bottom m/z interval for calculating the scan ratios.
    im_start (float): Start of the ion mobility interval for calculating scan ratios.
    im_end (float): End of the ion mobility interval for calculating the scan ratios.
    scans (int): Number of scans.
    divisions (int, optional): Number of divisions/intervals for the m/z range. Default is 100.

    Returns:
    list: Window width ratio of diagonal synchro scans.

    This function aims to calculate isolation window widths where each isolation window covers the same
    number of precursors. For this, it divides the scan area in diagonal stripes (divisions), calculates 
    the precursor content/coverage of each of this stripe, and combines them according to the number of 
    scans so that they contain a balanced number of precursors. The isolation window widths are 
    normalized and returned as a scan ratio list. It uses the 'method_evaluator.coverage()' function to 
    calculate the coverage of the library against the scan area definition.
    Note that this function is defined elsewhere.
    """

    list_coverage_per_bin = list()
    list_start_low = list()

    temp_low = start_low
    step = (end_low - start_low)/divisions

    print("Please wait - this step can take several minutes")
    for num in range(divisions):
        # calculate the scan area for each division/bin
        list_start_low.append(temp_low)
        top_limit_IM = [(temp_low, im_start), (temp_high, im_end)]
        temp_low += step
        temp_high += step
        bottom_limit_IM = [(temp_low, im_start), (temp_high, im_end)]

        # calculate the peptide coverage for each division/bin
        lib = method_evaluator.coverage(library, top_limit_IM, bottom_limit_IM)
        coverage_perc=len(lib[lib["insideIM"]==True].drop_duplicates(['Peptide', 'Charge']))/len(lib.drop_duplicates(['Peptide', 'Charge']))
        list_coverage_per_bin.append(coverage_perc)

    coverage_per_bin = sum(list_coverage_per_bin)/scans
    list_mz_diff = list()
    list_temp = list()

    # define how many scans have to be defined for a equal number of precursors per division/bin
    sum_coverage = 0
    for pos in range(len(list_coverage_per_bin)):
        if sum_coverage >= coverage_per_bin:
            temp_value = list_start_low[pos] - start_low
            list_mz_diff.append(temp_value)
            list_temp.append([start_low, list_start_low[pos]])
            start_low = list_start_low[pos]
            sum_coverage = list_coverage_per_bin[pos]
        else:
            sum_coverage += list_coverage_per_bin[pos]

    temp_value = end_low - list_temp[-1][-1]
    list_mz_diff.append(temp_value)

    return list(np.array(list_mz_diff)/min(list_mz_diff))


def calculate_scan_area(
    dict_scan_area,
    dict_im_limits,
):
    """Calculates the exact scan area for a synchro-PASEF method in dependence to the ion mobility range.

    Parameters:
    dict_scan_area (dict): A dictionary containing the keys 'bottom_limit_IM' and 'top_limit_IM' with 
        respective values representing the parameters that define the bottom and top linear functions
        framing the scan area.
    dict_im_limits (dict): A dictionary containing the lower and upper ion mobility (IM) limits with the 
        keys 'low_limit_IM' and 'up_limit_IM'.

    Returns:
    pd.DataFrame: Returns a DataFrame representing the scan area bound by the lower and upper m/z and 
    ion mobility (IM) limits. The structure of the DataFrame is as follows:
    “lower_IM”: Lower ion mobility limit.
    “lower_mz_1”: Lower m/z limit for the lower ion mobility limit.
    “upper_mz_1”: Upper m/z limit for the lower ion mobility limit.
    “upper_IM”: Upper ion mobility limit.
    “lower_mz_2”: Lower m/z limit for the upper ion mobility limit.
    “upper_mz_2”: Upper m/z limit for the upper ion mobility limit.

    This function calculates the line equations of the scan area using the provided scan area dictionary 
    and the 'creating_dict' function. Then it calculates the intersection points of these equations with 
    the ion mobility limits using the 'calc_IM_intersection' function. Finally, it generates and returns 
    a DataFrame containing the calculated intersection points defining the scan area.

    Note: 'creating_dict' and 'calc_IM_intersection' functions are assumed to be defined elsewhere.
    """
    # calculate the line equations of the scan area
    dict_equation_bottom = creating_dict(dict_scan_area['bottom_limit_IM'])
    dict_equation_top = creating_dict(dict_scan_area['top_limit_IM'])


    # calculate the intersections of the scan area and the IM limits; output: df
    count = 1

    IM_1, Intersect_1 = calc_IM_intersection(dict_equation_top, dict_im_limits['low_limit_IM'])
    IM_1, Intersect_2 = calc_IM_intersection(dict_equation_bottom, dict_im_limits['low_limit_IM'])
    IM_2, Intersect_3 = calc_IM_intersection(dict_equation_top, dict_im_limits['up_limit_IM'])
    IM_2, Intersect_4 = calc_IM_intersection(dict_equation_bottom, dict_im_limits['up_limit_IM'])
    scan_area_definitions =  [
        IM_1,
        Intersect_1,
        Intersect_2,
        IM_2,
        Intersect_3,
        Intersect_4,
    ]

    df_scan_area = pd.DataFrame(scan_area_definitions).T
    df_scan_area.columns=[
        "lower_IM",
        "lower_mz_1",
        "upper_mz_1",
        "upper_IM",
        "lower_mz_2",
        "upper_mz_2",
        ]


    return df_scan_area


def check_slope(params_list):
    """Fix floating-point precision errors in parameter lists"""
    for i, params in enumerate(params_list):
        col0, col1, col2, col3, col4 = params
        slope = (col4 - col1) / (col3 - col0)
        print(slope)
        print(params)
    

def generate_isolation_windows(
    df_scan_area,
    scans,
    scan_ratios,
    scan_mode,
    window_modification,
    window_overlap,
    no_of_combined_scans,
    window_pattern,
):
    """Generates isolation window definitions for a synchro-PASEF method file.

    Parameters:
    df_scan_area (pd.DataFrame): DataFrame representing the scan area bound by the m/z and ion 
        mobility (IM) limits.
    scans (int): Number of synchro scans.
    scan_ratios (list): A list representing the ratios of the isolation widths of diagonal scans.
    scan_mode (str): Acquisition scheme type: classical_synchro-PASEF, highly_accurate_synchro-PASEF,
        individual_synchro-PASEF.
    window_modification (optional): Parameter for modifying the scan window, e.g. window overlap or 
        staggered window scheme.
    window_overlap (optional): Parameter for overlapping the scan window (0.2 means 20% of smallest window).
    no_of_combined_scans (optional): Number of combined scans for a staggered and classical-synchro-PASEF 
        scheme.
    window_pattern (optional): Pattern of the windows for an individually defined synchro-PASEF method.

    Returns:
    pd.DataFrame: Returns a DataFrame with isolation window definitions for a synchro-PASEF method file. 
    The DataFrame contains columns such as mass position start, mass position end for each IM position.

    This function starts by creating a simple acquisition scheme without considering the scan area using 
    the 'create_basic_scheme_with_scans' function. Then, using the 'window_calculator' function, it 
    calculates the m/z start and width for lower and upper ion mobility limits based on the scan area. 
    Ultimately, the result is converted into a DataFrame and returned.
    Note: 'create_basic_scheme_with_scans' and 'window_calculator' functions are assumed to be defined 
    elsewhere.
    """

    dict_isolation_windows = create_basic_scheme_with_scans(
        scans,
        scan_mode,
        window_modification,
        no_of_combined_scans,
        window_pattern,
        window_overlap
    )

    # generate window scheme
    #for lower IM limits
    mz_start_lower_IM, mz_width_lower_IM = window_calculator(
        df_scan_area,
        scan_ratios,
        dict_isolation_windows,  
        "1"
    )
    #for upper IM limits
    mz_start_upper_IM, mz_width_upper_IM = window_calculator(
        df_scan_area,
        scan_ratios,
        dict_isolation_windows,
        "2"
    )

    rounding_factor = 0
    list_method_parameters = list()
    for index in range(len(mz_start_lower_IM)):
        list_temp = [
            df_scan_area["lower_IM"].iloc[0], 
            np.round(mz_start_lower_IM[index], rounding_factor), 
            np.round(mz_start_lower_IM[index]+mz_width_lower_IM[index], rounding_factor),
            df_scan_area["upper_IM"].iloc[0],
            np.round(mz_start_upper_IM[index], rounding_factor)
        ]
        list_method_parameters.append(list_temp)
    
    # check_slope(list_method_parameters)

    df_method_parameters = pd.DataFrame(
        list_method_parameters,
        columns=[
            "mobility pos.1 [1/K0]",
            "mass pos.1 start [m/z]",
            "mass pos.1 end [m/z]",
            "mobility pos.2 [1/K0]",
            "mass pos.2 start [m/z]"
        ]
    )

    return df_method_parameters


def create_method_file(
    save_at,
    ms1_positions,
    df_method_parameters,
    scan_area,
    method_parameters,
    method_parameters_general
):
    """Generates a .txt method file for the scan mode synchro-PASEF and saves it at the specified location.

    Parameters:
    save_at (str): path, where the method file will be saved.
    ms1_positions (list): List of positions for full scan (MS1) spiked in between MS2 scans.
    df_method_parameters (pd.DataFrame): DataFrame with detailed isolation window coordinates.
    scan_area (list): Contains parameters for the top and bottom linear functions that frame the scan area.
    method_parameters (dict): A dictionary containing parameters of the acquisition method.
    method_parameters_general (dict): A dictionary containing more general parameters of the method.

    Returns:
    None

    The function generates a synchro-PASEF method file and writes the parameters used for method 
    generation also into a txt file such as scan area and method parameters. It tranforms the isolation 
    widnow coordinatess from DataFrame to a list. It adds full scans (MS1) at the specified positions and 
    finally writes all this information into a .txt file which will be saved in the specified location.
    """
    content_method_settings = list()
    content_method_settings.append("\n\n"+"#Method creation settings: \n"
                    +"#scan area: "+str(scan_area)
                    + ",\n"+"#method parameters: "+str(method_parameters)
                    + ",\n"+"#method parameters general: " + str(method_parameters_general)
                    + "\n\n")    


    results_input = df_method_parameters
    results = results_input.values.tolist()

    content_new = list()        
    for result in results:
        content_temp = "vista,"+str(result)+" \n"
        content_new.append(content_temp.replace("[",'').replace("]",''))

    count_ms1 = 0
    ms1_item = 'ms, -, -, -, -, - \n'
    for ms1 in ms1_positions:
        content_new.insert(ms1+count_ms1-1, ms1_item)
        count_ms1 += 1

    header ='type, mobility pos.1 [1/K0], mass pos.1 start [m/z], mass pos.1 end [m/z], mobility pos.2 [1/K0], mass pos.2 start [m/z] \n'       
    content_new.insert(0, header)  
    final_content = content_new + content_method_settings
    with open(save_at+'/synchroPasef.txt', 'w') as filehandle:
            filehandle.writelines(final_content)


def creating_dict(
    limit_IM_temp: list,
)-> dict:
    """Creates a dictionary with the linear equations derived from the provided coordinates in the 
    dimensions m/z and ion mobility.

    Parameters:
    limit_IM_temp (list): List of m/z and ion mobility coordinates.

    Returns:
    dict: Returns a dictionary in which the keys are the slopes and the values are the yolks of the 
        linear equations.

    This function iterates over pairs of points in the provided list, and each time uses 'lin_equ' 
    function to calculate the linear equation between those two points, then stores the equations in a 
    dictionary. This process continues until only one point is left in the list.
    Note: 'lin_equ' function is assumed to be defined elsewhere in the code.
    """
    dict_equation = dict()
    for i in range(len(limit_IM_temp)):
        if len(limit_IM_temp) > 1:
            dict_equation = lin_equ(dict_equation, limit_IM_temp[-2], limit_IM_temp[-1])
            limit_IM_temp.remove(limit_IM_temp[-1])
    return dict_equation


def calc_IM_intersection(
    dict_equation, 
    IM
):
    """Calculates the intersection point of the ion mobility (IM) limits with the linear equations that
    define the scan area to return m/z positions.

    Parameters:
    dict_equation (dict): A dictionary in which the keys are the slopes and values are the yolks of the 
        equations.
    IM (float): The ion mobility (IM) value at which the intersection point with a linear equation will 
        be calculated.

    Returns:
    tuple: Returns the ion mobility and mz intersection.

    This function iterates over the items in the passed dictionary containing linear equations, and for 
    each one calculates an intersection point with a vertical line at the provided ion mobility value. 
    The function returns the ion mobility value and the mz intersection.
    """

    for key, value in dict_equation.items():
        if key[2] <= IM and IM < key[3]:
            im_intersection = (IM - value[1]) / value[0]
        else:
            continue
        return IM, im_intersection
    

def create_basic_scheme_with_scans(
    scans,
    scan_mode,
    window_modification,
    no_of_combined_scans,
    window_pattern,
    window_overlap
):
    """Generates a basic scheme of isolation windows for the given scan mode and parameters without
    the scan area information.

    Parameters:
    scans (int): Number of synchro scans.
    scan_mode (str): There are three basic acquisition schemes for the scan mode synchro-PASEF: 
        'classical_synchro-PASEF', 'highly_accurate_synchro-PASEF' or 'individual_synchro-PASEF'.
    window_modification (str): Parameter for modifying the scan window:'None', 'window_overlap' 
        or 'staggered'.
    no_of_combined_scans (int): Number of combined scans for a staggered and classical-synchro-PASEF 
        scheme.
    window_overlap (float): Parameter for overlapping the scan window (0.2 means 20% of smallest window).
    window_pattern (list): Pattern of the windows for an individually defined synchro-PASEF method.

    Returns:
    dict: Dictionary with isolation window coordinates.

    This function creates a basic scheme of the scan isolation windows as per the provided scan mode 
    requirements and parameters. The scheme is a dictionary, in which the keys are the isolation window 
    numbers, and the values are the relevant isolation window coordinates.
    The function choose different isolation window scheme creations functions based on the provided 
    scan mode and window modification option.
    Assume that the following functions are defined elsewhere in the code: 'classic_synchropasef', 
    'classic_synchropasef_staggered', 'highly_accurate_synchropasef', 
    'highly_accurate_synchropasef_staggered', 'individual_synchropasef'.
    """
    dict_isolation_windows_temp = dict()
    if scan_mode == 'classical_synchro-PASEF':
        if window_modification == 'None':
            dict_isolation_windows = classic_synchropasef(
                dict_isolation_windows_temp,
                scans,
                window_overlap = 0.0
            )
        elif window_modification == 'overlap':
            dict_isolation_windows = classic_synchropasef(
                dict_isolation_windows_temp,
                scans,
                window_overlap
            )
        elif window_modification == 'staggered':
            dict_isolation_windows = classic_synchropasef_staggered(
                dict_isolation_windows_temp,
                scans,
                no_of_combined_scans,
            )
    elif scan_mode == 'highly_accurate_synchro-PASEF':
        if window_modification == 'None':
            dict_isolation_windows = highly_accurate_synchropasef(
                dict_isolation_windows_temp,
                scans,
                window_overlap = 0.0
            )
        elif window_modification == 'overlap':
            dict_isolation_windows = highly_accurate_synchropasef(
                dict_isolation_windows_temp,
                scans,
                window_overlap
            )
        elif window_modification == 'staggered':
            dict_isolation_windows = highly_accurate_synchropasef_staggered(
                dict_isolation_windows_temp,
                scans,
            )
    elif scan_mode == 'individual_synchro-PASEF':
        dict_isolation_windows = individual_synchropasef(
            dict_isolation_windows_temp,
            window_pattern,
        )
    else:
        print("Error: Scan mode not provided")

    return dict_isolation_windows


def window_calculator(
    df,
    scan_ratios,
    dict_isolation_windows,
    im_position, #bottom = 1, top = 2
):
    """
    Calculates the starting points (mz_start) and the widths (mz_width) of different synchro scans for 
    a synchro-PASEF method. This function is executed once to calculate the m/z positions and widths for 
    the lower ion mobility position of the scan area parallelogram and separatly for the upper one.

    Parameters:
    df: The input dataframe containing the upper_mz and lower_mz values for each im_position. These 
        coordinates define the scan area of the acquisition method.
        E.g.    lower_IM  lower_mz_1  upper_mz_1  upper_IM  lower_mz_2  upper_mz_2
             0       0.7  158.457425  458.457425       1.3  1047.20167  1347.20167
    scan_ratios: The list of scan ratios determines the window lengths of the synchro scans in the 
        acquisition scheme.
    dict_isolation_windows: A dictionary containing the normalized mz values and widths from the basic
        acquisition scheme. E.g. {'window0': {'mz': [0.0], 'width': [0]}, 'window1': {'mz': [1.0], ...]}}
    im_position (int): Defines the ion mobility position in the scan area polygone for that specific
        calculation. Can either be 1 (bottom) or 2 (top).

    Returns:
    mz_start (list): The starting points of the scan windows in the m/z dimension.
    mz_width (list): The widths of the scan windows in the m/z dimension.
    """
    total_length = round(df["upper_mz_"+im_position].iloc[0]-df["lower_mz_"+im_position].iloc[0], 4)
    mz_value = df["lower_mz_"+im_position].iloc[0]
    division_length = total_length/sum(scan_ratios)

    scan_ratio_pos = 0
    list_scan_ratios = list()
    for bin_size in scan_ratios:
        list_scan_ratios.append(scan_ratio_pos)
        scan_ratio_pos += division_length*bin_size

    mz_start = list()
    mz_width = list()
    for key in dict_isolation_windows.keys():
        position_in_list_scan_ratios = math.ceil(dict_isolation_windows[key]['mz'][0])
        window_overlap = (position_in_list_scan_ratios-dict_isolation_windows[key]['mz'][0])*division_length
        mz_start.append(mz_value + list_scan_ratios[position_in_list_scan_ratios] - window_overlap)
        mz_width_temp = list()
        for width in dict_isolation_windows[key]['width']:
            mz_width_temp.append(scan_ratios[width]*division_length)
        mz_width.append(sum(mz_width_temp))

    return mz_start, mz_width


def lin_equ(dict_equation, l1, l2):
    """
    Calculates the slope (m) and y-intercept (c) of a line defined by two points (l1, l2) and stores the 
    result in a dictionary. The line is encoded in the format l=(x,y).

    Parameters:
    dict_equation (dict): The dictonary is empty when parsed into this function.
    l1, l2 (tuple): Tuple that represents a single point on the 2D plane. Each tuple should contain two 
        floats, the x and y coordinates.

    Returns:
    dict_equation (dict): Dictionary that will be updated with the new line equation. The key is a tuple 
    of four numbers (x1, x2, y1, y2) representing two points on the line. The value is a tuple of two 
    numbers (m, c) representing the slope and y-intercept. 
    """
    m = float((l2[1] - l1[1])) / float(l2[0] - l1[0])
    c = (l2[1] - (m * l2[0]))

    dict_equation[(l1[0], l2[0], l1[1], l2[1])] = (m, c)

    return dict_equation


def classic_synchropasef(
    dict_isolation_windows_temp,
    scans,
    window_overlap
):
    """
    Updates a provided dictionary with normalized mz and width information for each synchro scan. It 
    generates a basic scheme of a  classical synchro-PASEF method, possibly with overlap.

    The dictionary is updated with the scan numbers as key, and each key is mapped to a sub-dictionary 
    with 'mz' and 'width' keys. For each scan, 'mz' is calculated as the scan number minus the product of 
    the window_overlap and a counter, while 'width' is simply the scan number.

    Parameters:
    dict_isolation_windows_temp (dict): Existing dictionary to be updated with the new window information. 
    scans (int): The total number of synchro scans to perform calculations for.
    window_overlap (int or float): The overlap of each window, used when calculating 'mz'.

    Returns:
    dict_isolation_windows_temp (dict): The original dictionary updated with the calculated window 
    information.
    """
    count = 0
    for i in range(scans):
        dict_isolation_windows_temp['window'+str(i)] = {'mz': [i-window_overlap*count], 'width': [i]}
        count += 1
    return dict_isolation_windows_temp


def classic_synchropasef_staggered(
    dict_isolation_windows_temp,
    scans,
    no_of_combined_scans,
):
    """
    Updates a provided dictionary with normalized mz and width information for each synchro scan. It 
    generates a basic scheme of a classical, but staggered synchro-PASEF method.

    The dictionary is updated with the scan numbers as key, and each key is mapped to a sub-dictionary 
    with 'mz' and 'width' keys. For each scan, 'mz' is simply the scan number, while 'width' is a range of 
    scan numbers starting from the current scan number and extending for the amount of 
    'no_of_combined_scans'. Or in other words, no_of_combined_scans indicates how many windows are combined
    to one.

    Parameters:
    dict_isolation_windows_temp (dict): Existing dictionary to be updated with the new window information. 
    scans (int): The total number of synchro scans to perform calculations for.
    no_of_combined_scans (int): The number of scans to combine into a single isolation window.

    Returns:
    dict_isolation_windows_temp (dict): The original dictionary updated with the calculated window 
    information.
    """
    for i in range(scans-no_of_combined_scans+1):
        dict_isolation_windows_temp['window'+str(i)] = {
            'mz': [i],
            'width': [i + window_pos for window_pos in range(no_of_combined_scans)]
        }
    return dict_isolation_windows_temp


def highly_accurate_synchropasef(
    dict_isolation_windows_temp,
    scans,
    window_overlap
):
    """
    Updates a provided dictionary with normalized mz and width information for each synchro scan. It 
    generates a basic scheme of a highly accurate synchro-PASEF method, possibly with overlap.

    The dictionary is updated with the scan numbers as key, and each key is mapped to a sub-dictionary 
    with 'mz' and 'width' keys. For each pair of scans, 'mz' for the first is simply 0, and  'mz' for 
    the second is calculated as the scan number plus 1 minus the window_overlap.
    'width' for both is a range of scan numbers starting from the current scan number and encompassing 
    all remaining scans.

    Parameters:
    dict_isolation_windows_temp (dict): Existing dictionary to be updated with the new window information. 
    scans (int): The total number of synchro scans to perform calculations for.
    window_overlap (int or float): The overlap of each window, used when calculating 'mz' for the 
    second window of each pair.

    Returns:
    dict_isolation_windows_temp (dict): The original dictionary updated with the calculated window 
    information.
    """
    count = 0
    for i in range(0, scans-1):
        dict_isolation_windows_temp['window'+str(count)] = {
            'mz': [0],
            'width': [window_pos for window_pos in range(i+1)]
        }
        count+=1
        dict_isolation_windows_temp['window'+str(count)] = {
            'mz': [i+1-window_overlap],
            'width': [i+1+window_pos for window_pos in range(scans-i-1)]
        }
        count+=1
    return dict_isolation_windows_temp


def highly_accurate_synchropasef_staggered(
    dict_isolation_windows_temp,
    scans,
):
    """
    Updates a provided dictionary with normalized mz and width information for each scan. It generates 
    a basic scheme of a staggered highly accurate synchro-PASEF method.

    The dictionary is updated with the scan numbers as key, and each key is mapped to a sub-dictionary 
    with 'mz' and 'width' keys. For each pair of scans, 'mz' for the first is simply 0, and  'mz' for 
    the second is calculated as the scan number plus 1.
    'width' for both is a range of scan numbers starting from the current scan number and encompassing 
    all remaining scans, while overlapping by one window width.

    Parameters:
    dict_isolation_windows_temp (dict): Existing dictionary to be updated with the new window information. 
    scans (int): The total number of synchro scans to perform calculations for.

    Returns:
    dict_isolation_windows_temp (dict): The original dictionary updated with the calculated window 
    information.
    """
    count = 0
    for i in range(0, scans-2):
        dict_isolation_windows_temp['window'+str(count)] = {
            'mz': [0],
            'width': [window_pos for window_pos in range(i+2)]
        }
        count+=1
        dict_isolation_windows_temp['window'+str(count)] = {
            'mz': [i+1],
            'width': [i+1+window_pos for window_pos in range(scans-i-1)]
        }
        count+=1
    return dict_isolation_windows_temp


def individual_synchropasef(
    dict_isolation_windows_temp,
    window_pattern,
):
    """
    Updates a provided dictionary with normalized mz and width information for each scan. The 
    function generates a basic scheme of a synchro-PASEF method based on an individual pattern of windows 
    provided.

    The dictionary is updated with the scan numbers as key, and each key is mapped to a 
    sub-dictionary with 'mz' and 'width' keys. For each window in the pattern, 'mz' is the first value 
    of the window range, and 'width' is a range from the first to the second value in the window range.

    Parameters:
    dict_isolation_windows_temp (dict): Existing dictionary to be updated with the new window 
    information. 
    window_pattern (list): A list of tuples, where each tuple represents a range of synchro scans.

    Returns:
    dict_isolation_windows_temp (dict): The original dictionary updated with the calculated window 
    information.
    """
    count = 0
    for i in window_pattern:
        dict_isolation_windows_temp['window'+str(count)] = {
            'mz': [i[0]],
            'width': [window_pos for window_pos in range(i[0], i[1]+1)]
        }
        count+=1

    return dict_isolation_windows_temp
