from pydiaid.synchropasef import plots

import pandas as pd
import numpy as np
# import statistics

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


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


def coverage_calculated_at_center(library, scan_area_definition):
    """
    Calculates the fraction of unique peptides covered at the ion mobility apex in a given proteomics 
    library.

    It first determines if a peptide is within the defined scan area (calculation based on the ion 
    mobility apex), filters the resulting dataframe to include only those within the scan area, and then 
    entries. The function finally returns the fraction of unique peptide entries that are covered within 
    the scan area.

    Parameters:
    library (pd.DataFrame): The protein library data frame to calculate coverage for. It should have at 
        least the columns 'Peptide' and 'Charge'.
        
    scan_area_definition (list or tuple): A four-element list or tuple defining the boundaries of the scan 
        area.

    Returns:
    float: The fractional coverage of unique peptides in the library that are within the defined scan area. The 
           calculations are based on the ion mobility apex of the peptides.
    """
    lib = coverage(library, scan_area_definition[:2], scan_area_definition[2:])
    library_coverage_center = lib[lib["insideIM"]==True]
    library_coverage_center_wo_2x = library_coverage_center.drop_duplicates(['Peptide', 'Charge'])
    return len(library_coverage_center_wo_2x)/len(library.drop_duplicates(['Peptide', 'Charge']))


def calculate_coverage_total_per_scan_per_charge_state(
    df_parameters_final,
    library
):
    """
    Calculates the coverage of unique peptides for all scans and per scan for distinct peptide charge 
    states of a given proteomics library. The calculations are based on the ion mobility apex of the 
    peptides.

    It gathers necessary scan area parameters from the acquisition method definitions, determines if a 
    peptide is within the computed scan area, filters the resulting dataframe to include only covered 
    peptides, and finally returns a summary of covered peptides (including different charge states) for 
    all scans and per scan in percentage terms.

    Parameters:
    df_parameters_final (pd.DataFrame): Acquisition method parameters defining the scan areas.
    library (pd.DataFrame): Protein library to calculate coverage for. Columns should include 
        'Proteins', 'Peptide' and 'Charge'.

    Returns:
    pd.DataFrame: A summary table covering the calculate coverage for all precursors and different 
        charge states per scan.
    pd.DataFrame: A dataframe holding the method parameters including the scan area definitions.
    """
    df_temp = df_parameters_final[
        ["mobility pos.1 [1/K0]",
        "mass pos.1 start [m/z]",
        "mass pos.1 end [m/z]",
        "mobility pos.2 [1/K0]",
        "mass pos.2 start [m/z]"]
    ][df_parameters_final["type"]=="vista"].astype(float).reset_index(drop=True)
    
    # calculate coverage for all scans
    scan_area_definition_temp = [
        (min(df_temp["mass pos.1 start [m/z]"]), df_temp["mobility pos.1 [1/K0]"].iloc[0]), 
        (min(df_temp["mass pos.2 start [m/z]"]), df_temp["mobility pos.2 [1/K0]"].iloc[0]),
        (max(df_temp["mass pos.1 end [m/z]"]), df_temp["mobility pos.1 [1/K0]"].iloc[0]), 
        (max(df_temp["mass pos.2 start [m/z]"]+(df_temp["mass pos.1 end [m/z]"]-df_temp["mass pos.1 start [m/z]"])), df_temp["mobility pos.2 [1/K0]"].iloc[0])
    ]
    lib = coverage(library, scan_area_definition_temp[:2], scan_area_definition_temp[2:])
    library_coverage_center = lib[lib["insideIM"]==True]
    
    dict_lib_precursors = dict()
    dict_lib_precursors["all synchro scans"] = library_coverage_center
    
    # calculate coverage per scan
    dict_lib_precursors_per_scan = coverage_per_scan(
        df_temp,
        library,
    )
    dict_lib_precursors.update(dict_lib_precursors_per_scan)
    
    dict_evaluation = {
        "scan":[],
        "unique proteins in the library": len(library.drop_duplicates(['Proteins'])),
        "unique precursors in the library": len(library.drop_duplicates(['Peptide', 'Charge'])),
        "No. of covered proteins":[],
        "No. of covered precursors":[],
        "No. of covered, singly charged precursors": [],
        "No. of covered, doubly charged precursors": [],
        "No. of covered, triply charged precursors": [],
        "No. of covered, quadruply charged precursors": [],
        "proteins covered [%]": [],
        "precursors covered [%] (counting IM peak at center)": [],
        "singly charged precursors covered [%]": [],
        "doubly charged precursors covered [%]": [],
        "triply charged precursors covered [%]": [],
        "quadruply charged precursors covered [%]": [],  
    }
    
    for item in dict_lib_precursors.keys():
        library_per_scan_temp = dict_lib_precursors[item]
        dict_evaluation["scan"].append(item)
    
        filter_for = ["Proteins"]
        filtered_temp = len(library_per_scan_temp.drop_duplicates(filter_for))
        dict_evaluation["No. of covered proteins"].append(filtered_temp)
        dict_evaluation["proteins covered [%]"].append(round(filtered_temp/len(library.drop_duplicates(filter_for))*100, 1))
    
        filter_for = ["Peptide", "Charge"]
        filtered_temp = len(library_per_scan_temp.drop_duplicates(filter_for))
        dict_evaluation["No. of covered precursors"].append(filtered_temp)
        dict_evaluation["precursors covered [%] (counting IM peak at center)"].append(round(filtered_temp/len(library.drop_duplicates(filter_for))*100, 1))
    
        for key, value in {"singly":1, "doubly":2, "triply":3, "quadruply":4}.items():
            try:
                filtered_temp = len(library_per_scan_temp[library_per_scan_temp["Charge"]==value].drop_duplicates(filter_for))
                dict_evaluation["No. of covered, "+key+" charged precursors"].append(filtered_temp)
                dict_evaluation[key+" charged precursors covered [%]"].append(round(filtered_temp/len(library[library["Charge"]==value].drop_duplicates(filter_for))*100, 1)) 
            except Exception as e:
                dict_evaluation[key+" charged precursors covered [%]"].append(0)
                print(e, key+" charged precursours are not present") 
                
    df_coverage = pd.DataFrame(dict_evaluation).T
    df_coverage.columns = df_coverage.iloc[0].values
    return df_coverage[1:], df_temp


def calculate_slicing_and_coverage_in_total(library, df_temp):
    """
    Calculates both slicing coverage and total coverage for all scans and per scan for distinct peptide 
    charge states of a given proteomics library.

    The function requires a dataframe representing a library of proteins and the acquisition method 
    parameters. It calculates the slicing and overall coverage for all scans first, then calculates the 
    same per scan and aggregats the results.

    Parameters:
    library (pd.DataFrame): Protein library to calculate coverage for. Columns should include 'Proteins', 
    'Peptide' and 'Charge'.
    df_temp (pd.DataFrame): Acquisition method parameters defining the scan areas.

    Returns:
    pd.DataFrame: A dataframe presenting the slicing coverage for all precursors and different charge 
    states across all scans and per scan. The rows of the data frame represent different charge states 
    and peptides, while the columns present the coverage information for all synchro scans and per synchro
    scan.
    """
    # calculate slicing
    number_of_scans = len(df_temp)
    scans = range(1, number_of_scans+1)

    list_dict_slicing = list()
    dict_slicing = calculate_slicing_and_coverage_per_scan(library, df_temp, scans)
    list_dict_slicing.append(dict_slicing)

    for scan in scans:
        dict_slicing_temp = calculate_slicing_and_coverage_per_scan(library, df_temp, [scan])
        list_dict_slicing.append(dict_slicing_temp)

    df_slicing = pd.DataFrame(list_dict_slicing, index = ["all synchro scans"]+['#synchro scan '+ str(scan) for scan in range(1, number_of_scans+1)]).T

    return df_slicing


def coverage(library_temp, top_limit_IM, bottom_limit_IM, im = "IM"):
    """
    Calculates the coverage of peptides using the provided scan area limits. This function can calculate 
    the coverage based on the ion mobility apex but also every other ion mobility peak position.
   
    The function takes in a dataframe representing a library of proteins, coordinates representing the 
    scan area are within the parameters 'top_limit_IM' and 'bottom_limit_IM' and a default parameter 
    indicating the ion peak position ('IM' for apex by default). It utilizes these inputs to form a 
    polygon using the provided coordinates/limits. It checks for every peptide in the library if it 
    resides within this polygon region and adds a new column to the library dataframe indicating this 
    information.
   
    Parameters:
    library_temp (pd.DataFrame): Protein library to calculate coverage for. Columns should include 
        'Proteins', 'Peptide' and 'Charge'.
    top_limit_IM (tuple): The coordinates defining the top limit of the scan area.
    bottom_limit_IM (tuple): The coordinates defining the bottom limit of the scan area.
    im (str, optional): The ion mobility peak position to consider for the calculations, 'IM' for apex
        by default.
   
    Returns:
    pd.DataFrame: The originally provided library dataframe with an additional column "inside"+im 
    indicating whether a peptide resides within the polygon region or not.
    """
    polygon = Polygon([top_limit_IM[0], bottom_limit_IM[0], bottom_limit_IM[1], top_limit_IM[1]])
    library_temp["inside"+im] = library_temp.apply(lambda x: polygon.contains(Point(x["mz"], x[im])), axis = 1)

    return library_temp


def coverage_per_scan(
    df_temp,
    library_temp,
):
    """
    Calculates the coverage of peptides for each scan within the defined scan areas.

    The function loops through the provided dataframe containing acquisition method parameters
    defining scan areas (df_temp) and defines the scan area based on the parameters. It then 
    calculates the coverage of peptides within this area by calling the 'coverage' function.
    The peptides that reside within the defined area are then added to a dictionary with 
    scan number as the key. 

    Parameters:
    df_temp (pd.DataFrame): Acquisition method parameters defining the scan areas.
    library_temp (pd.DataFrame): Protein library to calculate coverage for. Columns should include 
        'Proteins', 'Peptide', and 'Charge'.

    Returns:
    dict: A dictionary where keys are scan numbers and values are dataframes of peptides present
    within the defined scan area for that particular scan.
    """
    dict_coverage_per_scan = dict()
    for index, row in df_temp.iterrows():
        scan_area_definition_temp = [
            (row["mass pos.1 start [m/z]"], row["mobility pos.1 [1/K0]"]), 
            (row["mass pos.2 start [m/z]"], row["mobility pos.2 [1/K0]"]),
            (row["mass pos.1 end [m/z]"], row["mobility pos.1 [1/K0]"]), 
            (row["mass pos.2 start [m/z]"]+(row["mass pos.1 end [m/z]"]-row["mass pos.1 start [m/z]"]), row["mobility pos.2 [1/K0]"])
        ]
        lib = coverage(library_temp, scan_area_definition_temp[:2], scan_area_definition_temp[2:])
        library_coverage_center = lib[lib["insideIM"]==True]
        dict_coverage_per_scan["#synchro scan "+str(index+1)] = library_coverage_center
    return dict_coverage_per_scan


def calculate_slicing_and_coverage_per_scan(library, df_temp, scans):
    """
    Calculates slicing and coverage percentages and a combination theirof for selected scans of a given 
    proteomics library. The calculations require the ion mobility apex and ion mobility peak width of 
    peptides.

    The function loops through the given scans. First it retrieves the linear fuctions that frame the scan 
    area for each scan (get_scan_lines), calculates if a peak was partitioned/sliced 
    (finding_mobility_edges), if the peptide is covered considering the complete ion mobility peak 
    (precursor_covered_with_width) and calculates the fragment coverage alias the proportion of the 
    covered ion mobility peak width. The results are aggregated into a dictionary which is returned.

    Parameters:
    library (pd.DataFrame): Protein library to calculate coverage for. Columns should include 'Proteins', 
    'Peptide', 'Charge', 'mz', and 'IM'.
    df_temp (pd.DataFrame): Acquisition method parameters defining the scan areas.
    scans (list): List of scan numbers on which to perform the slicing coverage calculation.

    Returns:
    dict: A dictionary holding the slicing and coverage results for selected scans. The dictionary 
    includes percentages for sliced, covered, and covered & sliced precursors for each scan, as well as 
    the average fragment coverage.
    """ 
    dict_slicing = dict()
    library_temp = library.copy()

    # get scan area equations scan specific 
    mobility_equation_dict = plots.get_scan_lines(df_temp)
    filtered_mobility_equation_dict = {}
    for filter_item in scans:
        filter_string = str(filter_item)
        filtered_dict_temp = {k:v for (k,v) in mobility_equation_dict.items() if filter_string in k}
        filtered_mobility_equation_dict.update(filtered_dict_temp)

    library_temp["mobility_edges"] = library_temp.apply(
        lambda x: find_mobility_edges(
            mz=x["mz"],
            im=x["IM"],
            im_width=x["IMlength"],
            mobility_equation_dict=filtered_mobility_equation_dict,
        ),
        axis=1,
    )

    library_temp["number_of_mobility_edges"] = library_temp.mobility_edges.apply(len)
    library_temp["scans"] = library_temp.mobility_edges.apply(lambda x: np.unique([int(i[4]) for i in x.keys()]))

    percentage_sliced_precursors_fixed_width = (library_temp.number_of_mobility_edges > 0).sum() / len(library_temp)

    library_temp['IM_top'] = library_temp['IM'] + (library_temp["IMlength"])/2
    library_temp['IM_bottom'] = library_temp['IM'] - (library_temp["IMlength"])/2
    library_temp['covered'] = library_temp.apply(lambda x: precursor_covered_with_width(x['mz'], x['IM_top'], x['IM_bottom'], mobility_equation_dict, scans), axis = 1)

    percise_precursor_coverage_fixed_width = len(library_temp[library_temp['covered'] == 'yes']) / len(library_temp)
    try:
        precursor_sliced_and_covered_fixed_width = len(library_temp[(library_temp['covered'] == 'yes') & (library_temp['number_of_mobility_edges'] > 0)]) / len(library_temp[library_temp['covered'] == 'yes'])
        dict_slicing["covered and sliced precursors [%] (counting complete IM peak width)"] = round(100 * precursor_sliced_and_covered_fixed_width, 2)
    except Exception as e:
        precursor_sliced_and_covered_fixed_width = 0
        print(e, "Precursours that are covered and sliced are not present") 
    library_temp['fragment_coverage'] = library_temp.apply(lambda x: calculate_fragment_coverage(x, scans), axis=1)

    dict_slicing["average fragment coverage [%] (counting complete IM peak width)"] = round(library_temp['fragment_coverage'].mean(), 4)*100
    dict_slicing["sliced precursors [%] (counting complete IM peak width)"] = round(100 * percentage_sliced_precursors_fixed_width, 2)
    dict_slicing["covered precursors [%] (counting complete IM peak width)"] = round(100 * percise_precursor_coverage_fixed_width, 2)
    return dict_slicing


def find_mobility_edges(
    mz,
    im,
    im_width,
    mobility_equation_dict
):
    """
    Finds if a peptide was sliced/partitioned along the ion mobility width (= mobility edges). The 
    calculation is based on linear functions that represent the scan area limits and the mz, ion 
    mobility (im) position and ion mobility width (im_width). The parameters of the linear function such
    as slope and intercept are defined in mobility_equation_dict.

    This function takes the m/z (mass to charge ratio) value for each precursor, applies the scan area
    equations (provided as a dictionary, mobility_equation_dict) for each precursor and calculates the 
    slicing position based on the product of m/z value and slope, added to the intercept.

    If the slicing position cuts the precursor based on the given peptide information (ion mobility and 
    ion mobility width, the slicing position value is added to the dictionary with the corresponding 
    precursor key.

    Parameters:
    mz (float): The m/z value for a specific precursor.
    im (float): The ion mobility value for a specific precursor.
    im_width (float): The ion mobility width for a specific precursor that is used to determine if a 
        precursor is sliced.
    mobility_equation_dict (dict): The dictionary attempting to find the slicing positions (mobility edges) 
        for different precursors by using the linear functions of the scan area in the form of slope and 
        intercept values.

    Returns:
    dict_precursor_slicing (dict): A dictionary that holds precursor keys and corresponding slicing 
    positions (where precursor is cut in the ion mobility dimension).
    """
    dict_precursor_slicing = dict()
    for key, (slope, intercept) in mobility_equation_dict.items():
        slicing_position = mz * slope + intercept
        if slicing_position_cuts_precursor_in_im(
            slicing_position,
            im,
            im_width,
        ):
            dict_precursor_slicing[key] = slicing_position
    return dict_precursor_slicing


def precursor_covered_with_width(mz, im_top, im_bottom, mobility_equation_dict, scans):
    """
    Determines whether a precursor is covered within the scan area using the complete ion mobility peak
    width of each peptide. The scan area is defined by linear equations that define the top and 
    bottom limits of the scans parallelogram.

    This function uses the given m/z (mass to charge ratio) value for the precursor, the ion mobility 
    limits of each peak (im_top and im_bottom), as well as the linear function parameters for different 
    scans stored in a dictionary (mobility_equation_dict), to determine whether the precursor is within 
    these limits.

    The limits for ion mobility are determined by applying the linear function equations to the m/z 
    value. The function then checks if the top and bottom ion mobility part of the precursor are within 
    these limits.

    If they are, then the precursor is considered as covered and the function returns 'yes'. If they 
    are not, it returns 'no'.

    Parameters:
    mz (float): The m/z value for a specific precursor.
    im_top (float): The upper limit of the ion mobility peak width for the precursor.
    im_bottom (float): The lower limit of the ion mobility peak width for the precursor.
    mobility_equation_dict (dict): A dictionary holding the linear equations for each scan.
    scans (iterable): An iterable indicating the scans that the function should consider.

    Returns:
    str: Whether the precursor is covered ('yes' or 'no').
    """
    top_border = mobility_equation_dict['scan'+str(min(scans))+'top'][0] * mz + mobility_equation_dict['scan'+str(min(scans))+'top'][1]
    bottom_border = mobility_equation_dict['scan'+str(max(scans))+'bottom'][0] * mz + mobility_equation_dict['scan'+str(max(scans))+'bottom'][1]
    if top_border > im_top > bottom_border and  top_border > im_bottom > bottom_border:
        return 'yes'
    elif top_border > im_top > bottom_border or  top_border > im_bottom > bottom_border:
        return 'yes'
    elif im_top > top_border > im_bottom or  im_top > bottom_border > im_bottom:
        return 'yes'
    else:
        return 'no'


def calculate_fragment_coverage(
    item,
    scans,
):
    """"
    Calculates the fragment coverage for precursors.

    This function computes the fragment coverage based on slicing position (= mobility edges), which 
    should be already calculated and stored in the item dataframe passed to the function. It does so 
    by checking if there are any mobility edges for a given item, and if so, it calculates the amount 
    of coverage.

    The fragment coverage is calculated per scan considering the difference between the top and bottom 
    ion mobility peak position for each precursor. It uses these to compute a factor of 'not covered' 
    value which is then subtracted from 1 to get the 'covered' value.

    If there are no slicing positions for the precursor, this function checks whether the precursor is 
    covered. If it is, it returns 1, indicating full coverage.

    Parameters:
    item (dict): A dataframe containing details about the fragment to calculate coverage for. Must 
        include the columns 'mobility_edges', 'IM_top', 'IM_bottom' and 'covered'.
    scans (iterable): An iterable indicating the scans that the function should consider.

    Returns:
    float: The calculated coverage of the fragment, represented as a decimal percentage (1 equals to 
    100% and 0 equals to 0%).
    """
    item_edges = item['mobility_edges']
    if len(item_edges)>0:
        item_list_keys = list(item_edges.keys())
        item_list_values = list(item_edges.values())
        not_covered = 0
        if str(scans[0])+"top" in item_list_keys[0]:
            not_covered += abs(item_list_values[0] - item["IM_top"])/(item["IM_top"]-item["IM_bottom"])
        elif str(scans[-1])+"bottom" in item_list_keys[-1]:
            not_covered += abs(item_list_values[-1] - item["IM_bottom"])/(item["IM_top"]-item["IM_bottom"])
        else:
            not_covered += 0
        return 1-not_covered
    elif len(item_edges) == 0:
        if item["covered"] == "yes":
            return 1
        else:
            next


def slicing_position_cuts_precursor_in_im(
    slicing_position,
    im,
    im_width,
):
    """
    Determines if a slicing position intersects with the provided ion mobility peak width/range.

    Parameters:
    slicing_position (float): The calculated slicing position.
    im (float): The given ion mobility peak posiiton.
    im_width (float): The ion mobility peak width.

    Returns:
    bool: True if the slicing position cuts the precursor in the ion mobility width, False otherwise.
    """
    if slicing_position <= im - im_width / 2:
        return False
    if slicing_position >= im + im_width / 2:
        return False
    return True

