from diaid_pasef import loader
from diaid_pasef import method_creator
from diaid_pasef import method_evaluation
from diaid_pasef import graphs

# for data manipulation:
import pandas as pd

# for creating new folders:
import os

# todo: accelerate slow steps: loading of proteomics library, kernel density
# estimation, Bayesian optimization
# todo: change file locations to general solution


def create_folder(
    paths: list,
) -> None:
    """This function creates a folder structure where all output files are
        printed.

    Parameters:
    paths (list): a list of all path names that should be created.
    """
    for create_item in paths:
        if os.path.exists(create_item):
            next
        else:
            os.makedirs(create_item)


def library_information(method_conf, save_at):

    library = loader.load_library(
        method_conf["input"]["library_name"],
        method_conf["input"]["analysis_software"],
        method_conf["input"]["PTM"]
        )

    # 1st calculate kernel density coordinates:
    xi, yi, zi = graphs.kernel_density_calculation(
        library,
        method_conf["graphs"]["numbins"]
        )

    # plot important plots of the library to understand the method creation
    # step
    graphs.plot_precursor_distribution_as_histogram(
        library, method_conf["graphs"],
        save_at +
        '/input_library/Histogram_precursor_distribution_in_library.pdf'
        )
    graphs.plot_density(
        xi,
        yi,
        zi,
        method_conf["graphs"],
        save_at + '/input_library/Kernel_density_distribution_library.pdf')

    # save library specific information
    dict_charge_of_precursor = method_evaluation.calculate_percentage_multiple_charged_precursors(library)
    pd.DataFrame(
        {
            "column name": list(dict_charge_of_precursor.keys()),
            "column value": list(dict_charge_of_precursor.values())
        }
        ).to_csv(
            save_at
            + '/input_library/percentage_of_multiple_charged_precursors.csv',
            index=False
            )

    return xi, yi, zi, library


def precursor_within_scan_area(library, method_parameters, save_at):
    # save parameters for method evaluation as .csv
    dict_precursors_within_scan_area = method_evaluation.calculate_precursor_within_scan_area(
        library,
        method_parameters["mz"],
        method_parameters["ion_mobility"]
    )
    print(dict_precursors_within_scan_area)

    pd.DataFrame({
        "column name": list(dict_precursors_within_scan_area.keys()),
        "column value": list(dict_precursors_within_scan_area.values())
    }).to_csv(
        save_at + '/final_method/precursors_within_scan_area.csv',
        index=False)


def create_final_method(
    library: pd.DataFrame,
    method_parameters: dict,
    dim: list,
    method_conf: dict,
) -> pd.DataFrame:
    """This function calculates the final diaPASEF window scheme with the
        optimized scan area coordinates and writes the .txt diaPASEF parameter
        file as input file for timsControl. By default, the scan area
        coordinates are up-shifted by 0.022 1/K0 since the quadrupole has a
        dead time at the beginning of the diaPASEF window where no signal is
        acquired. The up-shift should prevent that this dead area lies within
        an area with a high precursor density.

    Parameters:
    library (pd.DataFrame): a pre-filtered data frame with unified column names
        containing all required precursor information.
    method_parameters (dict): dictionary, which includes all input parameters
        for creating a diaPASEF method (e.g., m/z-range, ion mobility-range,
        number of diaPASEF scans, number of ion mobility steps, the overlap of
        the diaPASEF windows, a shift of the final window scheme in the ion
        mobility axis to account for the blind spot of the changing quadrupole
        at the top of each diaPASEF window).
    dim (list): y-coordinates of the scan area = A1, A2, B1, B2. x-coordinate
        for A1 and B1 is the lower m/z-range value, x-coordinate for A2 and B2
        is the upper m/z-range value.
    method_conf (dict): this dictionary contains all input parameters for all
        sub-functions.

    Returns:
    pd.DataFrame: data frame that contains the scan type (PASEF), scan number,
        and the corresponding diaPASEF window coordinates for each window per
        scan.
    """
    file_name_method = method_conf["input"]["save_at"] + '/final_method/diaPASEF_method.txt'
    # create the diaPASEF method scheme; dim[0],dim[1],dim[2],dim[3] == A1, A2,
    # B1, B2.
    df_parameters_final = method_creator.method_creation(
        library["mz"],
        method_parameters,
        dim[0] + method_parameters["shift_of_final_method"],
        dim[1] + method_parameters["shift_of_final_method"],
        dim[2] + method_parameters["shift_of_final_method"],
        dim[3] + method_parameters["shift_of_final_method"]
    )
    method_creator.create_parameter_dataframe(
        df_parameters_final,
        file_name_method
    )
    return df_parameters_final


def final_method_information(
    df_parameters_final,
    xi,
    yi,
    zi,
    method_conf,
    save_at,
    library,
    method_parameters,
    dim
):
    # plot created method on top of kernel density estimation of library
    graphs.plot_density_and_method(
        df_parameters_final,
        xi,
        yi,
        zi,
        method_conf["graphs"],
        save_at +'/final_method/Kernel_density_distribution_and_final_method.pdf'
    )

    # save parameters for method evaluation as .csv
    dict_evaluation_of_final_method = method_evaluation.coverage(
        df_parameters_final,
        library
    )

    if dim is None:
        next
    else:
        dict_evaluation_of_final_method["final A1, A2, B1, B2 values"] = str([
            dim[0] + method_parameters["shift_of_final_method"],
            dim[1] + method_parameters["shift_of_final_method"],
            dim[2] + method_parameters["shift_of_final_method"],
            dim[3] + method_parameters["shift_of_final_method"]
            ]
        )

    pd.DataFrame({
        "column name": list(dict_evaluation_of_final_method.keys()),
        "column value": list(dict_evaluation_of_final_method.values())
    }).to_csv(
        save_at + '/final_method/parameters_to_evaluate_method.csv',
        index=False)


def evaluate_for_multiple_charged_prec(
    method_conf,
    save_at,
    library,
    df_parameters_final
):

    for charge in range(2, 5, 1):
        xi_charge, yi_charge, zi_charge = graphs.kernel_density_calculation_multiple_charge(
            library,
            method_conf["graphs"]["numbins"],
            charge)
        graphs.plot_density(
            xi_charge,
            yi_charge,
            zi_charge,
            method_conf["graphs"],
            save_at + '/evaluation_multiple_charged_prec/Kernel_density_distribution_library_'+ str(charge) +'x_charged_precursors.pdf')
        graphs.plot_density_and_method(
            df_parameters_final,
            xi_charge,
            yi_charge,
            zi_charge,
            method_conf["graphs"],
            save_at +'/evaluation_multiple_charged_prec/Kernel_density_distribution_and_final_method_'+ str(charge) +'x_charged_precursors.pdf'
        )
