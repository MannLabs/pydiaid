from diaid_pasef import loader
from diaid_pasef import method_creator
from diaid_pasef import method_evaluation
from diaid_pasef import graphs

# for data manipulation:
import pandas as pd

# for creating new folders:
import os

# importing for scientific and numeric manipulations:
import numpy as np
import skopt

# todo: accelerate slow steps: loading of proteomics library, kernel density
# estimation, Bayesian optimization
# todo: function for creating specific method: specific A1, A2, B1, B2
# todo: function for evaluation of own method: plot method on to of mz/IM
# plain, calculate coverage

# todo: make this a parameter
# DefaultDiaParameters  = [....]
def library_plus_evaluate_method(
    method_conf: dict
) -> None:
    method_parameters = method_conf["method_parameters"]
    save_at = method_conf["input"]["save_at"]

    folder_paths = [
        save_at,
        save_at+'/input_library',
        save_at+'/final_method'
        ]
    create_folder(folder_paths)

    xi, yi, zi, library = library_information(method_conf, save_at)

    df_parameters_final = pd.read_csv(
        method_conf["input"]["diaPASEF_method_only_used_for_method_evaluation"],
        skiprows=4,
        names=["MS Type", "Cycle Id", "Start IM", "End IM", "Start Mass", "End Mass", "CE"]
    )

    # save input_parameter as .csv
    pd.DataFrame(
        {
        "column name": list(method_conf["input"].keys()) +
            list(method_conf["method_parameters"].keys()) +
            list(method_conf["graphs"].keys()),
        "column value": list(method_conf["input"].values()) +
            list(method_conf["method_parameters"].values()) +
            list(method_conf["graphs"].values())
        }
    ).to_csv(save_at + '/final_method/input_parameters.csv', index=False)

    final_method_information(df_parameters_final, xi, yi, zi, method_conf, save_at, library, method_parameters, None)

    print("Method evaluated")


def library_plus_create_methods(
    method_conf: dict
) -> None:
    method_parameters = method_conf["method_parameters"]
    save_at = method_conf["input"]["save_at"]
    dim = method_parameters["scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF"]

    folder_paths = [
        save_at,
        save_at+'/input_library',
        save_at+'/final_method'
        ]
    create_folder(folder_paths)

    xi, yi, zi, library = library_information(method_conf, save_at)

    df_parameters_final = create_final_method(
        library,
        method_parameters,
        dim,
        method_conf
        )

    # save input_parameter as .csv
    pd.DataFrame(
        {
        "column name": list(method_conf["input"].keys()) +
            list(method_conf["method_parameters"].keys()),
        "column value": list(method_conf["input"].values()) +
            list(method_conf["method_parameters"].values())
        }
    ).to_csv(save_at + '/final_method/input_parameters.csv', index=False)

    print("Method created")


def run_all(
    method_conf: dict,
) -> None:
    """This function carries out all sub-functions step by step: loading of the
    proteomics library, creating a folder for the output information,
    calculation of the kernel density estimation for the density plots,
    Bayesian optimizationof diaPASEF method parameters by trying multiple scan
    area coordinates using a Gaussian process, writing all input and output
    information in csv files, plot with diaPASEF windows overlaid on a kernel
    density estimation of precursors.

    Parameters:
    method_conf (dict): this dictionary contains all input parameters for all
        sub-functions.

    """
    method_parameters = method_conf["method_parameters"]
    optimizer_parameters = method_conf["optimizer"]
    save_at = method_conf["input"]["save_at"]

    folder_paths = [
        save_at, save_at+'/optimization_plots',
        save_at+'/input_library',
        save_at+'/final_method'
        ]
    create_folder(folder_paths)

    xi, yi, zi, library = library_information(method_conf, save_at)

    opt_result = optimization(library, method_parameters, xi, yi, zi, method_conf, optimizer_parameters)

    df_parameters_final = create_final_method(
        library,
        method_parameters,
        opt_result.x,
        method_conf
        )

    # save input_parameter as .csv
    pd.DataFrame(
        {
        "column name": list(method_conf["input"].keys()) +
            list(method_conf["method_parameters"].keys()) +
            list(method_conf["graphs"].keys()) +
            list(method_conf["optimizer"].keys()),
        "column value": list(method_conf["input"].values()) +
            list(method_conf["method_parameters"].values()) +
            list(method_conf["graphs"].values()) +
            list(method_conf["optimizer"].values())
        }
    ).to_csv(save_at + '/final_method/input_parameters.csv', index=False)

    final_method_information(df_parameters_final, xi, yi, zi, method_conf, save_at, library, method_parameters, opt_result)


def final_method_information(df_parameters_final, xi, yi, zi, method_conf, save_at, library, method_parameters, opt_result):
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
    dict_precursors_within_mz = method_evaluation.calculate_precursor_within_mz_range(
        library,
        method_parameters["mz"]
    )
    dict_precursors_coverage = method_evaluation.coverage(
        df_parameters_final,
        library
    )
    dict_evaluation_of_final_method = {
        **dict_precursors_within_mz,
        **dict_precursors_coverage
    }

    if opt_result is None:
        next
    else:
        dict_evaluation_of_final_method["final A1, A2, B1, B2 values"] = str([
            opt_result.x[0] + method_parameters["shift_of_final_method"],
            opt_result.x[0] + opt_result.x[1] + method_parameters["shift_of_final_method"],
            opt_result.x[0] + opt_result.x[2] + method_parameters["shift_of_final_method"],
            opt_result.x[0] + opt_result.x[3] + method_parameters["shift_of_final_method"]]
        )
    pd.DataFrame({
        "column name": list(dict_evaluation_of_final_method.keys()),
        "column value": list(dict_evaluation_of_final_method.values())
    }).to_csv(
        save_at + '/final_method/parameters_to_evaluate_method.csv',
        index=False)


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


def optimization(library, method_parameters, xi, yi, zi, method_conf, optimizer_parameters):
    # Optimize method parameters Bayesian optimization using Gaussian process,
    # values should follow a multivariate gaussian curve
    opt_result = skopt.gp_minimize(
        lambda dim: single_optimization_run(
            library,
            method_parameters,
            optimizer_parameters["evaluation_parameter"],
            dim,
            xi,
            yi,
            zi,
            method_conf
        ), [
            (optimizer_parameters["YA1"][0], optimizer_parameters["YA1"][1]),
            (optimizer_parameters["YA2"][0]-optimizer_parameters["YA1"][0], optimizer_parameters["YA2"][1]-optimizer_parameters["YA1"][0]),
            (optimizer_parameters["YB1"][0]-optimizer_parameters["YA1"][0], optimizer_parameters["YB1"][1]-optimizer_parameters["YA1"][0]),
            (optimizer_parameters["YB2"][0]-optimizer_parameters["YA1"][0], optimizer_parameters["YB2"][1]-optimizer_parameters["YA1"][0])
        ],
        optimizer_parameters["n_random_starts"],
        optimizer_parameters["n_calls"],  # n to test
        optimizer_parameters["n_start"],  # random starts
        optimizer_parameters["initial_points"],  # initial points
        )

    print("########")
    print("BEST RESULT")
    print("INPUT: " + str([opt_result.x[0], opt_result.x[0]+opt_result.x[1], opt_result.x[0]+opt_result.x[2], opt_result.x[0]+opt_result.x[3]]))
    print("OUTPUT: " + str(1.0 / opt_result.fun))
    print("########")

    write_tested_parameters_into_csv(
        opt_result,
        optimizer_parameters,
        method_conf
        )

    return opt_result


def single_optimization_run(
    library: pd.DataFrame,
    method_parameters: dict,
    evaluation_parameter: str,
    dim: list,
    xi: np.ndarray,
    yi: np.ndarray,
    zi: np.ndarray,
    method_conf: dict,
) -> int:
    """This function contains one method creation step. First a diaPASEF scheme
        is calculated based on the scan area coordinates, second the evaluation
        parameter is calculated to estimate the optimization potential, last
        the diaPASEF method scheme is printed on top of the kernel density
        estimation of the used proteomics library to follow the optimization
        process visually.

    Parameters:
    library (pd.DataFrame): a pre-filtered data frame with unified column names
    containing all required precursor
        information.
    method_parameters (dict): dictionary, which includes all input parameters
        for creating a diaPASEF method (e.g., m/z-range, ion mobility-range,
        number of diaPASEF scans, number of ion mobility steps, the overlap of
        the diaPASEF windows, a shift of the final window scheme in the ion
        mobility axis to account for the blind spot of the evaluation_parameter
        (str): parameter, which was selected for optimization (e.g., all
        precursors, all doubly charged precursors).
    dim (list): y-coordinates of the scan area = A1, A2, B1, B2. x-coordinate
        for A1 and B1 is the lower m/z-range value, x-coordinate for A2 and B2
        is the upper m/z-range value.
    xi, yi, zi (numpy.ndarray): coordinates of the kernel density estimation
        where zi indicates the density.
    method_conf (dict): this dictionary contains all input parameters for all
        sub-functions.

    Returns:
    int: the value of the parameter, that was selected for optimization (e.g.,
        all precursors, all doubly charged precursors).
    """
    # create the diaPASEF method scheme; dim[0],dim[1],dim[2],dim[3] == A1, A2,
    # B1, B2. Calculate the optimization parameter.
    df_parameters_final = method_creator.method_creation(
        library["mz"],
        method_parameters,
        dim[0],
        dim[0]+dim[1],
        dim[0]+dim[2],
        dim[0]+dim[3]
    )
    dict_precursors_coverage = method_evaluation.coverage(
        df_parameters_final,
        library
    )
    result = 1.0 / float(dict_precursors_coverage[evaluation_parameter])
    print("RUN WITH: " + str([dim[0], dim[0]+dim[1], dim[0]+dim[2], dim[0]+dim[3]]) + " | RESULT: " + str(1.0 / result))

    # plot the created diaPASEF methods on top of the kernel density estimation
    file_name_optimization_plot = method_conf["input"]["save_at"] + '/optimization_plots/Optimization_plot_A1_' + str(dim[0]) + "_A2_ " +str(dim[0]+dim[1]) + "_B1_ " + str(dim[0]+dim[2]) + "_B2_ " + str(dim[0]+dim[3]) + "_result_" + str(1.0 / result) + ".png"
    graphs.plot_density_and_method(
        df_parameters_final,
        xi,
        yi,
        zi,
        method_conf["graphs"],
        file_name_optimization_plot
    )
    return result


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
        dim[0] + dim[1] + method_parameters["shift_of_final_method"],
        dim[0] + dim[2] + method_parameters["shift_of_final_method"],
        dim[0] + dim[3] + method_parameters["shift_of_final_method"]
    )
    method_creator.create_parameter_dataframe(
        df_parameters_final,
        file_name_method
    )
    return df_parameters_final


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


def write_tested_parameters_into_csv(
    opt_result: pd.DataFrame,  # scipy.optimize.optimize.OptimizeResult
    optimizer_parameters: str,
    method_conf: dict,
) -> None:
    """This function writes all important parameters into a .csv that were used
        during the method optimization step.

    Parameters:
    opt_result (pd.DataFrame): a dictionary containing all input and output
        information used and created by the function skopt.gp_minimize.
    optimizer_parameters (str): parameter, which was selected for optimization
        (e.g., all precursors, all doubly charged precursors).
    method_conf (dict): this dictionary contains all input parameters for all
        sub-functions.

    """
    df_tested_parameters = pd.DataFrame({
        "Tested_A1_A2_B1_B2": opt_result.x_iters,
        "result": opt_result.func_vals
    })
    df_tested_parameters[optimizer_parameters["evaluation_parameter"]] = df_tested_parameters["result"].apply(lambda x: 1.0 / x)
    df_tested_parameters.to_csv(
        method_conf["input"]["save_at"] + '/optimization_plots/tested_polygon_parameters.csv',
        index=False)
