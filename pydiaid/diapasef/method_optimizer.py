from pydiaid.diapasef  import method_creator
from pydiaid.diapasef  import method_evaluation
from pydiaid.diapasef  import graphs

# for data manipulation
import pandas as pd

# importing for scientific and numeric manipulations
import numpy as np
import skopt


def optimization(
    library: pd.DataFrame,
    method_parameters: dict,
    xi: np.ndarray,
    yi: np.ndarray,
    zi: np.ndarray,
    method_conf: dict,
    optimizer_parameters: dict,
) -> 'scipy.optimize.optimize.OptimizeResult':
    """ Optimize dia-PASEF method parameters: Bayesian optimization using a
    Gaussian process. The output values should follow a multivariate gaussian
    curve.

    Parameters:
    library (pd.DataFrame): a pre-filtered data frame with unified column names
        containing all required precursor information.
    method_parameters (dict): dictionary, which includes all input parameters
        for creating a dia-PASEF method.
    xi, yi, zi (numpy.ndarray): coordinates of the kernel density estimation
        where zi indicates the density
    method_conf (dict): this dictionary contains all input parameters for all
        sub-functions.
    optimizer_parameters (dict): dictonary, which includes all input parameters
        for optimizing a dia-PASEF method in dependence of a proteomics library.

    Returns:
    scipy.optimize.optimize.OptimizeResult: this parameter contains important
    information regarding the optimization process (e.g., start settings,
    tested parameters and test results)
    """

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
            (
            optimizer_parameters["YA2"][0]-optimizer_parameters["YA1"][0],
            optimizer_parameters["YA2"][1]-optimizer_parameters["YA1"][0]
            ),
            (
            optimizer_parameters["YB1"][0]-optimizer_parameters["YA1"][0],
            optimizer_parameters["YB1"][1]-optimizer_parameters["YA1"][0]
            ),
            (
            optimizer_parameters["YB2"][0]-optimizer_parameters["YA1"][0],
            optimizer_parameters["YB2"][1]-optimizer_parameters["YA1"][0]
            )
        ],
        #n_random_starts = optimizer_parameters["n_random_starts"],
        n_calls = optimizer_parameters["n_calls"],  # n to test
        #optimizer_parameters["n_start"],  # random starts
        n_initial_points = optimizer_parameters["initial_points"],  # initial points
        )

    optimized_results = [opt_result.x[0],
        opt_result.x[0]+opt_result.x[1],
        opt_result.x[0]+opt_result.x[2],
        opt_result.x[0]+opt_result.x[3]]

    print("########")
    print("BEST RESULT")
    print("INPUT: " + str(
        [opt_result.x[0],
        opt_result.x[0]+opt_result.x[1],
        opt_result.x[0]+opt_result.x[2],
        opt_result.x[0]+opt_result.x[3]]
        ))
    print("OUTPUT: " + str(1.0 / opt_result.fun))
    print("########")

    write_tested_parameters_into_csv(
        opt_result,
        optimizer_parameters,
        method_conf
        )

    # return opt_result
    return optimized_results


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
    """This function contains one optimization step. First a dia-PASEF scheme
        is calculated based on the scan area coordinates, second the evaluation
        parameter is calculated to estimate the optimization potential, last
        the dia-PASEF method scheme is printed on top of the kernel density
        estimation of the used proteomics library to follow the optimization
        process visually.

    Parameters:
    library (pd.DataFrame): a pre-filtered data frame with unified column names
        containing all required precursor information.
    method_parameters (dict): dictionary, which includes all input parameters
        for creating a dia-PASEF method.
    evaluation_parameter (str): parameter, which was selected for optimization
        (e.g., all precursors, all doubly charged precursors).
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
    #library_reduced = library[library['Charge'] != 1]

    df_parameters_final = method_creator.method_creation(
        #library_reduced["mz"],
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
    print("RUN WITH: " + str(
    [dim[0],
    dim[0]+dim[1],
    dim[0]+dim[2],
    dim[0]+dim[3]]
    ) + " | RESULT: " + str(1.0 / result))

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


def write_tested_parameters_into_csv(
    opt_result: pd.DataFrame,  # scipy.optimize.optimize.OptimizeResult
    optimizer_parameters: str,
    method_conf: dict,
) -> None:
    """This function writes all important parameters that were used during the
    method optimization step into a .csv file.

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


def find_matching_filename_index(filenames, parameters):
    """
    Find the index of a filename that matches the given parameters.
    
    Args:
        filenames (list): List of filenames to search through
        parameters (list): List of parameter values [A1, A2, B1, B2]
    
    Returns:
        int: Index of matching filename or -1 if no match found
    """
    # Create the search pattern
    param_names = ['A1', 'A2', 'B1', 'B2']
    search_pattern = '_'.join(f"{name}_{value}" for name, value in zip(param_names, parameters))
    
    # Search through filenames
    for index, filename in enumerate(filenames):
        # Remove potential spaces in filename for consistent comparison
        cleaned_filename = filename.replace(' ', '')
        if search_pattern in cleaned_filename:
            return index
            
    return -1
