#!python


# external
import click

# local
import diaid_pasef
from diaid_pasef import loader
from diaid_pasef import main
from diaid_pasef import method_optimizer

# for data manipulation:
import pandas as pd

# for saving input data as json:
import json

# for creating new folders:
import os

BASE_PATH = os.path.dirname(__file__)
LIB_PATH = os.path.join(BASE_PATH, "lib")
DEFAULT_FILE = os.path.join(
    LIB_PATH,
    "default_parameters.json"
)


@click.group(
    context_settings=dict(
        help_option_names=['-h', '--help'],
    ),
    invoke_without_command=True
)
@click.pass_context
@click.version_option(diaid_pasef.__version__, "-v", "--version")
def run(ctx, **kwargs):
    name = f"diAID-PASEF {diaid_pasef.__version__}"
    click.echo("*" * (len(name) + 4))
    click.echo(f"* {name} *")
    click.echo("*" * (len(name) + 4))
    if ctx.invoked_subcommand is None:
        click.echo(run.get_help(ctx))


@run.command("gui", help="Start graphical user interface.")
def gui():
    import diaid_pasef.gui
    diaid_pasef.gui.run()


@run.command("optimize", help="Optimize a dia-PASEF method.")
@click.option(
    "parameter_file_name",
    "-p",
    help=f"Parameter file (check out {DEFAULT_FILE} for an example)",
    default=DEFAULT_FILE,
    required=True,
)
def optimize_window_scheme(parameter_file_name):
    import diaid_pasef.main
    import json
    print(f"Using parameter file {parameter_file_name}")
    with open(parameter_file_name, "r") as infile:
        method_conf = json.load(infile)
    run_all(method_conf)


@run.command("create", help="Create a specific dia-PASEF method.")
@click.option(
    "parameter_file_name",
    "-p",
    help=f"Parameter file (check out {DEFAULT_FILE} for an example)",
    default=DEFAULT_FILE,
    required=True,
)
def create_diaPASEF_method(parameter_file_name,):
    import diaid_pasef.main
    import json
    print(f"Using parameter file {parameter_file_name}")
    with open(parameter_file_name, "r") as infile:
        method_conf = json.load(infile)
    library_plus_create_methods(method_conf)


@run.command("evaluate", help="Evaluate a dia-PASEF method.")
@click.option(
    "parameter_file_name",
    "-p",
    help=f"Parameter file (check out {DEFAULT_FILE} for an example)",
    default=DEFAULT_FILE,
    required=True,
)
def evaluate_diaPASEF_method(parameter_file_name,):
    import diaid_pasef.main
    import json
    print(f"Using parameter file {parameter_file_name}")
    with open(parameter_file_name, "r") as infile:
        method_conf = json.load(infile)
    library_plus_evaluate_method(method_conf)


@run.command("charge", help="Evaluate a dia-PASEF method for multiple charged precursors.")
@click.option(
    "parameter_file_name",
    "-p",
    help=f"Parameter file (check out {DEFAULT_FILE} for an example)",
    default=DEFAULT_FILE,
    required=True,
)
def evaluate_multiple_charged_prec(parameter_file_name,):
    import diaid_pasef.main
    import json
    print(f"Using parameter file {parameter_file_name}")
    with open(parameter_file_name, "r") as infile:
        method_conf = json.load(infile)
    multiple_charged_prec(method_conf)


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
        save_at,
        save_at+'/optimization_plots',
        save_at+'/input_library',
        save_at+'/final_method'
        ]
    main.create_folder(folder_paths)

    # save input parameters as json file for reusing
    out_file = open(save_at+"/input_parameters.json", "w")
    json.dump(method_conf, out_file, indent = 4)
    out_file.close()

    xi, yi, zi, library = main.library_information(
        method_conf,
        save_at
        )

    main.precursor_within_scan_area(
        library,
        method_parameters,
        save_at
        )

    opt_result = method_optimizer.optimization(
        library,
        method_parameters,
        xi,
        yi,
        zi,
        method_conf,
        optimizer_parameters
        )

    dim = [
        opt_result.x[0],
        opt_result.x[0] + opt_result.x[1],
        opt_result.x[0] + opt_result.x[2],
        opt_result.x[0] + opt_result.x[3]
        ]

    df_parameters_final = main.create_final_method(
        library,
        method_parameters,
        dim,
        method_conf
        )

    main.final_method_information(
        df_parameters_final,
        xi,
        yi,
        zi,
        method_conf,
        save_at, library,
        method_parameters,
        dim
        )


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
    main.create_folder(folder_paths)

    # save input parameters as json file for reusing
    out_file = open(save_at+"/input_parameters.json", "w")
    json.dump(method_conf, out_file, indent = 4)
    out_file.close()

    df_parameters_final = pd.read_csv(
        method_conf["input"]["diaPASEF_method_only_used_for_evaluate"],
        skiprows=4,
        names=["MS Type", "Cycle Id", "Start IM", "End IM", "Start Mass", "End Mass", "CE"]
        )

    xi, yi, zi, library = main.library_information(
        method_conf,
        save_at
        )

    main.final_method_information(
        df_parameters_final,
        xi,
        yi,
        zi,
        method_conf,
        save_at,
        library,
        method_parameters,
        None
        )

    print("Method evaluated")


def library_plus_create_methods(
    method_conf: dict
) -> None:
    method_parameters = method_conf["method_parameters"]
    save_at = method_conf["input"]["save_at"]
    dim = method_parameters["scan_area_A1_A2_B1_B2_only_used_for_create"]

    folder_paths = [
        save_at,
        save_at+'/final_method'
        ]
    main.create_folder(folder_paths)

    # save input parameters as json file for reusing
    out_file = open(save_at+"/input_parameters.json", "w")
    json.dump(method_conf, out_file, indent = 4)
    out_file.close()

    library = loader.load_library(
        method_conf["input"]["library_name"],
        method_conf["input"]["analysis_software"],
        method_conf["input"]["PTM"]
        )

    main.precursor_within_scan_area(
        library,
        method_parameters,
        save_at
        )

    main.create_final_method(
        library,
        method_parameters,
        dim,
        method_conf
        )

    print("Method created")


def multiple_charged_prec(
    method_conf: dict
) -> None:

    save_at = method_conf["input"]["save_at"]

    folder_paths = [
        save_at,
        save_at+'/evaluation_multiple_charged_prec'
        ]
    main.create_folder(folder_paths)

    # save input parameters as json file for reusing
    out_file = open(save_at+"/input_parameters.json", "w")
    json.dump(method_conf, out_file, indent = 4)
    out_file.close()

    df_parameters_final = pd.read_csv(
        method_conf["input"]["diaPASEF_method_only_used_for_evaluate"],
        skiprows=4,
        names=["MS Type", "Cycle Id", "Start IM", "End IM", "Start Mass", "End Mass", "CE"]
        )

    library = loader.load_library(
        method_conf["input"]["library_name"],
        method_conf["input"]["analysis_software"],
        method_conf["input"]["PTM"]
        )

    main.evaluate_for_multiple_charged_prec(
        method_conf,
        save_at,
        library,
        df_parameters_final
    )
