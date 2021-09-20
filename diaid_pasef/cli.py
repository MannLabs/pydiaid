#!python


# external
import click

# local
import diaid_pasef

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


@run.command("create", help="Create window scheme.")
@click.option(
    "parameter_file_name",
    "-p",
    help=f"Parameter file (check out {DEFAULT_FILE} for an example)",
    default=DEFAULT_FILE,
    required=True,
)
def create_window_scheme(parameter_file_name):
    import diaid_pasef.main
    import json
    print(f"Using parameter file {parameter_file_name}")
    with open(parameter_file_name, "r") as infile:
        method_conf = json.load(infile)
    diaid_pasef.main.run_all(method_conf)
