#!python


# external
import click

# local
import pydiaid


@click.group(
    context_settings=dict(
        help_option_names=['-h', '--help'],
    ),
    invoke_without_command=True
)
@click.pass_context
@click.version_option(pydiaid.__version__, "-v", "--version")
def run(ctx, **kwargs):
    name = f"py_diAID {pydiaid.__version__}"
    click.echo("*" * (len(name) + 4))
    click.echo(f"* {name} *")
    click.echo("*" * (len(name) + 4))
    if ctx.invoked_subcommand is None:
        click.echo(run.get_help(ctx))


@run.command("gui", help="Start graphical user interface.")
def gui():
    import pydiaid.gui
    pydiaid.gui.run()
