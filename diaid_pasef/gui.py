#!python
import sys
import os
import time
import logging

# visualization
import panel as pn
import bokeh.server.views.ws

# modules
import diaid_pasef

BASE_PATH = os.path.dirname(__file__)
IMG_PATH = os.path.join(BASE_PATH, "img")
STYLE_PATH = os.path.join(BASE_PATH, "style")
DOCS_PATH = os.path.join(BASE_PATH, "docs")


class HeaderWidget(object):
    """This class creates a layout for the header of the dashboard with the name of the tool and all links to the MPI website, the MPI Mann Lab page and the GitHub repo.

    Parameters
    ----------
    title : str
        The name of the tool.
    img_folder_path : str
        The path to the folder with all necessary image files.
    github_url : str
        The path to the GitHub repository of the project.

    Attributes
    ----------
    header_title : pn.pane.Markdown
        A Panel Markdown pane that returns the title of the tool.
    mpi_biochem_logo : pn.pane.PNG
        A Panel PNG pane that embeds a png image file of the MPI Biochemisty logo and makes the image clickable with the link to the official website of the department.
    mpi_logo : pn.pane.JPG
        A Panel JPG pane that embeds a jpg image file of the MPI Biochemisty logo and makes the image clickable with the link to the official website.
    github_logo : pn.pane.PNG
        A Panel PNG pane that embeds a png image file of the GitHub logo and makes the image clickable with the link to the GitHub repository of the project.
    """
    def __init__(
        self,
        title,
        img_folder_path,
        github_url
    ):
        self.header_title = pn.pane.Markdown(
            f'# {title}',
            sizing_mode='stretch_width',
        )
        self.biochem_logo_path = os.path.join(
            img_folder_path,
            "mpi_logo.png"
        )
        self.mpi_logo_path = os.path.join(
            img_folder_path,
            "max-planck-gesellschaft.jpg"
        )
        self.github_logo_path = os.path.join(
            img_folder_path,
            "github.png"
        )
        self.mpi_biochem_logo = pn.pane.PNG(
            self.biochem_logo_path,
            link_url='https://www.biochem.mpg.de/mann',
            width=60,
            height=60,
            align='start'
        )
        self.mpi_logo = pn.pane.JPG(
            self.mpi_logo_path,
            link_url='https://www.biochem.mpg.de/en',
            height=62,
            embed=True,
            width=62,
            margin=(5, 0, 0, 5),
            css_classes=['opt']
        )
        self.github_logo = pn.pane.PNG(
            self.github_logo_path,
            link_url=github_url,
            height=70,
            align='end'
        )

    def create(self):
        layout = pn.Row(
            self.mpi_biochem_logo,
            self.mpi_logo,
            self.header_title,
            self.github_logo,
            height=73,
            sizing_mode='stretch_width'
        )
        return layout


class MainWidget(object):
    """This class creates a layout for the main part of the dashboard with the description of the tool and a button to download the manual for the project's GUI.

    Parameters
    ----------
    description : str
        The short description of the tool.
    manual_path : str
        The path to the GUI manual.

    Attributes
    ----------
    project_description : pn.pane.Markdown
        A Panel Markdown pane that shows the description of the project.
    manual : pn.widgets.FileDownload
        A Panel FileDownload widget that allows to download the GUI manual of the tool.
    """
    def __init__(
        self,
        description,
        manual_path
    ):
        self.project_description = pn.pane.Markdown(
            description,
            margin=(10, 0, 10, 0),
            css_classes=['main-part'],
            align='start',
            width=460
        )
        self.manual = pn.widgets.FileDownload(
            file=manual_path,
            label='Download Manual',
            button_type='default',
            align='center',
            auto=True,
            height=31,
            width=200,
            margin=(0, 20, 0, 0)
        )

    def create(self):
        layout = pn.Row(
            self.project_description,
            pn.layout.HSpacer(width=500),
            self.manual,
            background='#eaeaea',
            align='center',
            sizing_mode='stretch_width',
            height=190,
            margin=(10, 8, 10, 8),
            css_classes=['background']
        )
        return layout


class LoadLibraryCard(object):
    pass


class CalculateMzIMCard(object):
    pass


class OptimizationCard(object):
    pass


class FinalMethodCard(object):
    pass


def get_css_style(
    file_name="dashboard_style.css",
    directory=STYLE_PATH
):
    file = os.path.join(
        directory,
        file_name
    )
    with open(file) as f:
        return f.read()

def init_panel():
    pn.extension(raw_css=[get_css_style()])
    # pn.extension('plotly')


class GUI(object):
    # TODO: move to alphabase and docstring

    def __init__(
        self,
        name,
        img_folder_path,
        github_url,
        run_in_background=False,
        automatic_close=True,
    ):
        self.name = name
        self.tab_counter = 0
        self.header = HeaderWidget(
            name,
            img_folder_path,
            github_url
        )
        self.layout = pn.Column(
            self.header.create(),
            sizing_mode='stretch_width',
            min_width=1270
        )
        self.run_in_background = run_in_background
        self.automatic_close = automatic_close

    def start_server(self, run_in_background=False):
        if self.automatic_close:
            bokeh_ws_handler = bokeh.server.views.ws.WSHandler
            self.bokeh_server_open = bokeh_ws_handler.open
            bokeh_ws_handler.open = self.__open_browser_tab(
                self.bokeh_server_open
            )
            self.bokeh_server_on_close = bokeh_ws_handler.on_close
            bokeh_ws_handler.on_close = self.__close_browser_tab(
                self.bokeh_server_on_close
            )
        self.server = self.layout.show(threaded=True, title=self.name)
        if not run_in_background:
            self.server.join()
        elif not self.run_in_background:
            self.server.join()

    def __open_browser_tab(self, func):
        def wrapper(*args, **kwargs):
            self.tab_counter += 1
            return func(*args, **kwargs)
        return wrapper

    def __close_browser_tab(self, func):
        def wrapper(*args, **kwargs):
            self.tab_counter -= 1
            return_value = func(*args, **kwargs)
            if self.tab_counter == 0:
                self.stop_server()
            return return_value
        return wrapper

    def stop_server(self):
        logging.info("Stopping server...")
        self.server.stop()
        if self.automatic_close:
            bokeh_ws_handler = bokeh.server.views.ws.WSHandler
            bokeh_ws_handler.open = self.bokeh_server_open
            bokeh_ws_handler.on_close = self.bokeh_server_on_close


class DiAIDPasefGUI(GUI):
    # TODO: docstring

    def __init__(self, start_server=False):
        super().__init__(
            name=f"diAID-PASEF {diaid_pasef.__version__}",
            img_folder_path=IMG_PATH,
            github_url='https://github.com/MannLabs/diaid_pasef',
        )
        self.project_description = """### Include the description."""
        self.manual_path = os.path.join(
            DOCS_PATH,
            "manual.pdf"
        )
        self.main_widget = MainWidget(
            self.project_description,
            self.manual_path
        )

        # ERROR/WARNING MESSAGES
        self.error_message_upload = "The selected file can't be uploaded. Please check the instructions for data uploading."

        # self.library = LoadLibraryCard()
        self.layout += [
            self.main_widget.create(),
            # self.LoadLibraryCard.create(),
        ]
        if start_server:
            self.start_server()


def run():
    # set style
    init_panel()
    DiAIDPasefGUI(start_server=True)


if __name__ == '__main__':
    run()
