#!python

import os
import logging

# visualization
import panel as pn
import bokeh.server.views.ws

import holoviews as hv


#local
# import pydiaid.synchropasef as pydiaid
from pydiaid.synchropasef.gui_synchropasef import SynchroCardWidget
from pydiaid.diapasef.gui_diapasef import DiaCardWidget
from pydiaid.oadia.gui_orbitrapastraldia import OrbitrapAstralCardWidget
# from pydiaid.quadrupolecalibration.gui_quad import QuadCardWidget

import warnings
warnings.filterwarnings("ignore")

pn.extension(sizing_mode="stretch_width")


# paths
BASE_PATH = os.path.dirname(__file__)
STYLE_PATH = os.path.join(BASE_PATH, "style")
DOCS_PATH = os.path.join(BASE_PATH, "docs")
IMG_PATH = os.path.join(BASE_PATH, "img")
LATEST_GITHUB_INIT_FILE = "https://github.com/MannLabs/pydiaid/tree/main/pydiaid/__init__.py"


def get_css_style(
    file_name="dashboard_style.css",
    directory=STYLE_PATH
):
    """
    The get_css_style function is used to read and return the content of a specified CSS style file.

    Parameters:
    file_name (str): The name of the CSS file to be read. By default, it is set to “dashboard_style.css”.
    directory (str): The directory where the CSS file is located. By default, it is set to the constant 
        STYLE_PATH.

    Returns:
    A string containing the content of the specified CSS file.
    """
    file = os.path.join(
        directory,
        file_name
    )
    with open(file) as f:
        return f.read()


def init_panel():
    """
    The init_panel function is used to initialize the Panel and other extensions.

    This function initializes the Panel with custom CSS style from 'get_css_style' function and also 
    extends Panel to have 'bokeh', 'plotly', and 'tabulator' extensions.

    Returns:
    None. This function doesn't return anything. It is an initial setup process and all changes are 
    performed in the global context.
    """
    pn.extension(raw_css=[get_css_style()])
    hv.extension('bokeh')
    pn.extension('plotly')
    pn.extension('tabulator')


class HeaderWidget(object):
    """This class creates a layout for the header of the dashboard
    with the name of the tool and all links to the MPI website,
    the MPI Mann Lab page and the GitHub repo.
    Parameters
    ----------
    title : str
        The name of the tool.
    Attributes
    ----------
    header_title : pn.pane.Markdown
        A Panel Markdown pane that returns the title of the tool.
    mpi_biochem_logo : pn.pane.PNG
        A Panel PNG pane that embeds a png image file of the MPI
        Biochemisty logo and makes the image clickable with the
        link to the official website.
    mpi_logo : pn.pane.JPG
        A Panel JPG pane that embeds a jpg image file of the MPI
        Biochemisty logo and makes the image clickable with the
        link to the official website.
    github_logo : pn.pane.PNG
        A Panel PNG pane that embeds a png image file of the
        GitHub logo and makes the image clickable with the
        link to the GitHub repository of the project.
    """

    def __init__(
        self,
        title,
        img_folder_path,
        github_url
    ):
        self.layout = None
        self.header_title = pn.pane.Markdown(
            f'# {title}',
            sizing_mode='stretch_width',
            align='center'
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

    def create_layout(self):
        self.layout = pn.Row(
            self.mpi_biochem_logo,
            self.mpi_logo,
            pn.Spacer(sizing_mode='stretch_width'),
            self.header_title,
            pn.Spacer(sizing_mode='stretch_width'),
            self.github_logo,
            height=100,
            sizing_mode='stretch_width',
        )
        return self.layout


class GUI(object):
    """
    The GUI class creates a web interface for visualization with start/stop server functionality and tab 
    management.

    Instances of GUI are initialized with a name, GitHub URL, and flags to specify if the server should 
    run in background or close automatically.

    Attributes:
    name (str): Name of the GUI.
    tab_counter (int): Counter to track number of open tabs.
    header: Instance of HeaderWidget.
    layout: Structure of the GUI, comprises the header at initialization.
    run_in_background (bool): Specifies if the server should run in the background. Default is False.
    automatic_close (bool): Specifies if the server should close automatically. Default is True.
    server: Initialized later in the start_server function.

    Methods:
    init: Constructor for the LoadLibraryCard class.
    start_server: Starts the server.
    __open_browser_tab: Intercepts opening of new tabs, incrementing tab_counter.
    __close_browser_tab: Intercepts closing of tabs, decrementing tab_counter and stopping server if 
        tab_counter reaches 0.
    stop_server: Stops the server and resets tab management if automatic_close is True.
    """
    def __init__(
        self,
        name,
        github_url,
        run_in_background=False,
        automatic_close=True,
    ):
        self.name = name
        self.tab_counter = 0
        self.header = HeaderWidget(
            name,
            IMG_PATH,
            github_url
        )
        self.layout = pn.Column(
            self.header.create_layout(),
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


# this code was taken from the AlphaTims Python package (https://github.com/MannLabs/alphatims/blob/master/alphatims/utils.py) and modified
def check_github_version(silent=False) -> str:
    """Checks and logs the current version of py_diAID.
    Check if the local version equals the py_diAID GitHub main branch.
    This is only possible with an active internet connection and
    if no credentials are required for GitHub.
    Parameters
    ----------
    silent : str
        Use the logger to display the obtained conclusion.
        Default is False.
    Returns
    -------
    : str
        The version on the py_diAID GitHub main branch.
        "" if no version can be found on GitHub
    """
    import urllib.request
    import urllib.error
    try:
        with urllib.request.urlopen(LATEST_GITHUB_INIT_FILE) as version_file:
            for line in version_file.read().decode('utf-8').split("\n"):
                if line.startswith("__version__"):
                    github_version = line.split()[2][1:-1]
                    if not silent:
                        if github_version != pydiaid.__version__:
                            logging.info(
                                f"You are currently using py_diAID version "
                                f"{pydiaid.__version__}. "
                                f"However, the latest version of py_diAID on "
                                f"GitHub is {github_version}. Checkout "
                                "https://github.com/MannLabs/pydiaid.git "
                                "for instructions on how to update py_diAID"
                                "..."
                            )
                            logging.info("")
                        else:
                            logging.info(
                                "Current py_diAID version is up-to-date "
                                "with GitHub."
                            )
                            logging.info("")
                    return github_version
    except IndexError:
        logging.info(
            "Could not check GitHub for the latest py_diAID release."
        )
        logging.info("")
        return ""
    except urllib.error.URLError:
        logging.info(
            "Could not check GitHub for the latest py_diAID release."
        )
        logging.info("")
        return ""


class MainWidget(object):
    """This class create a layout for the main part of the
    dashboard with the description of the tool and a button
    to download the manual for the project's GUI.
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
        A Panel FileDownload widget that allows to download the GUI
        manual of the tool.
    """
    def __init__(
        self,
        description,
        manual_path
    ):
        self.layout = None
        self.project_description = pn.pane.Markdown(
            description,
            margin=(20, 0, 10, 0),
            css_classes=['main-part'],
            align='start',
            width=620
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
        self.download_new_version_button = pn.widgets.Button(
            button_type="danger",
            align='center',
            height=31,
            width=200,
            margin=(20, 20, 0, 0)
        )

    def create_layout(self):
        latest_github_version = check_github_version(
           silent=False
        )
        if latest_github_version and \
                latest_github_version != pydiaid.__version__:
            self.download_new_version_button.name = f"Download version {latest_github_version}"
            download_new_version_button = self.download_new_version_button
            download_new_version_button.js_on_click(
                code="""window.open("https://github.com/MannLabs/pydiaid/releases/latest")"""
            )
        else:
            download_new_version_button = None

        self.layout = pn.Row(
            pn.Spacer(sizing_mode='stretch_width'),
            self.project_description,
            pn.Spacer(sizing_mode='stretch_width'),
            pn.Column(
                self.manual,
                download_new_version_button,
                align='center',
            ),
            background='#eaeaea',
            align='center',
            sizing_mode='stretch_width',
            height=265,
            margin=(10, 8, 10, 8),
            css_classes=['background']
        )
        return self.layout


class BaseWidget(object):
    """
    The BaseWidget class serves as the base for all widget classes, providing basic structure and 
    functionality.

    Each instance of BaseWidget is initialized with a name and contains an update event for managing 
    dependant updates.

    Attributes:
    name (str): Name of the widget.
    __update_event: An integer input widget to keep track of update events.
    depends: A Panel Depends object that wraps the value of the __update_event.
    active_depends: A Panel Depends object that wraps the value of the __update_event and watches for 
        changes.

    Methods:
    init: Constructor for the LoadLibraryCard class.
    trigger_dependancy: Increments the value of __update_event by 1, effectively triggering a change 
        event for all dependants.
    """
    def __init__(self, name):
        self.name = name
        self.__update_event = pn.widgets.IntInput(value=0)
        self.depends = pn.depends(self.__update_event.param.value)
        self.active_depends = pn.depends(
            self.__update_event.param.value,
            watch=True
        )

    def trigger_dependancy(self):
        self.__update_event.value += 1


class TabsWidget(object):
    """
    The TabsWidget class used to create a Tabs layout with different widgets.
    This class creates a Tabs layout with different widget layouts. These layouts are specified in the 
    create_layout method.

    Attributes:
    layout: Updated in the create_layout method with the Panel with Tabs layout.

    Methods:
    init: Constructor for the LoadLibraryCard class.
    create_layout: Creates and returns the panel layout for the class as tabs.

    Note:
    The create_layout method takes an optional parameter 'tab_list'. 
    """
    def __init__(self):
        self.layout = None

    def create_layout(
        self,
        tab_list=None
    ):
        self.layout = pn.Tabs(
            background='#eaeaea',
            tabs_location='above',
            margin=(10, 8, 10, 8),
            sizing_mode='stretch_width',
            css_classes=['background'],
            )
        self.layout.extend(
        [
            ('Optimal dia-PASEF', DiaCardWidget().create_layout()),
            ('Synchro-PASEF', SynchroCardWidget().create_layout()),
            ('Orbitrap Astral DIA', OrbitrapAstralCardWidget().create_layout())
        ]
        )
        return self.layout


class pydiaidGui(GUI):
    """
    PydiaidGui class extends GUI to create the user interface for the py_diAID tool.
    This class is responsible for creating and rendering the user interface for the py_diAID tool. It 
    includes a description, a manual, and tabs for different functionalities.

    Attributes:
    project_description (str): Description of the py_diAID tool.
    manual_path (str): Path to the manual.
    main_widget: Instance of MainWidget.
    tabs: Instance of TabsWidget.

    Inherited attributes:
    name, layout, automatic_close, server from GUI class.

    Methods:
    Inherited methods start_server and stop_server from GUI class.
    """
    def __init__(self, start_server=False):
        super().__init__(
            name="py_diAID",
            github_url='https://github.com/MannLabs/pydiaid',
        )
        self.project_description = """#### py_diAID is a Python tool that automatically and optimally places DIA (Data-Independent Acquisition) window schemes for efficient precursor coverage. Using pre-acquired precursor information, it generates dia-PASEF, synchro-PASEF, and Orbitrap Astral DIA methods. The name diAID stands for Automated Isolation Design for DIA.\n <i style="font-size: 0.8em; display: block">  Please cite: Skowronek et al., Mann, MCP, 2022 for optimal dia-PASEF methods, \n Skowronek et al., Willems, Raether, Mann, MCP, 2023 for synchro-PASEF, \n Skowronek et al., Mann, Nat Protoc, 2025 for PASEF workflows and py_diAID.</i>"""

        self.manual_path = os.path.join(
            DOCS_PATH,
            "manual.pdf"
        )
        self.main_widget = MainWidget(
            self.project_description,
            self.manual_path
        )
        self.tabs = TabsWidget()

        self.header = HeaderWidget(
            "py_diAID",
            IMG_PATH,
            'https://github.com/MannLabs/pydiaid'
        )

        self.layout = pn.Column(
            self.header.create_layout(),
            self.main_widget.create_layout(),
            self.tabs.create_layout([
                ('Optimal dia-PASEF', DiaCardWidget().create_layout()),
                ('Synchro-PASEF', SynchroCardWidget().create_layout()),
                ('Orbitrap Astral DIA', OrbitrapAstralCardWidget().create_layout())
            ])
        )
        if start_server:
            self.start_server()


def run():
    """
    This function initializes the Panel using the init_panel method and starts the pydiAID GUI.
    Returns: None. 
    """
    init_panel()
    gui = pydiaidGui(start_server=False)
    
    # Run the GUI on the main thread
    pn.serve(gui.layout, show=True, title="py_diAID")

if __name__ == "__main__":
    run()