#!python
import sys
import os
import re
import time
import json
import logging
import platform
import pandas as pd
import numpy as np

# visualization
import panel as pn
import bokeh.server.views.ws

# modules
import diaid_pasef
import diaid_pasef.loader
import diaid_pasef.main
import diaid_pasef.graphs

import warnings
warnings.filterwarnings("ignore")

# paths
BASE_PATH = os.path.dirname(__file__)
IMG_PATH = os.path.join(BASE_PATH, "img")
STYLE_PATH = os.path.join(BASE_PATH, "style")
DOCS_PATH = os.path.join(BASE_PATH, "docs")
DEFAULT_FILE = os.path.join(
    BASE_PATH,
    'lib',
    "default_parameters.json"
)

with open(DEFAULT_FILE, "r") as infile:
    method_conf = json.load(infile)

if platform.system() == 'Windows':
    method_path_placeholder = 'D:\diaid_pasef\static\DIAParameterspy3TC.txt'
    library_path_placeholder = 'D:\diaid_pasef\static\AlphaPept_results.csv'
    save_folder_placeholder = 'D:\diaid_pasef\static'
else:
    method_path_placeholder = '/Users/diaid_pasef/static/DIAParameterspy3TC.txt'
    library_path_placeholder = '/Users/diaid_pasef/static/AlphaPept_results.csv'
    save_folder_placeholder = '/Users/diaid_pasef/static'


class BaseWidget(object):
    # TODO: docstring
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


class LoadLibraryCard(BaseWidget):
    # TODO: docstring
    def __init__(self):
        super().__init__(name="Data")
        self.library = None
        self.layout = None
        self.path_library = pn.widgets.TextInput(
            name='Specify the path to the library:',
            placeholder=library_path_placeholder,
            value='/Users/eugeniavoytik/Projects/diaid_pasef/diaid_pasef/static/AlphaPept_results.csv/',
            width=900,
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )
        self.path_save_folder = pn.widgets.TextInput(
            name='Save the output to the following folder:',
            value='/Users/eugeniavoytik/Projects/diaid_pasef/diaid_pasef/static/',
            placeholder=save_folder_placeholder,
            width=900,
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )
        self.ptm = pn.widgets.TextInput(
            name='Specify the PTM:',
            value='None',
            placeholder='Phospho',
            width=680,
            sizing_mode='stretch_width',
            margin=(15, 15, 15, 15)
        )
        self.analysis_software = pn.widgets.Select(
            name='Analysis software',
            value='MaxQuant',
            options=['AlphaPept', 'MaxQuant', 'MS_Fragger', 'Spectronaut_single-shot', 'Spectronaut_library'],
            width=200,
            margin=(15, 15, 15, 15)
        )
        # UPLOAD DATA
        self.upload_button = pn.widgets.Button(
            name='Upload library',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )
        self.upload_progress = pn.indicators.Progress(
            # max=1,
            # value=1,
            active=False,
            bar_color='light',
            width=250,
            align='center',
            margin=(-10, 0, 30, 0)
        )
        self.import_error = pn.pane.Alert(
            alert_type="danger",
            sizing_mode='stretch_width',
            # object='test warning message',
            margin=(30, 15, 5, 15),
        )

    def create(self):
        self.layout = pn.Card(
            pn.Row(
                pn.Column(
                    self.path_library,
                    self.path_save_folder,
                    pn.Row(
                        self.ptm,
                        self.analysis_software,
                    ),
                    margin=(10, 30, 10, 10),
                ),
                pn.Spacer(sizing_mode='stretch_width'),
                pn.Column(
                    self.upload_button,
                    self.upload_progress,
                    # self.import_error,
                    align='center',
                    margin=(100, 40, 0, 0),
                )
            ),
            pn.layout.Divider(
                sizing_mode='stretch_width',
                margin=(0, 10, -20, 10),
            ),
            pn.Row(
                None, None, None
            ),
            title='Load Library',
            collapsed=True,
            collapsible=True,
            header_background='#eaeaea',
            background ='white',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            # height=470,
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

        dependances = {
            self.path_library: [self.update_parameters, 'value'],
            self.path_save_folder: [self.update_parameters, 'value'],
            self.ptm: [self.update_parameters, 'value'],
            self.analysis_software: [self.update_parameters, 'value'],
            # self.upload_button: [self.upload_data, 'clicks'],

        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1],
                # onlychanged=False
            )
        self.upload_button.param.watch(
            self.upload_data,
            'clicks',
            onlychanged=True
        )
        return self.layout


    def update_parameters(self, event):
        global method_conf
        convertion_dict = {
            self.path_library.name: "library_name",
            self.path_save_folder.name: "save_at",
            self.ptm.name: "PTM",
            self.analysis_software.name: "analysis_software",
        }
        method_conf['input'][convertion_dict[event.obj.name]] = event.new


    def upload_data(self, *args):
        self.upload_progress.active = True
        self.library = diaid_pasef.loader.load_library(
            self.path_library.value,
            "AlphaPept",
            method_conf["input"]["PTM"]
        )
        self.upload_progress.active = False
        folder_paths = [
            self.path_save_folder.value,
            os.path.join(
                self.path_save_folder.value,
                'input_library'
            ),
            os.path.join(
                self.path_save_folder.value,
                'final_method'
            ),
        ]
        diaid_pasef.main.create_folder(folder_paths)
        self.xi, self.yi, self.zi = diaid_pasef.graphs.kernel_density_calculation(
            self.library,
            method_conf["graphs"]["numbins"]
        )
        self.layout[2][0] = pn.pane.Matplotlib(
            diaid_pasef.graphs.plot_precursor_distribution_as_histogram(
                self.library,
                method_conf["graphs"],
                os.path.join(
                    self.path_save_folder.value,
                    'input_library',
                    'Histogram_precursor_distribution_in_library.png'
                ),
                gui=True
            ),
            tight=True
        )
        self.layout[2][1] = pn.pane.Matplotlib(
            diaid_pasef.graphs.plot_density(
                self.xi,
                self.yi,
                self.zi,
                method_conf["graphs"],
                os.path.join(
                    self.path_save_folder.value,
                    'input_library',
                    'Kernel_density_distribution_library.png'
                ),
                gui=True
            ),
            tight=True
        )
        dict_charge_of_precursor = diaid_pasef.method_evaluation.calculate_percentage_multiple_charged_precursors(self.library)
        mult_charged_precursor_info = pd.DataFrame(
            {
                "column name": list(dict_charge_of_precursor.keys()),
                "column value": list(dict_charge_of_precursor.values())
            }
        )
        mult_charged_precursor_info.to_csv(
            os.path.join(
                self.path_save_folder.value,
                'input_library',
                'percentage_of_multiple_charged_precursors.csv'
            ),
            index=False
        )
        self.layout[2][2] = pn.Column(
            pn.pane.Markdown(
                '### Percentage of multiple charged precursors',
                 align='center'
            ),
            pn.widgets.DataFrame(
                mult_charged_precursor_info,
                autosize_mode='fit_viewport',
                auto_edit=False
            ),
            margin=(0, 50),
            sizing_mode='stretch_width',
            align='center'
        )
        self.trigger_dependancy()
        self.upload_progress.active = False


class SpecifyParametersCard(object):
    # TODO: docstring
    def __init__(self):
        self.mz = pn.widgets.EditableRangeSlider(
            name='M/z range',
            start=100,
            end=2000,
            value=(300, 1250),
            step=50,
            margin=(15, 15, 0, 15),
            width=430,
        )
        self.ion_mobility = pn.widgets.EditableRangeSlider(
            name='IM range',
            start=0.5,
            end=2.0,
            value=(0.6, 1.6),
            step=0.1,
            margin=(15, 15, 0, 15),
            width=430,
        )
        self.num_dia_pasef_scans = pn.widgets.IntInput(
            name='Number of DIA PASEF scans',
            start=1,
            end=20,
            value=12,
            step=1,
            margin=(15, 15, 0, 15),
            width=430,
        )
        self.im_steps = pn.widgets.IntInput(
            name='Number of IM steps',
            start=1,
            end=10,
            value=2,
            step=1,
            margin=(15, 15, 0, 15),
            width=430,
        )
        self.overlap = pn.widgets.IntInput(
            name='Overlap',
            start=0,
            end=10,
            value=0,
            step=1,
            margin=(15, 15, 0, 15),
            width=430,
        )
        self.shift_of_final_method = pn.widgets.FloatInput(
            name='Shift of final method',
            start=0.0,
            end=1.0,
            value=0.22,
            step=0.01,
            margin=(15, 15, 0, 15),
            width=430,
        )
        self.calculate_button = pn.widgets.Button(
            name='Calculate',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )


    def create(self):
        self.layout = pn.Card(
            pn.Row(
                pn.Column(
                    pn.WidgetBox(
                        pn.Row(
                            self.mz,
                            self.ion_mobility,
                            sizing_mode='stretch_width',
                        ),
                        pn.Row(
                            self.num_dia_pasef_scans,
                            self.im_steps,
                            sizing_mode='stretch_width'
                        ),
                        pn.Row(
                            self.overlap,
                            self.shift_of_final_method,
                            sizing_mode='stretch_width'
                        ),
                        sizing_mode='stretch_width',
                        margin=(20, 10, 30, 10),
                        height=220
                    ),
                    margin=(10, 30, 10, 10),
                ),
                pn.Spacer(sizing_mode='stretch_width'),
                pn.Column(
                    self.calculate_button,
                    # self.upload_progress,
                    # self.import_error,
                    align='center',
                    margin=(100, 130, 0, 0),
                )
            ),
            pn.layout.Divider(
                sizing_mode='stretch_width',
                margin=(-20, 10, -20, 10),
            ),
            pn.Row(
                None
            ),
            title='Specify Method Parameters',
            collapsed=True,
            collapsible=True,
            header_background='#eaeaea',
            background ='white',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            # height=470,
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )
        dependances = {
            self.mz: [self.update_parameters, 'value'],
            self.ion_mobility: [self.update_parameters, 'value'],
            self.num_dia_pasef_scans: [self.update_parameters, 'value'],
            self.im_steps: [self.update_parameters, 'value'],
            self.overlap: [self.update_parameters, 'value'],
            self.shift_of_final_method: [self.update_parameters, 'value'],
            self.calculate_button: [self.run_calculation, 'clicks'],
        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1],
                onlychanged=True
            )
        return self.layout

    def update_parameters(self, event):
        global method_conf
        convertion_dict = {
            self.mz.name: "mz",
            self.ion_mobility.name: "ion_mobility",
            self.num_dia_pasef_scans.name: "num_dia_pasef_scans",
            self.im_steps.name: "im_steps",
            self.overlap.name: "overlap",
            self.shift_of_final_method.name: "shift_of_final_method",
        }
        if isinstance(event.new, tuple):
            method_conf['method_parameters'][convertion_dict[event.obj.name]] = list(event.new)
        else:
            method_conf['method_parameters'][convertion_dict[event.obj.name]] = event.new

    def run_calculation(self, *args):
        pass



class OptimizationCard(BaseWidget):
    # TODO: docstring
    def __init__(self, data):
        super().__init__(name="Optimization")
        self.data = data
        self.opt_result_x = [0, 0, 0, 0]
        self.n_calls = pn.widgets.IntInput(
            name='Number of calls',
            start=1,
            end=10,
            value=5,
            step=1,
            margin=(15, 15, 0, 15),
            width=430,
        )
        self.n_start = pn.widgets.IntInput(
            name='n_start',
            start=1,
            end=10,
            value=1,
            step=1,
            margin=(15, 15, 0, 15),
            width=430,
        )
        self.initial_points = pn.widgets.IntInput(
            name='initial_points',
            start=1,
            end=10,
            value=1,
            step=1,
            margin=(15, 15, 0, 15),
            width=430,
        )
        self.YA1 = pn.widgets.EditableRangeSlider(
            name='YA1 range',
            start=0.2,
            end=1.0,
            value=(0.4, 0.6),
            step=0.1,
            margin=(15, 15, 0, 15),
            width=430,
        )
        self.YA2 = pn.widgets.EditableRangeSlider(
            name='YA2 range',
            start=0.6,
            end=1.2,
            value=(0.8, 0.95),
            step=0.05,
            margin=(15, 15, 0, 15),
            width=430,
        )
        self.YB1 = pn.widgets.EditableRangeSlider(
            name='YB1 range',
            start=0.6,
            end=1.5,
            value=(0.96, 1.2),
            step=0.1,
            margin=(15, 15, 0, 15),
            width=430,
        )
        self.YB2 = pn.widgets.EditableRangeSlider(
            name='YB2 range',
            start=1.0,
            end=1.8,
            value=(1.25, 1.5),
            step=0.05,
            margin=(15, 15, 0, 15),
            width=430,
        )
        self.evaluation_parameter = pn.widgets.Select(
            name='Evaluation parameter',
            value='No. of covered precursors',
            options=['No. of covered precursors'],
            width=430,
            margin=(15, 15, 0, 15)
        )
        self.optimize_button = pn.widgets.Button(
            name='Optimize',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )
        self.optimize_spinner = pn.indicators.LoadingSpinner(
            value=False,
            bgcolor='light',
            color='secondary',
            align='center',
            margin=(0, 0, 0, 20),
            width=35,
            height=35
        )
        self.scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF = pn.widgets.LiteralInput(
            name='Scan area A1/A2/B1/B2 (only used for specific diaPASEF)',
            value=[0,0,0,0],
            type=list,
            margin=(15, 15, 0, 15),
            width=900
        )

    def create(self):
        self.layout = pn.Card(
            pn.Row(
                pn.Column(
                    pn.WidgetBox(
                        pn.Row(
                            self.n_calls,
                            self.n_start,
                            sizing_mode='stretch_width',
                        ),
                        pn.Row(
                            self.initial_points,
                            self.evaluation_parameter,
                            sizing_mode='stretch_width'
                        ),
                        pn.Row(
                            self.YA1,
                            self.YA2,
                            sizing_mode='stretch_width'
                        ),
                        pn.Row(
                            self.YB1,
                            self.YB2,
                            sizing_mode='stretch_width'
                        ),
                        sizing_mode='stretch_width',
                        margin=(20, 10, 30, 10),
                        height=270
                    ),
                    margin=(10, 30, 10, 10),
                ),
                pn.Spacer(sizing_mode='stretch_width'),
                pn.Row(
                    self.optimize_button,
                    self.optimize_spinner,
                    # self.import_error,
                    align='start',
                    margin=(100, 10, 0, 0),
                )
            ),
            pn.layout.Divider(
                sizing_mode='stretch_width',
                margin=(-20, 10, -20, 10),
            ),
            pn.Column(
                None,
                pn.Row(
                    None,
                    None,
                ),
                align='center',
            ),
            title='Optimization',
            collapsed=False,
            collapsible=True,
            header_background='#eaeaea',
            background ='white',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            # height=470,
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

        dependances = {
            self.n_calls: [self.update_parameters, 'value'],
            self.n_start: [self.update_parameters, 'value'],
            self.initial_points: [self.update_parameters, 'value'],
            self.evaluation_parameter: [self.update_parameters, 'value'],
            self.YA1: [self.update_parameters, 'value'],
            self.YA2: [self.update_parameters, 'value'],
            self.YB1: [self.update_parameters, 'value'],
            self.YB2: [self.update_parameters, 'value'],
            self.optimize_button: [self.run_optimization, 'clicks'],
        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1],
                onlychanged=True
            )
        return self.layout

    def update_parameters(self, event):
        global method_conf
        convertion_dict = {
            self.n_calls.name: "n_calls",
            self.n_start.name: "n_start",
            self.initial_points: "initial_points",
            self.evaluation_parameter.name: "evaluation_parameter",
            self.YA1.name: "YA1",
            self.YA2.name: "YA2",
            self.YB1.name: "YB1",
            self.YB2.name: "YB2",
        }
        if isinstance(event.new, tuple):
            method_conf['optimizer'][convertion_dict[event.obj.name]] = list(event.new)
        else:
            method_conf['optimizer'][convertion_dict[event.obj.name]] = event.new

    def run_optimization(self, *args):
        self.optimize_spinner.value = True

        self.folder_path = [
            os.path.join(
                method_conf["input"]["save_at"],
                'optimization_plots'
            )
        ]
        diaid_pasef.main.create_folder(self.folder_path)

        self.opt_result = diaid_pasef.main.optimization(
            self.data.library,
            method_conf["method_parameters"],
            self.data.xi,
            self.data.yi,
            self.data.zi,
            method_conf,
            method_conf["optimizer"]
        )
        self.scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF.value = self.opt_result.x

        self.filenames_plots =  diaid_pasef.loader.get_file_names_from_directory(
            self.folder_path[0],
            'png'
        )

        self.player = pn.widgets.Player(
            start=0,
            end=len(self.filenames_plots)-1,
            value=0,
            loop_policy='loop',
            align='center',
            margin=(0, 0, 50, 0)
        )

        opt_plot_df = diaid_pasef.loader.create_opt_plot_df(
            self.filenames_plots[self.player.value]
        )

        self.kde_plot = pn.pane.PNG(
            object=os.path.join(
                self.folder_path[0],
                self.filenames_plots[self.player.value]
            ),
            height=345,
            width=460,
            align='center',
        )
        self.kde_plot_table = pn.widgets.DataFrame(
            opt_plot_df,
            autosize_mode='fit_viewport',
            margin=(0, 0, 0, 100),
            align='center',
            auto_edit=False
        )

        self.layout[2][0] = self.player
        self.layout[2][1][0] = self.kde_plot
        self.layout[2][1][1] = self.kde_plot_table

        self.player.param.watch(
            self.update_kde_plot_df,
            'value'
        )

        self.trigger_dependancy()
        self.optimize_spinner.value = False

    def update_kde_plot_df(self, event):
        self.kde_plot_table.value = diaid_pasef.loader.create_opt_plot_df(
            self.filenames_plots[event.new]
        )
        self.kde_plot.object = os.path.join(
            self.folder_path[0],
            self.filenames_plots[event.new]
        )
        self.layout[2][1][0] = self.kde_plot
        self.layout[2][1][1] = self.kde_plot_table


class CreateMethodCard(BaseWidget):
    # TODO: docstring
    def __init__(self, data, opt_widget):
        super().__init__(name="Method_creation")
        self.data = data
        self.opt_widget = opt_widget
        self.create_button = pn.widgets.Button(
            name='Create',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )
        self.path_method = pn.widgets.TextInput(
            name='Specify the path to the method file:',
            placeholder=method_path_placeholder,
            value='/Users/diaid_pasef/static/DIAParameterspy3TC.txt',
            width=900,
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )

    def create(self):
        self.layout = pn.Card(
            pn.Row(
                pn.Column(
                    pn.WidgetBox(
                        pn.Row(
                            self.opt_widget.scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF,
                            sizing_mode='stretch_width'
                        ),
                        sizing_mode='stretch_width',
                        margin=(20, 10, 30, 10),
                        height=90
                    ),
                    margin=(10, 30, 10, 10),
                ),
                pn.Spacer(sizing_mode='stretch_width'),
                pn.Row(
                    self.create_button,
                    # self.import_error,
                    align='center',
                    margin=(100, 130, 0, 0),
                )
            ),
            pn.layout.Divider(
                sizing_mode='stretch_width',
                margin=(-20, 10, -20, 10),
            ),
            pn.Row(
                None,
                None,
                align='center'
            ),
            title='Create Method',
            collapsed=False,
            collapsible=True,
            header_background='#eaeaea',
            background ='white',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            # height=470,
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

        dependances = {
            self.opt_widget.scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF: [self.update_parameters, 'value'],
            self.create_button: [self.create_method, 'clicks'],
        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1],
                onlychanged=True
            )
        return self.layout

    def update_parameters(self, event):
        global method_conf
        convertion_dict = {
            self.opt_widget.scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF.name: "scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF",
        }
        method_conf['method_parameters'][convertion_dict[event.obj.name]] = event.new

    def create_method(self, event):
        df_parameters_final = diaid_pasef.main.create_final_method(
            self.data.library,
            method_conf["method_parameters"],
            self.opt_widget.scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF.value,
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
        ).to_csv(
            os.path.join(
                method_conf["input"]["save_at"],
                'final_method',
                'input_parameters.csv'
            ),
            index=False
        )

        diaid_pasef.main.final_method_information(
            df_parameters_final,
            self.data.xi,
            self.data.yi,
            self.data.zi,
            method_conf,
            method_conf["input"]["save_at"],
            self.data.library,
            method_conf["method_parameters"],
            self.opt_widget.opt_result
        )
        final_method_path = os.path.join(
            method_conf["input"]["save_at"],
            'final_method',
            'diaPASEF_method.txt'
        )
        self.path_method.value = final_method_path
        self.layout[2][0] = pn.pane.PNG(
            object=os.path.join(
                method_conf["input"]["save_at"],
                'final_method',
                'Kernel_density_distribution_and_final_method.png'
            ),
            height=345,
            width=460,
            align='center',
        )
        self.layout[2][1] = pn.widgets.DataFrame(
            pd.read_csv(
                final_method_path,
                header=None,
                skiprows=3,
                names=["MS Type", "Cycle Id", "Start IM", "End IM", "Start Mass", "End Mass"]
            ),
            autosize_mode='fit_viewport',
            margin=(0, 0, 0, 100),
            align='center',
            auto_edit=False,
            height=240,
        )
        self.trigger_dependancy()

class EvaluateMethodCard(object):
    # TODO: docstring
    def __init__(self, method_creation):
        self.method_creation = method_creation
        self.evaluate_button = pn.widgets.Button(
            name='Evaluate',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )

    def create(self):
        self.layout = pn.Card(
            pn.Row(
                pn.Column(
                    self.method_creation.path_method,
                    sizing_mode='stretch_width',
                    margin=(20, 10, 30, 10),
                ),
                pn.Spacer(sizing_mode='stretch_width'),
                pn.Row(
                    self.evaluate_button,
                    align='center',
                    margin=(100, 130, 0, 0),
                )
            ),
            pn.layout.Divider(
                sizing_mode='stretch_width',
                margin=(-20, 10, -20, 10),
            ),
            pn.Row(
                None
            ),
            title='Method Evaluation',
            collapsed=False,
            collapsible=True,
            header_background='#eaeaea',
            background ='white',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            # height=470,
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

        dependances = {
            self.method_creation.path_method: [self.update_parameters, 'value'],
            # self.evaluate_button: [, 'clicks'],
        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1],
                onlychanged=True
            )
        return self.layout

    def update_parameters(self, event):
        global method_conf
        convertion_dict = {
            self.method_creation.path_method.name: "diaPASEF_method_only_used_for_method_evaluation",
        }
        method_conf['input'][convertion_dict[event.obj.name]] = event.new


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

        self.data = LoadLibraryCard()
        self.method_parameters = SpecifyParametersCard()
        self.optimization = OptimizationCard(self.data)
        self.method_creation = CreateMethodCard(self.data, self.optimization)
        self.method_evaluation = EvaluateMethodCard(self.method_creation)
        self.layout += [
            self.main_widget.create(),
            self.data.create(),
            self.method_parameters.create(),
            self.optimization.create(),
            self.method_creation.create(),
            self.method_evaluation.create(),
        ]
        if start_server:
            self.start_server()


def run():
    # set style
    init_panel()
    DiAIDPasefGUI(start_server=True)


if __name__ == '__main__':
    run()
