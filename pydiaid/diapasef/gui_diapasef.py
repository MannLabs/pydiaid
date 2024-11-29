# from bokeh.io import curdoc

#!python
import os
import json
import glob

# import logging
import platform
import pandas as pd

# visualization
import panel as pn
import imageio.v3 as iio

# modules
import pydiaid.diapasef as pydiaid
import pydiaid.loader as  loader
import pydiaid.diapasef.main as main
import pydiaid.diapasef.graphs as graphs
import pydiaid.diapasef.method_optimizer as method_optimizer
import pydiaid.diapasef.method_evaluation as method_evaluation
import pydiaid.synchropasef.plots as syP_plots

import warnings
warnings.filterwarnings("ignore")

# paths
BASE_PATH = os.path.dirname(__file__)
IMG_PATH = os.path.join(BASE_PATH, "img")
STYLE_PATH = os.path.join(BASE_PATH, "style")
DOCS_PATH = os.path.join(BASE_PATH, "docs")
DEFAULT_FILE = os.path.join(
    BASE_PATH,
    'static',
    "default_parameters.json"
)

with open(DEFAULT_FILE, "r") as infile:
    method_conf = json.load(infile)

if platform.system() == 'Windows':
    method_path_placeholder = 'D:\pydiaid\pydiaid\diapasef\static\diaPASEF_method.txt'
    library_path_placeholder = 'D:\pydiaid\pydiaid\diapasef\static\AlphaPept_results.csv'
    save_folder_placeholder = 'D:\pydiaid\pydiaid\diapasef\static'
else:
    method_path_placeholder = '/Users/pydiaid/diapasef/static/diaPASEF_method.txt'
    library_path_placeholder = '/Users/pydiaid/diapasef/static/AlphaPept_results.csv'
    save_folder_placeholder = '/Users/pydiaid/diapasef/static'


class BaseWidget(object):
    """
    BaseWidget class that initializes a base widget with given name and contains methods to track
    updates and trigger dependancies on update_event.

    Attributes:
        name (str): the name of the widget.
        __update_event (IntInput): an integer input widget that tracks the number of updates.
        depends (depends): dependent on the value of update_event, used for setting up callback functions.
        active_depends (depends): similar to depends, but also execute callbacks immediately upon 
            instantiation.
    Methods:
    __init__: The constructor for BaseWidget class.
    trigger_dependancy: Increases the value of __update_event by 1 to trigger the dependency.
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


class LoadLibraryCard(BaseWidget):
    """
    LoadLibraryCard class that extends BaseWidget class. The class represents a layout to load a 
    specified proteomics library, generates library related plots and analyze the presence of multiple 
    charged precursors.

    Attributes:
    library (None): The library to be loaded for evaluation.
    layout (None): The layout of the widget.
    path_library (pn.widgets.TextInput): Widget for entering the file path of the proteomics library.
    path_save_folder (pn.widgets.TextInput): Widget for entering the destination folder path for 
        saving output.
    ptm (pn.widgets.LiteralInput): Widget for specifying the Post-Translational Modification under analysis.
    analysis_software (pn.widgets.Select): Widget for selecting the analysis software.
    plot_mz (pn.widgets.EditableRangeSlider): Widget for defining the plot range for m/z.
    plot_im (pn.widgets.EditableRangeSlider): Widget for defining the ion mobility plot range.
    numbins (pn.widgets.IntInput): Widget for defining the number of bins for a kernel density plot.
    window_transparency (pn.widgets.FloatInput): Widget for defining the transparency of the dia window.
    window_frame_color (pn.widgets.Select): Widget for selecting the window frame color.
    window_color (pn.widgets.Select): Widget for selecting the window color.
    upload_button (pn.widgets.Button): Button to upload the library.
    upload_progress (pn.indicators.Progress): Indicator to show the progress of the data upload.

    Methods:
    init: Constructor for the LoadLibraryCard class.
    create_layout: Creates and returns the widget layout.
    update_parameters: Updates the method configuration parameters based on user input.
    update_parameters_plotting: Updates the plot configuration parameters based on user input.
    upload_data: Uploads the specified library, performs data analysis and generates library related 
        plots.
    """
    def __init__(self, description):
        super().__init__(name="Data")
        self.library = None
        self.layout = None
        self.path_library = pn.widgets.TextInput(
            name='Specify the path to the library:',
            placeholder=library_path_placeholder,
            value=os.path.join(os.path.dirname(__file__), "static", "AlphaPept_results.csv"),
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )
        self.path_save_folder = pn.widgets.TextInput(
            name='Save the output to the following folder:',
            # value=method_conf["input"]["save_at"],
            placeholder=save_folder_placeholder,
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )
        self.ptm = pn.widgets.LiteralInput(
            name='Specify the PTM:',
            # value=method_conf["input"]["PTM"],
            placeholder="['Phospho']",
            sizing_mode='stretch_width',
            margin=(15, 15, 15, 15)
        )
        self.analysis_software = pn.widgets.Select(
            name='Analysis software',
            value=method_conf["input"]["analysis_software"],
            options=['AlphaPept', 'MaxQuant', 'MSFragger', 'Spectronaut single-run', 'Spectronaut library', "DIANN single-run", "DIANN library", 'AlphaPeptDeep library', 'AlphaDIA'],
            min_width=200,
            margin=(15, 15, 15, 15)
        )
        # PLOTS
        self.plot_mz = pn.widgets.EditableRangeSlider(
            name='Plot m/z-range [Da]',
            start=100,
            end=1700,
            value=tuple(method_conf['graphs']['plot_mz']),
            step=50,
            margin=(15, 15, 0, 15),
        )
        self.plot_im = pn.widgets.EditableRangeSlider(
            name='Plot ion mobility range [1/K0]',
            start=0.6,
            end=1.6,
            value=tuple(method_conf['graphs']['plot_IM']),
            step=0.1,
            margin=(15, 15, 0, 15),
        )
        self.numbins = pn.widgets.IntInput(
            name='Number of bins',
            start=1,
            end=300,
            value=method_conf['graphs']['numbins'],
            step=1,
            margin=(15, 15, 0, 15),
        )
        self.window_transparency = pn.widgets.FloatInput(
            name='Transparency',
            start=0.1,
            end=1,
            value=method_conf['graphs']['window_transparency'],
            step=0.1,
            margin=(15, 15, 0, 15),
        )
        self.window_frame_color = pn.widgets.Select(
            name='Frame color',
            value=method_conf['graphs']['window_frame_color'],
            options=['black', 'grey'],
            margin=(15, 15, 0, 15),
        )
        self.window_color = pn.widgets.Select(
            name='Color',
            value=method_conf['graphs']['window_color'],
            options=['yellow', 'green', 'black', 'grey', 'white'],
            margin=(15, 15, 0, 15),
        )
        self.load_library_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            # css_classes=['main-part'],
            align='start',
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
            active=False,
            bar_color='light',
            width=250,
            align='center',
            margin=(0, 0, 30, 0)
        )
        self.import_error = pn.pane.Alert(
            alert_type="danger",
            # sizing_mode='stretch_width',
            # object='test warning message',
            margin=(30, 15, 5, 15),
        )

    def create_layout(self):
        plot_widgets = pn.Column(
            pn.Row(
                pn.Column(self.plot_mz), 
                pn.Column(self.plot_im)),
            pn.Row(self.numbins, self.window_transparency, sizing_mode='stretch_width'),
            pn.Row(self.window_frame_color, self.window_color, sizing_mode='stretch_width'),
            pn.Spacer(height=20), 
            sizing_mode='stretch_width',
        )

        self.layout = pn.Card(
            self.load_library_descr,
            pn.Row(
                pn.Column(
                    self.path_library,
                    self.path_save_folder,
                    pn.Row(self.analysis_software, self.ptm),
                    pn.WidgetBox(
                        plot_widgets,
                        sizing_mode='stretch_width',
                        margin=(20, 10, 50, 10),
                    ),
                    margin=(10, 30, 10, 10),
                    sizing_mode='stretch_width',
                ),
                pn.Column(
                    pn.Spacer(),  
                    self.upload_button,
                    self.upload_progress,
                    pn.Spacer(),  
                    width=400,
                    align='center',
                ),
                sizing_mode='stretch_width',
            ),
            pn.layout.Divider(
            ),
            pn.Row(
                None, None, None
            ),
            title='Load Library',
            collapsed=False,
            collapsible=True,
            header_background='#eaeaea',
            background ='white',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

        dependances = {
            self.path_library: [self.update_parameters, 'value'],
            self.path_save_folder: [self.update_parameters, 'value'],
            self.ptm: [self.update_parameters, 'value'],
            self.analysis_software: [self.update_parameters, 'value'],
            self.plot_mz: [self.update_parameters_plotting, 'value'],
            self.plot_im: [self.update_parameters_plotting, 'value'],
            self.numbins: [self.update_parameters_plotting, 'value'],
            self.window_transparency: [self.update_parameters_plotting, 'value'],
            self.window_frame_color: [self.update_parameters_plotting, 'value'],
            self.window_color: [self.update_parameters_plotting, 'value'],
            self.upload_button: [self.upload_data, 'clicks'],

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
            self.path_library.name: "library_name",
            self.path_save_folder.name: "save_at",
            self.ptm.name: "PTM",
            self.analysis_software.name: "analysis_software",
        }
        method_conf['input'][convertion_dict[event.obj.name]] = event.new

    def update_parameters_plotting(self, event):
        global method_conf
        convertion_dict = {
            self.plot_mz.name: "plot_mz",
            self.plot_im.name: "plot_IM",
            self.numbins.name: "numbins",
            self.window_frame_color.name: "window_frame_color",
            self.window_color.name: "window_color",
            self.window_transparency.name: "window_transparency",
        }
        method_conf['graphs'][convertion_dict[event.obj.name]] = event.new

    def upload_data(self, *args):
        self.upload_progress.active = True
        self.library = loader.load_library(
            self.path_library.value,
            method_conf["input"]["analysis_software"],
            method_conf["input"]["PTM"]
        )

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
        main.create_folder(folder_paths)
        self.xi, self.yi, self.zi = graphs.kernel_density_calculation(
            self.library,
            method_conf["graphs"]["numbins"]
        )
        self.layout[3][0] = pn.Column(
            pn.pane.Markdown(
                '### Histogram: Precursor distribution in m/z',
                 align='center'
            ),
            pn.pane.Matplotlib(
            graphs.plot_precursor_distribution_as_histogram(
                self.library,
                method_conf["graphs"],
                os.path.join(
                    self.path_save_folder.value,
                    'input_library',
                    'Histogram_precursor_distribution_in_library.png'
                ),
                gui=True
            ),
            # margin=(20, 0),
            tight=True
        )
        )
        self.layout[3][1] = pn.Column(
            pn.pane.Markdown(
                '### Precursor cloud plotted across m/z and ion mobility',
                 align='center'
            ),
            pn.pane.Matplotlib(
            graphs.plot_density(
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
            tight=True,
        )
        )
        dict_charge_of_precursor = method_evaluation.calculate_percentage_multiple_charged_precursors(self.library)
        mult_charged_precursor_info = pd.DataFrame(
            {
                "charge state of precursors": list(dict_charge_of_precursor.keys()),
                "ratio [%]": list(dict_charge_of_precursor.values())
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
        self.layout[3][2] = pn.Column(
            pn.pane.Markdown(
                '### Percentage of multiple charged precursors',
                 align='center'
            ),
            pn.widgets.Tabulator(
                mult_charged_precursor_info,
                layout='fit_data_table', 
                width=400
            ),
            sizing_mode='stretch_width',
        )
        self.trigger_dependancy()
        self.upload_progress.active = False


class SpecifyParametersCard(BaseWidget):
    """Specifies a Parameters Card Class.

    This class, which inherits BaseWidget class, initializes a widget with name 'Parameters'. It also 
    contains methods to create the layout of the widget, update parameters and run calculation for the
    precursor content within the scan area.

    Attributes:
    data (str): The input data.
    mz (EditableRangeSlider): Slider for selecting mass to charge ratio range.
    ion_mobility (EditableRangeSlider): Slider for selecting ion mobility range
    num_dia_pasef_scans (IntInput): Selects the number of dia-PASEF scans, ranging from 1 to 40.
    im_steps (IntInput): Selects the number of ion mobility windows per dia-PASEF scan, ranging from 
        1 to 20.
    overlap (IntInput): Selects the isolation window overlap value in Dalton, ranging from 0 to 10.
    shift_of_final_method (FloatInput): Adjusts the shift of the final acquisition scheme (in the
        ion mobility dimension), ranging from -0.5 to 1.0.
    max_width (IntInput): Sets maximum allowed window width [Da].
    spec_param_table (DataFrame): Specificies the parameters table within the Widget.
    specify_parameter_descr (Markdown): The description for specifying parameters.
    calculate_button (Button): A button widget to initiate calculation of parameters.

    Methods:
    init: The constructor for SpecifyParametersCard class.
    create_layout: Creates the layout for the Parameters widget.
    update_parameters: Updates the global method configuration with new parameters specified by users.
    run_calculation: Performs calculations based on the selected parameters and updates the layout with 
        the result.
    """
    def __init__(self, data, description):
        super().__init__(name="Parameters")
        self.data = data
        self.mz = pn.widgets.EditableRangeSlider(
            name='m/z-range [Da]',
            start=100,
            end=1700,
            value=tuple(method_conf['method_parameters']['mz']),
            step=50,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.ion_mobility = pn.widgets.EditableRangeSlider(
            name='Ion mobility range [1/K0]',
            start=0.6,
            end=1.6,
            value=tuple(method_conf['method_parameters']['ion_mobility']),
            step=0.05,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.num_dia_pasef_scans = pn.widgets.IntInput(
            name='Number of dia-PASEF scans',
            start=1,
            end=40,
            value=method_conf['method_parameters']['num_dia_pasef_scans'],
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.im_steps = pn.widgets.IntInput(
            name='Number of ion mobility windows / dia-PASEF scan',
            start=1,
            end=20,
            value=method_conf['method_parameters']['im_steps'],
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.overlap = pn.widgets.IntInput(
            name='Isolation window overlap [Da]',
            start=0,
            end=10,
            value=method_conf['method_parameters']['overlap'],
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.shift_of_final_method = pn.widgets.FloatInput(
            name='Shift of the final acquisition scheme (in IM dimension) [1/K0]',
            start=-0.5,
            end=1.0,
            value=method_conf['method_parameters']['shift_of_final_method'],
            step=0.01,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.spec_param_table = pn.widgets.DataFrame(
            autosize_mode='fit_viewport',
            align='center',
            auto_edit=False,
        )
        self.specify_parameter_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start',
        )
        self.calculate_button = pn.widgets.Button(
            name='Calculate',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )

    def create_layout(self):
        parameter_widgets = pn.Column(
            pn.Row(pn.Column(self.mz), pn.Column(self.ion_mobility)),
            pn.Row(self.num_dia_pasef_scans, self.im_steps, sizing_mode='stretch_width'),
            pn.Row(self.overlap, self.shift_of_final_method, sizing_mode='stretch_width'),
            pn.Spacer(height=20),
            sizing_mode='stretch_width',
        )

        self.layout = pn.Card(
            self.specify_parameter_descr,
            pn.Row(
                pn.Column(
                    pn.WidgetBox(
                        parameter_widgets,
                        sizing_mode='stretch_width',
                        margin=(20, 10, 50, 10),
                    ),
                    margin=(10, 30, 10, 10),
                    sizing_mode='stretch_width',
                ),
                pn.Column(
                    pn.Spacer(),
                    self.calculate_button,
                    pn.Spacer(),
                    width=400,
                    align='center',
                ),
                sizing_mode='stretch_width',
            ),
            pn.layout.Divider(),
            pn.Row(
                None,
            ),
            title='Specify Method Parameters',
            collapsed=False,
            collapsible=True,
            header_background='#eaeaea',
            background='white',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
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
        dict_precursors_within_scan_area = method_evaluation.calculate_precursor_within_scan_area(
            self.data.library,
            method_conf['method_parameters']["mz"],
            method_conf['method_parameters']["ion_mobility"]
        )

        df_precursors_within_scan_area = pd.DataFrame(
            {
                "precursors within the scan area [%]": list(dict_precursors_within_scan_area.values())
            }
        )

        df_precursors_within_scan_area.to_csv(
            os.path.join(
                method_conf["input"]["save_at"],
                'final_method',
                'precursors_within_scan_area.csv'
            ),
            index=False
        )
        
        self.layout[3][0] = pn.Column(
            pn.pane.Markdown(
                '### Precursors within scan area',
                align='center'
            ),
            pn.widgets.Tabulator(
                df_precursors_within_scan_area,
                layout='fit_data_table', 
                width=400
            ),
            margin=(20, 50),
            sizing_mode='stretch_width',
        )



class OptimizationCard(BaseWidget):
    """Optimization card class that inherits BaseWidget.

    This class initializes the Optimization card widget with methods to manipulate the layout, update 
    parameters and run the optimization process for dia-PASEF methods.

    Attributes:
    data (str): Input data provided.
    n_calls (IntInput): An integer input widget field to accept the number of iteration optimization 
        steps. Recommended 100.
    initial_points (IntInput): An integer input widget field to accept the number of initial points for
        the optimization. Recommended 20.
    YA1 & YA2 (EditableRangeSlider): Editable range sliders to select the 'A1' and 'A2' range respectively.
    YB1 & YB2 (EditableRangeSlider): Editable range sliders to select the 'B1' and 'B2' range respectively.
    evaluation_parameter (Select): A select widget field to accept the evaluation parameter.
    optimize_button (Button): A button that initiates the optimization process.
    optimize_spinner (LoadingSpinner): A loading spinner widget that indicates optimization process.
    scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF (LiteralInput): Literal input field to define a 
        specific scan area for a dia-PASEF method.
    optimization_descr (Markdown): The description for the “Optimization” widget.

    Methods:
    init: The constructor for OptimizationCard class.
    create_layout: Creates the layout for the “Optimization” widget.
    update_parameters: Updates the global method configuration with new parameters specified by users.
    run_optimization: Activates the optimization process based on the selected parameters.
    update_kde_plot_df: This updates the Kernel Density Estimation (KDE) plot DataFrame on change in 
        player value.
    """
    def __init__(self, data, description):
        super().__init__(name="Optimize Method")
        self.data = data
        self.opt_result_x = [0, 0, 0, 0]
        self.n_calls = pn.widgets.IntInput(
            name='Number of iterative optimization steps',
            start=1,
            end=200,
            value=method_conf['optimizer']['n_calls'],
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.initial_points = pn.widgets.IntInput(
            name='Number of starting points',
            start=1,
            end=20,
            value=method_conf['optimizer']['initial_points'],
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.YA1 = pn.widgets.EditableRangeSlider(
            name='A1 range',
            start=0.2,
            end=2.0,
            value=tuple(method_conf['optimizer']['YA1']),
            step=0.1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.YA2 = pn.widgets.EditableRangeSlider(
            name='A2 range',
            start=0.2,
            end=2.0,
            value=tuple(method_conf['optimizer']['YA2']),
            step=0.05,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.YB1 = pn.widgets.EditableRangeSlider(
            name='B1 range',
            start=0.2,
            end=2.0,
            value=tuple(method_conf['optimizer']['YB1']),
            step=0.1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.YB2 = pn.widgets.EditableRangeSlider(
            name='B2 range',
            start=0.2,
            end=2.0,
            value=tuple(method_conf['optimizer']['YB2']),
            step=0.05,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.evaluation_parameter = pn.widgets.Select(
            name='Evaluation parameter',
            value=method_conf['optimizer']['evaluation_parameter'],
            options=[
                "No. of covered proteins",
                'No. of covered precursors',
                "No. of covered, doubly charged precursors",
                "No. of covered, triply charged precursors",
                "No. of covered, quadruply charged precursors"
            ],
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.optimize_button = pn.widgets.Button(
            name='Optimize',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )
        self.optimize_progress = pn.indicators.Progress(
            active=False,
            bar_color='light',
            width=250,
            align='center',
            margin=(0, 0, 30, 0)
        )
        self.scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF = pn.widgets.LiteralInput(
            name='Scan area A1/A2/B1/B2',
            value=[0,0,0,0],
            type=list,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.optimization_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start',
        )

    def create_layout(self):
        optimization_widgets = pn.Column(
            pn.Row(self.n_calls, self.initial_points, sizing_mode='stretch_width'),
            pn.Row(self.evaluation_parameter, sizing_mode='stretch_width'),
            pn.Row(self.YA1, self.YA2, sizing_mode='stretch_width'),
            pn.Row(self.YB1, self.YB2, sizing_mode='stretch_width'),
            pn.Spacer(height=20),
            sizing_mode='stretch_width',
        )

        self.layout = pn.Card(
            self.optimization_descr,
            pn.Row(
                pn.Column(
                    pn.WidgetBox(
                        optimization_widgets,
                        sizing_mode='stretch_width',
                        margin=(20, 10, 50, 10),
                    ),
                    margin=(10, 30, 10, 10),
                    sizing_mode='stretch_width',
                ),
                pn.Column(
                    pn.Spacer(),
                    self.optimize_button,
                    self.optimize_progress,
                    pn.Spacer(),
                    width=400,
                    align='center',
                ),
                sizing_mode='stretch_width',
            ),
            pn.layout.Divider(),
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
            background='white',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

        dependances = {
            self.n_calls: [self.update_parameters, 'value'],
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
            self.initial_points.name: "initial_points",
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
        self.optimize_progress.active = True

        self.folder_path = [
            os.path.join(
                method_conf["input"]["save_at"],
                'optimization_plots'
            )
        ]
        main.create_folder(self.folder_path)

        self.opt_result = method_optimizer.optimization(
            self.data.library,
            method_conf["method_parameters"],
            self.data.xi,
            self.data.yi,
            self.data.zi,
            method_conf,
            method_conf["optimizer"]
        )

        self.scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF.value = self.opt_result

        self.filenames_plots = loader.get_file_names_from_directory(
            self.folder_path[0],
            'png'
        )

        self.start_position_player = method_optimizer.find_matching_filename_index(
            self.filenames_plots, 
            self.scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF.value
            )

        self.player = pn.widgets.Player(
            start=0,
            end=len(self.filenames_plots)-1,
            value=self.start_position_player,
            loop_policy='loop',
            align='center',
            margin=(20, 0)
        )

        opt_plot_df = loader.create_opt_plot_df(
            self.filenames_plots[self.player.value]
        )

        self.kde_plot = pn.pane.PNG(
            object=os.path.join(
                self.folder_path[0],
                self.filenames_plots[self.player.value]
            ),
            # height=345,
            # width=460,
            align='center',
        )
        self.kde_plot_table = pn.widgets.Tabulator(
            opt_plot_df,
            margin=(0, 0, 20, 100),
            layout='fit_data_table', 
            width=350,
        )

        self.layout[3][0] = pn.Column(
            pn.pane.Markdown(
                '### Optimization Results',
                align='center'
            ),
            self.player,
            sizing_mode='stretch_width',
        )
        self.layout[3][1][0] = self.kde_plot
        self.layout[3][1][1] = self.kde_plot_table

        self.player.param.watch(
            self.update_kde_plot_df,
            'value'
        )

        self.trigger_dependancy()
        self.optimize_progress.active = False

    def update_kde_plot_df(self, event):
        self.kde_plot_table.value = loader.create_opt_plot_df(
            self.filenames_plots[event.new]
        )
        self.kde_plot.object = os.path.join(
            self.folder_path[0],
            self.filenames_plots[event.new]
        )
        self.layout[3][1][0] = self.kde_plot
        self.layout[3][1][1] = self.kde_plot_table


class CreateMethodCard(BaseWidget):
    """CreateMethod card that inherits the BaseWidget class.

    This class initializes the “Create Method” widget and provides methods to create the layout, 
    update the parameters and create a dia-PASEF method.

    Attributes:
    data (str): The input data.
    opt_widget (BaseWidget): An optimal widget inherited from the BaseWidget class.
    create_button (Button): Button widget to initiate the method creation process.
    path_method (TextInput): Input field for the path to the method file.
    create_method_descr (Markdown): The description for the “Create Method” widget.

    Methods:
    init: The constructor for CreateMethodCard class.
    create_layout: Arranges the structure layout of the “Create Method” widget.
    update_parameters: Updates the global method configuration according to user's specifications.
    create_method: Runs the method creation process and provides the Tabulator data frame of method 
        details.
    """
    def __init__(self, data, opt_widget, description):
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
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )
        self.create_method_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start',
        )

    def create_layout(self):
        create_method_widgets = pn.Column(
            pn.Row(
                self.opt_widget.scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF,
            ),
            pn.Spacer(height=20),
        )

        self.layout = pn.Card(
            self.create_method_descr,
            pn.Row(
                pn.Column(
                    pn.WidgetBox(
                        create_method_widgets,
                        sizing_mode='stretch_width',
                        margin=(20, 10, 50, 10),
                    ),
                    margin=(10, 30, 10, 10),
                    sizing_mode='stretch_width',
                ),
                pn.Column(
                    pn.Spacer(),
                    self.create_button,
                    pn.Spacer(),
                    width=400,
                    align='center',
                ),
                sizing_mode='stretch_width',
            ),
            pn.layout.Divider(),
            pn.Column(
                None,
                align='center',
                sizing_mode='stretch_width',
            ),
            title='Create Method',
            collapsed=False,
            collapsible=True,
            header_background='#eaeaea',
            background='white',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
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
            self.opt_widget.scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF.name: "scan_area_A1_A2_B1_B2_only_used_for_create",
        }
        method_conf['method_parameters'][convertion_dict[event.obj.name]] = event.new

    def create_method(self, event):
        df_parameters_final = main.create_final_method(
            self.data.library,
            method_conf["method_parameters"],
            self.opt_widget.scan_area_A1_A2_B1_B2_only_used_for_specific_diaPASEF.value,
            method_conf
        )

        # save input parameters as json file for reusing
        out_file = open(
                method_conf["input"]["save_at"]+'/input_parameters.json', "w")
        json.dump(method_conf, out_file, indent = 4)
        out_file.close()

        final_method_path = os.path.join(
            method_conf["input"]["save_at"],
            'final_method',
            'diaPASEF_method.txt'
        )
        self.path_method.value = final_method_path
        df_method = pd.read_csv(
                final_method_path,
            header=None,
            skiprows=3,
            names=["MS Type", "Cycle Id", "Start IM", "End IM", "Start Mass", "End Mass", "CE"]
        )
        self.layout[3] = pn.Column(
            pn.pane.Markdown(
                '### Created Method Parameters',
                align='center'
            ),
            pn.widgets.Tabulator(
                df_method,
                margin=(20, 0, 20, 100),
                layout='fit_data_table', 
                sizing_mode='stretch_width',
            ),
            sizing_mode='stretch_width',
            align='center',
        )
        self.trigger_dependancy()


class EvaluateMethodCard(object):
    """
    EvaluateMethodCard class which helps to evaluate a dia-PASEF method. It calculates the precursor
    coverage and displays the method on top of a kernel density plot.

    Attributes:
    data (any): This represents the data being used for evaluating the methods.
    method_creation (any): It represents the methods being used for creation.
    evaluate_button (pn.widgets.Button): It is a Primary button for evaluating the methods.
    evaluate_method_descr (pn.pane.Markdown): A portion inside the GUI, where the description of the 
        evaluation method is being posted in Markdown format.

    Methods:
    init(self, data, method_creation, description):The constructor for CreateMethodCard class.
    create_layout(self): This method creates the GUI layout for the EvaluateMethodCard. It sets up 
        a card with description, method creations path, button for evaluation and sets up listeners 
        for update events. 
    update_parameters(self, event): 
        This method updates the global 'method_conf' dictionary with new event values. The events are 
        only triggered when changes are observed in the parameters.
    evaluate_method(self, event): 
        This method is triggered after clicking the evaluate button. It will evaluate the method, 
        calculates the precursor_within_scan_area and coverage, and then displays a DataFrame with the 
        evaluation results and a plot image where the method is plotted on top of a kernel density 
        distribution of peptides.
    """
    def __init__(self, data, method_creation, description):
        self.data = data
        self.method_creation = method_creation
        self.evaluate_button = pn.widgets.Button(
            name='Evaluate',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )
        self.evaluate_method_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start',
        )

    def create_layout(self):
        evaluate_method_widgets = pn.Column(
            pn.Row(
                self.method_creation.path_method,
                sizing_mode='stretch_width'
            ),
            pn.Spacer(height=20),
            sizing_mode='stretch_width',
        )

        self.layout = pn.Card(
            self.evaluate_method_descr,
            pn.Row(
                pn.Column(
                    pn.WidgetBox(
                        evaluate_method_widgets,
                        sizing_mode='stretch_width',
                        margin=(20, 10, 50, 10),
                    ),
                    margin=(10, 30, 10, 10),
                    sizing_mode='stretch_width',
                ),
                pn.Column(
                    pn.Spacer(),
                    self.evaluate_button,
                    pn.Spacer(),
                    width=400,
                    align='center',
                ),
                sizing_mode='stretch_width',
            ),
            pn.layout.Divider(),
            pn.Column(
                None,
                align='center',
                sizing_mode='stretch_width',
            ),
            title='Evaluate Method',
            collapsed=False,
            collapsible=True,
            header_background='#eaeaea',
            background='white',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

        dependances = {
            self.method_creation.path_method: [self.update_parameters, 'value'],
            self.evaluate_button: [self.evaluate_method, 'clicks'],
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

    def evaluate_method(self, event):

        method_eval_dir = os.path.dirname(self.method_creation.path_method.value)

        df_parameters_final = pd.read_csv(
            method_conf["input"]["diaPASEF_method_only_used_for_method_evaluation"],
            skiprows=4,
            names=["MS Type", "Cycle Id", "Start IM", "End IM", "Start Mass", "End Mass", "CE"]
        )
        self.df_method = df_parameters_final

        main.final_method_information(
            df_parameters_final,
            self.data.xi,
            self.data.yi,
            self.data.zi,
            method_conf,
            method_conf["input"]["save_at"],
            self.data.library,
            method_conf["method_parameters"],
            method_conf["method_parameters"]["scan_area_A1_A2_B1_B2_only_used_for_create"]
        )

        dict_precursors_within_mz = method_evaluation.calculate_precursor_within_scan_area(
            self.data.library,
            method_conf["method_parameters"]["mz"],
            method_conf["method_parameters"]["ion_mobility"]
        )
        dict_precursors_coverage = method_evaluation.coverage(
            df_parameters_final,
            self.data.library,
        )
        dict_evaluation_of_final_method = {
            **dict_precursors_within_mz,
            **dict_precursors_coverage
        }

        if method_conf["method_parameters"]["scan_area_A1_A2_B1_B2_only_used_for_create"][0] != int:
            pass
        else:
            dict_evaluation_of_final_method["final A1, A2, B1, B2 values"] = str([
                method_conf["method_parameters"]["scan_area_A1_A2_B1_B2_only_used_for_create"][0] + method_conf['method_parameters']["shift_of_final_method"],
                method_conf["method_parameters"]["scan_area_A1_A2_B1_B2_only_used_for_create"][1] + method_conf['method_parameters']["shift_of_final_method"],
                method_conf["method_parameters"]["scan_area_A1_A2_B1_B2_only_used_for_create"][2] + method_conf['method_parameters']["shift_of_final_method"],
                method_conf["method_parameters"]["scan_area_A1_A2_B1_B2_only_used_for_create"][3] + method_conf['method_parameters']["shift_of_final_method"]]
            )
        final_df = pd.DataFrame({
            "evaluation parameter": list(dict_evaluation_of_final_method.keys()),
            "value": list(dict_evaluation_of_final_method.values())
        })
        final_df.to_csv(
            os.path.join(
                method_conf["input"]["save_at"],
                'final_method',
                'parameters_to_evaluate_method.csv'
            ),
            index=False
        )

        self.layout[3] = pn.Column(
            pn.Row(
                pn.Column(
                    pn.pane.Markdown(
                        '### Final method plotted on top of precursor cloud',
                        align='center'
                    ),
                    pn.pane.PNG(
                        object=os.path.join(
                            method_conf["input"]["save_at"],
                            'final_method',
                            'Kernel_density_distribution_and_final_method.png'
                        ),
                        # height=345,
                        # width=460,
                        align='center',
                        margin=(50, 0, 0, 0)
                    ),
                    sizing_mode='stretch_width',
                ),
                pn.Column(
                    pn.pane.Markdown(
                        '### Evaluation Results',
                        align='center'
                    ),
                    pn.widgets.Tabulator(
                        final_df,
                        margin=(20, 0, 20, 100),
                        sizing_mode='stretch_width',
                    ),
                    sizing_mode='stretch_width',
                ),
                sizing_mode='stretch_width',
            ),
            sizing_mode='stretch_width',
            align='center',
        )


class MakeGifCard(object):
    """
    MakeGifCard class includes properties and methods that provide users the ability to create a GIF of 
    a dia-PASEF method.

    Attributes:
    layout: Placeholder for the layout of the widget card.
    data: Provided data for the card.
    load_method: Details of the method loaded
    make_gif_descr (Markdown): A pane that displays the description.
    gif_button (Button): A button widget that allows users to create a GIF of dia-PASEF method.
    upload_progress (Progress): An indicator that displays a progress bar during GIF upload.
    gif_duration (FloatInput): Widget that allows user to input the time per frame in ms.
    im_steps (IntInput): Widget that allows user to input the ion mobility steps.
    scans_plotted_at_once (IntInput): Widget that allows user to input the windows highlighted at once.

    Methods:
    init: Constructor for the MakeGifCard class.
    create_layout(): Creates the panel layout for this widget card and sets up dependencies.
    run_gif(*args): Creates a GIF of a dia-PASEF method and displays it in the user interface.
    """
    def __init__(self, data, load_method, description):
        self.data = data
        self.load_method = load_method
        self.make_gif_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start',
        )
        self.gif_button = pn.widgets.Button(
            name='Create dia-PASEF method GIF',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )
        self.gif_progress = pn.indicators.Progress(
            active=False,
            bar_color='light',
            width=250,
            align='center',
            margin=(0, 0, 30, 0)
        )
        self.gif_duration = pn.widgets.FloatInput(
            name='Time per frame [ms]',
            start=0.001,
            end=1,
            value=0.001,  # Default value adjusted for dia-PASEF
            step=0.001,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.im_steps = pn.widgets.IntInput(
            name='Ion mobility steps',
            start=0,
            end=200,
            value=10,  # Default value adjusted for dia-PASEF
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.scans_plotted_at_once = pn.widgets.IntInput(
            name='Windows highlighted at once',
            start=0,
            end=20,
            value=1,
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )

    def create_layout(self):
        gif_settings = pn.Column(
            pn.Row(self.gif_duration, self.scans_plotted_at_once, sizing_mode='stretch_width'),
            pn.Row(self.im_steps, sizing_mode='stretch_width'),
            pn.Spacer(height=20),
            sizing_mode='stretch_width',
        )

        self.layout = pn.Card(
            self.make_gif_descr,
            pn.Row(
                pn.Column(
                    pn.WidgetBox(
                        gif_settings,
                        sizing_mode='stretch_width',
                        margin=(20, 10, 50, 10),
                    ),
                    margin=(10, 30, 10, 10),
                    sizing_mode='stretch_width',
                ),
                pn.Column(
                    pn.Spacer(),
                    self.gif_button,
                    self.gif_progress,
                    pn.Spacer(),
                    width=400,
                    align='center',
                ),
                sizing_mode='stretch_width',
            ),
            pn.layout.Divider(),
            pn.Column(
                None,
                align='center',
                sizing_mode='stretch_width',
            ),
            title='Make a GIF of a dia-PASEF Method',
            collapsed=True,
            collapsible=True,
            header_background='#eaeaea',
            background='white',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

        dependances = {
            self.gif_button: [self.run_gif, 'clicks'],
        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1],
                onlychanged=True
            )

        return self.layout

    def run_gif(self, *args):
        self.gif_progress.active = True

        save_at = os.path.join(
            self.data.path_save_folder.value,
            'gif',
        )

        # Clean up existing files
        if os.path.exists(save_at):
            files = glob.glob(os.path.join(save_at, "*"))
            for f in files:
                os.remove(f)
        else:
            os.makedirs(save_at)

        # Generate boxes for dia-PASEF
        boxes, df_filtered = graphs.generate_dia_boxes(
            self.load_method.df_method, 
            im_steps=self.im_steps.value
        )
        # Generate GIF frames
        syP_plots.generate_gif_single_windows(
            self.data.xi,
            self.data.yi,
            self.data.zi,
            method_conf["graphs"],
            boxes,
            range(1, len(df_filtered)+1),
            save_at+ "/",
            facecolor="#FF0098",
            window_color="white_grey",
            scans_plotted_at_once=self.scans_plotted_at_once.value,
        )

        # Compile frames into GIF
        images = []
        for file_name in sorted(os.listdir(save_at)):
            if file_name.endswith('.png'):
                file_path = os.path.join(save_at, file_name)
                images.append(iio.imread(file_path))

        gif_path = os.path.join(save_at, "diaPASEF.gif")
        iio.imwrite(gif_path, images, duration=self.gif_duration.value, loop=0)

        # Display the GIF in the interface
        self.layout[3] = pn.Column(
            pn.pane.Markdown('### dia-PASEF Method GIF', align='center'),
            pn.pane.GIF(
                gif_path,
                alt_text="Sorry, please try using a different browser",
                width=640, height=360
            ),
            sizing_mode='stretch_width',
            align='center',
        )

        self.gif_progress.active = False


class DiaCardWidget(object):
    """
    DiaCardWidget class describes the widgets and layout used in a Dashboard Interface for assessing 
    dia-PASEF methods. The class involves detailed description for each section including library 
    loading, parameter specification, method optimization, dia-PASEF method creation, and method 
    evaluation.

    Attributes:
    project_description (str): A string description about the project's functionality.
    load_library_description (str): A string description about loading the proteomics library.
    specify_parameter_description (str): A string description specifying the method parameters.
    optimization_description (str): A string description on how dia-PASEF method optimization is done.
    create_method_description (str): A string description for creating the dia-PASEF method.
    evaluate_method_description (str): A string description on evaluating the dia-PASEF method.
    manual_path (str): A string path to the documentation.
    error_message_upload (str): An error message displayed when file upload fails.
    data (LoadLibraryCard): An object representing the library loading component.
    method_parameters (SpecifyParametersCard): An object representing the component where parameters 
        are specified.
    optimization (OptimizationCard): An object representing the method optimization component.
    method_creation (CreateMethodCard): An object representing the method creation component.
    method_evaluation (EvaluateMethodCard): An object representing the method evaluation component.

    Methods:
    init: The constructor for DiaCardWidget class.
    create_layout: Creates the layout for the widget by combining all the given components.
    """
    def __init__(self):
        self.project_description = """#### py_diAID uses an Automated Isolation Design to generate optimal dia-PASEF methods with respect to the precursor density. It designs isolation windows with variable widths, which enable short acquisition cycles, while essentially covering the complete m/z-ion mobility-range."""
        self.load_library_description = "#### Load the library for the indicated analysis software to check the distribution of the precursors in the m/z-ion mobility plane."
        self.specify_parameter_description = "####  We found a strong correlation between a high theoretical and empirical precursor coverage. This result suggests using a scan area with a wide m/z-range and a narrow ion mobility range. Specify the number of dia-PASEF scans, which depend on the chromatographic peak width, and the number of ion mobility windows per dia-PASEF scan. We recommend two ion mobility windows per dia-PASEF scan."
        self.optimization_description = "#### py_diAID uses a Bayesian optimization following a Gaussian process to find the optimal scan area. We recommend 100 optimization steps and 20 starting points."
        self.create_method_description = "#### Create a dia-PASEF method with an optimal or an individually specified scan area."
        self.evaluate_method_description = "#### Evaluate the optimal dia-PASEF method or confirm if an already existing dia-PASEF method is suitable for your experiment."
        self.make_gif_description = "#### Visualize the scanning of your dia-PASEF method with a GIF."
        self.manual_path = os.path.join(
            DOCS_PATH,
            "manual.pdf"
        )

        # ERROR/WARNING MESSAGES
        self.error_message_upload = "The selected file can't be uploaded. Please check the instructions for data uploading."

        self.data = LoadLibraryCard(self.load_library_description)
        self.method_parameters = SpecifyParametersCard(
            self.data,
            self.specify_parameter_description
        )
        self.optimization = OptimizationCard(
            self.data,
            self.optimization_description
        )
        self.method_creation = CreateMethodCard(
            self.data,
            self.optimization,
            self.create_method_description
        )
        self.method_evaluation = EvaluateMethodCard(
            self.data,
            self.method_creation,
            self.evaluate_method_description
        )
        self.make_gif_method = MakeGifCard(
            self.data,
            self.method_evaluation,
            self.make_gif_description
            )

    def create_layout(
        self,
        tab_list=None
    ):
        return pn.Column(
            self.data.create_layout(),
            self.method_parameters.create_layout(),
            self.optimization.create_layout(),
            self.method_creation.create_layout(),
            self.method_evaluation.create_layout(),
            self.make_gif_method.create_layout(),
            sizing_mode='stretch_width',
        )
