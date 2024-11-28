import pandas as pd
import os
import platform
import json

# import alphatims.bruker
# import alphatims.utils

# visualization
import panel as pn

import imageio.v3 as iio
import glob

#local
import pydiaid.synchropasef as pydiaid
import pydiaid.synchropasef.loader_MS_parameter_file as loader_MS_parameter_file
import pydiaid.loader as loader
import pydiaid.synchropasef.method_creator as method_creator
import pydiaid.synchropasef.method_evaluator as method_evaluator
import pydiaid.synchropasef.plots as plots


# paths
BASE_PATH = os.path.dirname(__file__)
STYLE_PATH = os.path.join(BASE_PATH, "style")
DOCS_PATH = os.path.join(BASE_PATH, "docs")
IMG_PATH = os.path.join(BASE_PATH, "img")
LATEST_GITHUB_INIT_FILE = "https://github.com/MannLabs/pydiaid/tree/main/pydiaid/__init__.py"

DEFAULT_FILE = os.path.join(
    BASE_PATH,
    'static',
    "default_parameters.json"
)

with open(DEFAULT_FILE, "r") as infile:
    method_conf = json.load(infile)

if platform.system() == 'Windows':
    method_path_placeholder = 'D:\pydiaid\pydiaid\synchropasef\static\synchroPASEF_method.txt'
    library_path_placeholder = 'D:\pydiaid\pydiaid\static\evidence_MaxQuant_270223.txt'
    save_folder_placeholder = 'D:\pydiaid\pydiaid\static'
else:
    method_path_placeholder = '/Users/pydiaid/static/synchroPASEF_method.txt'
    library_path_placeholder = '/Users/pydiaid/static/evidence_MaxQuant_270223.txt'
    save_folder_placeholder = '/Users/pydiaid/static'


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
    `LoadLibraryCard` class extends `BaseWidget` with properties and methods for loading and visualizing 
    proteomics libraries.

    This class provides graphical user interface components to input parameters for loading a proteomics 
    library, including file path, analysis software and post-translational modifications (PTM), and 
    provides components for visualizing loaded data.

    Attributes:
        description (str): Textual description shown in the UI.
        library: Expected to be a DataFrame to hold the loaded library data.
        layout (Card): Panel Card that holds the widgets.
        path_library (TextInput): Widget for the user to specify the library path.
        path_save_folder (TextInput): Widget for the user to specify the save folder.
        ptm (LiteralInput): Widget for the user to specify the PTMs.
        analysis_software (Select): Widget for the user to select the analysis software used to
        plot_mz: This EditableRangeSlider widget controls the m/z-range [Da].
        plot_im: This EditableRangeSlider widget controls the ion mobility range [1/K0].
        numbins: An IntInput widget to set the number of bins for the kernel density plot.
        load_library_descr (Markdown): Widget displaying the description for library loading.
        upload_button (Button): Button to upload selected library.
        upload_progress (Progress): Indicator to show progress during library upload.
        import_error (Alert): Alerts the user if there is a mismatch between input format and 
            analysis software.

    Methods:
        init: Constructor for the LoadLibraryCard class.
        create_layout: Creates the panel layout for the class.
        update_parameters: Updates the parameters in method configuration regarding input data.
        update_parameters_plotting: Updates the parameters in method configuration regarding plots.
        upload_data: Uploads library data, creates relevant folders, calculates kernel density, generates 
            library plots and visualisations.
    """
    def __init__(self, description):
        super().__init__(name="Data")
        self.library = None
        self.layout = None
        self.path_library = pn.widgets.TextInput(
            name='Specify the path to the library:',
            placeholder=library_path_placeholder,
            value=os.path.join(os.path.dirname(__file__), "static","evidence_MaxQuant_270223.txt"),
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
            sizing_mode='stretch_width',
            margin=(15, 15, 15, 15)
        )
        self.plot_mz = pn.widgets.EditableRangeSlider(
            name='m/z-range [Da]',
            start=100,
            end=1700,
            value=tuple(method_conf['graphs']['plot_mz']),
            step=50,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.plot_im = pn.widgets.EditableRangeSlider(
            name='Ion mobility range [1/K0]',
            start=0.6,
            end=1.6,
            value=tuple(method_conf['graphs']['plot_IM']),
            step=0.1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.numbins = pn.widgets.IntInput(
            name='Number of bins',
            start=1,
            end=300,
            value=method_conf['graphs']['numbins'],
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.load_library_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start',
        )
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
            "Input format and analysis software do not agree",
            alert_type="danger",
            sizing_mode='stretch_width',
            margin=(30, 15, 5, 15),
        )

    def create_layout(self):
        plot_widgets = pn.Column(
            pn.Row(
                pn.Column(self.plot_mz), 
                pn.Column(self.plot_im)),
            pn.Row(self.numbins, sizing_mode='stretch_width'),
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
            pn.layout.Divider(),
            pn.Row(
                None, None, None
            ),
            title='Load Library',
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
            self.path_library: [self.update_parameters, 'value'],
            self.path_save_folder: [self.update_parameters, 'value'],
            self.ptm: [self.update_parameters, 'value'],
            self.analysis_software: [self.update_parameters, 'value'],
            self.plot_mz: [self.update_parameters_plotting, 'value'],
            self.plot_im: [self.update_parameters_plotting, 'value'],
            self.numbins: [self.update_parameters_plotting, 'value'],
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
            os.path.join(self.path_save_folder.value, 'input_library'),
            os.path.join(self.path_save_folder.value, 'final_method'),
            os.path.join(self.path_save_folder.value, 'gif'),
        ]
        method_creator.create_folder(folder_paths)
        self.xi, self.yi, self.zi = plots.kernel_density_calculation(
            self.library,
            method_conf["graphs"]["numbins"]
        )
        self.layout[3][0] = pn.Column(
            pn.pane.Markdown(
                '### Histogram: Precursor distribution in m/z',
                 align='center'
            ),
            pn.pane.Matplotlib(
            plots.plot_precursor_distribution_as_histogram(
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
        )
        plots.plot_precursor_distribution_as_histogram(
            self.library,
            method_conf["graphs"],
            os.path.join(
                self.path_save_folder.value,
                'input_library',
                'Histogram_precursor_distribution_in_library.pdf'
            )
        )
        self.layout[3][1] = pn.Column(
            pn.pane.Markdown(
                '### Precursor cloud plotted across m/z and ion mobility',
                 align='center'
            ),
            pn.pane.Matplotlib(
            plots.density_plot(
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
        )
        plots.density_plot(
            self.xi,
            self.yi,
            self.zi,
            method_conf["graphs"],
            os.path.join(
                self.path_save_folder.value,
                'input_library',
                'Kernel_density_distribution_library.pdf'
            )
        )
        dict_charge_of_precursor = method_evaluator.calculate_percentage_multiple_charged_precursors(self.library)
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
            margin=(0, 50),
            sizing_mode='stretch_width',
            # align='center'
        )
        self.trigger_dependancy()
        self.upload_progress.active = False


class DefineScanAreaCard(BaseWidget):
    """
    `DefineScanAreaCard` class extends `BaseWidget` with properties and methods related to the scan 
    area, including its width, charge state and related input parameters for generating and evaluating 
    scan areas.

    Attributes:
        data: Imported data as a DataFrame for scan area calculations.
        description (str): Textual description of the class shown in the UI.
        layout (Card): Panel Card that holds the widgets.
        scan_area_descr (Markdown): Label for the description input in the interface.
        scan_area_width (IntInput): Widget for the user to specify the scan area width.
        charge_state (EditableRangeSlider): Widget for the user to specify the charge states for the scan.
        calculate_scan_area_button (Button): Button to initiate scan area calculations.
        scan_area_definition (LiteralInput): Widget to yield the calculated scan area definitions.
        evaluate_scan_area_button (Button): Button to initiate scan area evaluations.
        scan_area_spinner (LoadingSpinner): Activity indicator during scan area calculation and evaluations.

    Methods:
        init: Constructor for the LoadLibraryCard class.
        create_layout: Creates the panel layout for the class.
        update_parameters: Updates the parameters in method configuration regarding scan areas.
        run_calculations: Runs calculations for the defined scan area.
        run_evaluations: Runs evaluations on the calculated scan area.
    """
    def __init__(self, data, description):
        super().__init__(name="ScanArea")
        self.data = data
        self.layout = None
        self.scan_area_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start',
        )
        self.scan_area_width = pn.widgets.IntInput(
            name='Scan area width (horizontal) [m/z]',
            start=0,
            end=500,
            value=method_conf["scan_area"]["scan_area_width"],
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.charge_state = pn.widgets.EditableRangeSlider(
            name='Charge state',
            start=1,
            end=4,
            value=tuple(method_conf["scan_area"]["charge_state"]),
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.calculate_scan_area_button = pn.widgets.Button(
            name='Calculate Scan Area',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )
        self.scan_area_definition = pn.widgets.LiteralInput(
            name='Scan area',
            value=method_conf["scan_area"]["scan_area"],
            type=list,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.evaluate_scan_area_button = pn.widgets.Button(
            name='Plot and Evaluate Scan Area',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )
        self.scan_area_spinner = pn.indicators.LoadingSpinner(
            value=False,
            bgcolor='light',
            color='secondary',
            align='center',
            # margin=(0, 40, 0, 10),
            width=35,
            height=35
        )

    def create_layout(self):
        scan_area_widgets = pn.Column(
            pn.Row(self.scan_area_width, self.charge_state, sizing_mode='stretch_width'),
            pn.Spacer(height=20), 
            sizing_mode='stretch_width',
        )

        scan_area_def_widgets = pn.Column(
            pn.Row(self.scan_area_definition, sizing_mode='stretch_width'),
            pn.Spacer(height=20), 
            sizing_mode='stretch_width',
        )

        self.layout = pn.Card(
            self.scan_area_descr,
            pn.Row(
                pn.Column(
                    pn.WidgetBox(
                        scan_area_widgets,
                        sizing_mode='stretch_width',
                        margin=(20, 10, 50, 10),
                    ),
                    margin=(10, 30, 10, 10),
                    sizing_mode='stretch_width',
                ),
                pn.Column(
                    pn.Spacer(),
                    self.calculate_scan_area_button,
                    self.scan_area_spinner,
                    pn.Spacer(),
                    width=400,
                    align='center',
                ),
                sizing_mode='stretch_width',
            ),
            pn.layout.Divider(),
            pn.Row(None, None, None),
            pn.Row(
                pn.Column(
                    pn.WidgetBox(
                        scan_area_def_widgets,
                        sizing_mode='stretch_width',
                        margin=(20, 10, 50, 10),
                    ),
                    margin=(10, 30, 10, 10),
                    sizing_mode='stretch_width',
                ),
                pn.Column(
                    pn.Spacer(),
                    self.evaluate_scan_area_button,
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
            title='Define Scan Area',
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
            self.scan_area_width: [self.update_parameters, 'value'],
            self.charge_state: [self.update_parameters, 'value'],
            self.calculate_scan_area_button: [self.run_calculations, 'clicks'],
            self.scan_area_definition: [self.update_parameters, 'value'],
            self.evaluate_scan_area_button: [self.run_evaluations, 'clicks'],
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
            self.scan_area_width.name: "scan_area_width",
            self.charge_state.name: "charge_state",
            self.scan_area_definition.name: "scan_area",
        }
        method_conf['scan_area'][convertion_dict[event.obj.name]] = event.new

    def run_calculations(self, *args):
        self.scan_area_spinner.value = True
        self.scan_area_definition.value = method_creator.caclulate_scan_area_definitions(
            self.scan_area_width.value,
            self.charge_state.value,
            self.data.library
        )
        self.scan_area_spinner.value = False

    def run_evaluations(self, *args):        
        self.layout[6][0] = pn.Column(
            pn.pane.Markdown(
                '### Scan area of method',
                 align='center'
            ),
            pn.pane.Matplotlib(
            plots.density_plot_plus_scan_area(
                self.data.xi,
                self.data.yi,
                self.data.zi,
                method_conf["graphs"],
                os.path.join(
                    self.data.path_save_folder.value,
                    'input_library',
                    'Scan_Area.png'
                ),
                self.scan_area_definition.value,
                gui=True
            ),
            tight=True
            ),
            margin=(20, 0, 20, 0),
        )

        plots.density_plot_plus_scan_area(
            self.data.xi,
            self.data.yi,
            self.data.zi,
            method_conf["graphs"],
            os.path.join(
                self.data.path_save_folder.value,
                'input_library',
                'Scan_Area.pdf'
            ),
            self.scan_area_definition.value
        )

        coverage_center = method_evaluator.coverage_calculated_at_center(self.data.library, self.scan_area_definition.value)
        covered_precursor_info = pd.DataFrame(
            {
                "metrix for coverage": ["center of IM peak"],
                "coverage": [coverage_center]
            }
        )
        self.layout[6][1] = pn.Column(
            pn.pane.Markdown(
                '### Percentage of covered precursors',
                 align='center'
            ),
            pn.widgets.Tabulator(
            covered_precursor_info,
            layout='fit_data_table', 
            width=400
            ),
            margin=(20, 50),
        )
        self.trigger_dependancy()
        covered_precursor_info.to_csv(
            os.path.join(
                self.data.path_save_folder.value,
                'final_method',
                'covered_precursor_info.csv'
            ),
            index=False
        )


class GenerateMethodCard(BaseWidget):
    """
    Generates a MethodCard that extends the BaseWidget. It allows users to define the settings for 
    their synchro-PASEF method. This includes determining the number of scans, setting the window type, 
    determining the ratios of the scans widths, and the scan mode. Provides a user interface for 
    inputing these settings and updates the `method_conf` global variable accordingly.

    Attributes:
    layout : Placeholder for the layout of the widget card.
    data : Provided data for the card.
    scan_area : Provided scan area details.
    generate_method_descr : Markdown pane that displays description.
    scans (IntInput): Widget that allows user to input the number of synchro scans.
    window_type (Select): Widget that allows user to select the type of window (equidistant and variable).
    scan_ratio (LiteralInput): Widget that allows user to input the scan width ratio.
    scan_mode (Select): Widget that allows user to select the scan mode (classical synchro-PASEF, 
        highly accurate synchro-PASEF or individual synchro-PASEF). 
    window_modification (Select): Widget that allows user to modify the window (overlap, staggered).
    window_overlap (FloatInput): Widget that allows user to input the overlap of the window.
    no_of_combined_scans (IntInput): Widget that allows user to input the number of combined scans for
        a staggerd and classical synchro-PASEF window scheme.
    window_pattern (LiteralInput): Widget to get window pattern.
    MS1_positions (LiteralInput): Widget to get MS1 positions.
    im_limits (EditableRangeSlider): Widget that allows the user to choose ion mobility limits for the 
        method.
    generate_method_button (Button): A button widget to generate synchro-PASEF method.
    method_error (Alert): Alert that shows a message if the method is unsupported.
    path_method (TextInput): Widget that allows user to specify the path to the method file.
    method_creation_spinner: An indicator that shows a loading spinner during method creation.

    Methods:
    init: Constructor for the LoadLibraryCard class.
    create_layout: Creates the panel layout for this widget card, including setting up dependencies for 
        widget changes.
    update_parameters: Updates the global variable `method_conf` based on the current widget settings.
    update_general_parameters: Updates the general parameters in the `method_conf` global variable based 
        on widget settings.
    run_generation: If the conditions are met, it uses the `method_creator` object to generate a method 
        and update the layout with the resulting plots.
    """
    def __init__(self, data, scan_area, description):
        super().__init__(name="GenerateMethod")
        self.data = data
        self.scan_area = scan_area
        self.generate_method_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start',
        )
        self.scans = pn.widgets.IntInput(
            name='Scans',
            start=0,
            end=50,
            value=method_conf["method_parameters"]["scans"],
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.window_type = pn.widgets.Select(
            name='Window type',
            value=method_conf["method_parameters"]["window_type"],
            options=["equidistant", "variable", "pre-definied"],
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.scan_ratio = pn.widgets.LiteralInput(
            name='Scan ratio * must have as many items as number of scans',
            value=method_conf["method_parameters"]["scan_ratio"],
            type=list,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.scan_mode = pn.widgets.Select(
            name='Acquisition scheme',
            value=method_conf["method_parameters"]["scan_mode"],
            options=["classical_synchro-PASEF", "highly_accurate_synchro-PASEF", "individual_synchro-PASEF"],
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.window_modification = pn.widgets.Select(
            name='Window modification',
            value=method_conf["method_parameters"]["window_modification"],
            options=['None', 'overlap', 'staggered'],
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.window_overlap = pn.widgets.FloatInput(
            name='Window overlap *0.2 means 20% of smallest window',
            start=0,
            end=10,
            value=method_conf["method_parameters"]["window_overlap"],
            step=0.1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.no_of_combined_scans = pn.widgets.IntInput(
            name='No of combined scans *required for staggered + classical synchro-PASEF',
            start=0,
            end=1000,
            value=method_conf["method_parameters"]["no_of_combined_scans"],
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.window_pattern = pn.widgets.LiteralInput(
            name='Window pattern *required for individual_synchro-PASEF',
            value=method_conf["method_parameters"]["window_pattern"],
            type=list,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.MS1_positions = pn.widgets.LiteralInput(
            name='MS1 positions',
            value=method_conf["method_parameters_general"]["MS1_positions"],
            type=list,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.im_limits = pn.widgets.EditableRangeSlider(
            name='IM limits for the method',
            start=0.5,
            end=1.7,
            value=tuple([method_conf["method_parameters_general"]["dict_im_limits"]['low_limit_IM'],
                         method_conf["method_parameters_general"]["dict_im_limits"]['up_limit_IM']]),
            step=0.001,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.generate_method_button = pn.widgets.Button(
            name='Generate synchro-PASEF method',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )
        self.method_error = pn.pane.Alert(
            "This method is not supported by timsControl <= 6.0",
            alert_type="danger",
            sizing_mode='stretch_width',
            align='center',
            margin=(10, 0, 0, 0)
        )
        self.path_method = pn.widgets.TextInput(
            name='Specify the path to the method file:',
            placeholder=method_path_placeholder,
            value=method_conf['input']['method_only_used_for_evaluate'],
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )
        self.method_creation_progress = pn.indicators.Progress(
            active=False,
            bar_color='light',
            width=250,
            align='center',
            margin=(0, 0, 30, 0)
        )

    def create_layout(self):
        general_settings = pn.Column(
            pn.Row(self.im_limits, self.MS1_positions, sizing_mode='stretch_width'),
            pn.Row(self.scans, self.window_type, sizing_mode='stretch_width'),
            pn.Row(self.window_modification, self.window_overlap, sizing_mode='stretch_width'),
            pn.Spacer(height=20),
            sizing_mode='stretch_width',
        )

        additional_settings = pn.Card(
            pn.Column(
                pn.Row(self.scan_ratio, self.no_of_combined_scans, sizing_mode='stretch_width'),
                pn.Row(self.scan_mode, self.window_pattern, sizing_mode='stretch_width'),
                pn.Spacer(height=20),
                sizing_mode='stretch_width',
            ),
            title='Additional Method Settings',
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

        self.results_placeholder = pn.Row()
        self.error_placeholder = pn.Row()

        self.layout = pn.Card(
            self.generate_method_descr,
            pn.Row(
                pn.Column(
                    pn.WidgetBox(
                        general_settings,
                        sizing_mode='stretch_width',
                        margin=(20, 10, 50, 10),
                    ),
                    additional_settings,
                    sizing_mode='stretch_width',
                ),
                pn.Column(
                    pn.Spacer(),
                    self.generate_method_button,
                    self.method_creation_progress,
                    self.error_placeholder,
                    pn.Spacer(),
                    width=400,
                    align='center',
                ),
                sizing_mode='stretch_width',
            ),
            pn.layout.Divider(),
            self.results_placeholder, 

            title='Generate Synchro-PASEF Method',
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
            self.scans: [self.update_parameters, 'value'],
            self.window_type: [self.update_parameters, 'value'],
            self.scan_ratio: [self.update_parameters, 'value'],
            self.scan_mode: [self.update_parameters, 'value'],
            self.window_modification: [self.update_parameters, 'value'],
            self.window_overlap: [self.update_parameters, 'value'],
            self.no_of_combined_scans: [self.update_parameters, 'value'],
            self.window_pattern: [self.update_parameters, 'value'],
            self.MS1_positions: [self.update_general_parameters, 'value'],
            self.im_limits: [self.update_general_parameters, 'value'],
            self.generate_method_button: [self.run_generation, 'clicks'],
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
            self.scans.name: "scans",
            self.window_type.name: "window_type",
            self.scan_ratio.name: "scan_ratio",
            self.scan_mode.name: "scan_mode",
            self.window_modification.name: "window_modification",
            self.window_overlap.name: "window_overlap",
            self.no_of_combined_scans.name: "no_of_combined_scans",
            self.window_pattern.name: "window_pattern",
        }
        method_conf["method_parameters"][convertion_dict[event.obj.name]] = event.new

    def update_general_parameters(self, event):
        global method_conf
        convertion_dict = {
            self.MS1_positions.name: "MS1_positions",
            self.im_limits.name: "dict_im_limits",
        }
        method_conf["method_parameters_general"][convertion_dict[event.obj.name]] = event.new

    def run_generation(self, *args):
        self.method_creation_progress.active = True
        self.results_placeholder.clear()
        self.error_placeholder.clear()

        if (self.window_type.value != "equidistant") or (self.scan_mode.value != "classical_synchro-PASEF"):
            self.error_placeholder.append(self.method_error)

        self.scan_ratio.value, self.scans.value = method_creator.create_method(
            self.data.path_save_folder.value + '/final_method',
            self.scan_area.scan_area_definition.value,
            {
                "scans": self.scans.value,
                "window_type": self.window_type.value,
                "scan_ratio": self.scan_ratio.value,
                "scan_mode": self.scan_mode.value,
                "window_modification": self.window_modification.value,
                "window_overlap": self.window_overlap.value,
                "no_of_combined_scans": self.no_of_combined_scans.value,
                "window_pattern": self.window_pattern.value,
            },
            {
                "MS1_positions": self.MS1_positions.value,
                "dict_im_limits": {'low_limit_IM': self.im_limits.value[0], 'up_limit_IM': self.im_limits.value[1]},
            },
            self.data.library
        )

        final_method_path = os.path.join(
            self.data.path_save_folder.value,
            'final_method',
            'synchroPasef.txt'
        )
        self.path_method.value = final_method_path

        self.results_placeholder.append(
            pn.Column(
                pn.pane.Markdown('### Method Generated', align='center'),
                pn.widgets.TextInput(value=final_method_path, disabled=True),
                sizing_mode='stretch_width',
                align='center',
            )
        )

        self.trigger_dependancy()
        self.method_creation_progress.active = False


class LoadMethodCard(BaseWidget):
    """
    LoadMethodCard class extends BaseWidget with properties and methods that enable users to upload a 
    previously generated synchro-PASEF method. This class provides an interface that allows users to 
    interact with previously generated synchro-PASEF methods.

    Attributes:
    layout : Placeholder for the layout of the widget card.
    data : Provided data for the card.
    method_creation : Details of the method creation.
    load_method_descr (Markdown) : A pane that displays the description.
    load_button (Button) : A button widget that allows users to load synchro-PASEF method.

    Methods:
    init: Constructor for the LoadLibraryCard class.
    create_layout(): This function creates the panel layout for this widget card and sets up 
        dependencies for widget changes.
    update_input_parameters(event): This function updates the 'method_conf' global variable based on 
        the current widget settings.
    run_load(*args): Function that uploads the method data specified by the user and displays it in the 
        user interface.
    """
    def __init__(self, data, method_creation, description):
        super().__init__(name="Load Method")
        self.data = data
        self.method_creation = method_creation
        self.load_method_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start',
        )
        self.load_button = pn.widgets.Button(
            name='Load synchro-PASEF method',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0),
        )

    def create_layout(self):
        load_method_widgets = pn.Column(
            self.method_creation.path_method,
            pn.Spacer(height=20),
            sizing_mode='stretch_width',
        )

        self.layout = pn.Card(
            self.load_method_descr,
            pn.Row(
                pn.Column(
                    pn.WidgetBox(
                        load_method_widgets,
                        sizing_mode='stretch_width',
                        margin=(20, 10, 50, 10),
                    ),
                    margin=(10, 30, 10, 10),
                    sizing_mode='stretch_width',
                ),
                pn.Column(
                    pn.Spacer(),
                    self.load_button,
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
            title='Load synchro-PASEF Method',
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
            self.load_button: [self.run_load, 'clicks'],
            self.method_creation.path_method: [self.update_input_parameters, 'value'],
        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1],
                onlychanged=True
            )

        return self.layout
    
    def update_input_parameters(self, event):
        global method_conf
        convertion_dict = {
            self.method_creation.path_method.name: "method_only_used_for_evaluate",
        }
        method_conf["input"][convertion_dict[event.obj.name]] = event.new

    def run_load(self, *args):

        df_parameters_final = loader_MS_parameter_file.load_MS_method_from_txt_file(
            self.method_creation.path_method.value
        )
        self.df_method = df_parameters_final

        self.layout[3] = pn.Column(
            pn.pane.Markdown('### Loaded synchro-PASEF Method', align='center'),
            pn.widgets.Tabulator(
                df_parameters_final,
                margin=(20, 50),
                layout='fit_data_table', 
                sizing_mode='stretch_width',
            ),
            sizing_mode='stretch_width',
            align='center',
        )


class EvaluateMethodCard(BaseWidget):
    """
    EvaluateMethodCard class extends BaseWidget and is used to evaluate synchro-PASEF methods.

    This class generates graphical user interfaces which allow users to evaluate the synchro-PASEF 
    method using pre-acquired library. It provides functionalities such as manually specifying the
    ramp time and ion mobility range of the method in the evaluation process.

    Attributes:
    data: Instance of LoadLibraryCard to hold the loaded library data.
    generate_method: Instance of GenerateMethodCard to hold generated method.
    load_method: Instance of LoadMethodCard to load the method.
    layout (Card): Panel Card that holds the widgets.
    evaluate_method_descr (Markdown): Widget for the user to read the description.
    im (EditableRangeSlider): Widget for the user to control the ion mobility (IM) range/limits.
    library: Expected to be a DataFrame to hold the loaded library data.
    ramp_steps: An IntInput widget to set the TIMS ramp steps.
    evaluate_button: A button widget to trigger the evaluation process.
    upload_progress: A progress indicator for the evaluation process.

    Methods:
    init: Constructor for the LoadLibraryCard class.
    create_layout: Creates the panel layout for the class.
    update_acquisition_parameters: Update acquisition parameters.
    update_layout_according_to_library: Update the layout according to whether 'IMlength' is in the 
        data library.
    run_evaluation: Run the evaluation of the method and update the layout with the results.
    """
    def __init__(self, data, generate_method, load_method, description):
        super().__init__(name="Evaluate")
        self.data = data
        self.generate_method = generate_method
        self.load_method = load_method
        self.layout = None
        self.evaluate_method_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start',
        )
        self.im = self.generate_method.im_limits
        self.library = self.data.library
        self.ramp_steps = pn.widgets.IntInput(
            name='TIMS ramp time *927 = 100ms, 463 = 50ms, 1390 = 150ms, 1854 = 200 ms',
            start=0,
            end=3000,
            value=method_conf["graphs"]["ramp_steps"],
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.evaluate_button = pn.widgets.Button(
            name='Evaluate synchro-PASEF method',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0),
        )
        self.evaluate_progress = pn.indicators.Progress(
            active=False,
            bar_color='light',
            width=250,
            align='center',
            margin=(0, 0, 30, 0)
        )

    def create_layout(self):
        evaluate_method_widgets = pn.Column(
            pn.Row(self.im, self.ramp_steps, sizing_mode='stretch_width'),
            pn.Spacer(height=20),
            sizing_mode='stretch_width',
        )

        self.layout = pn.Card(
            self.evaluate_method_descr,
            pn.Row(
                pn.Column(
                    None,
                    margin=(10, 30, 10, 10),
                    sizing_mode='stretch_width',
                ),
                pn.Column(
                    pn.Spacer(),
                    self.evaluate_button,
                    self.evaluate_progress,
                    pn.Spacer(),
                    width=400,
                    align='center',
                ),
                sizing_mode='stretch_width',
            ),
            pn.layout.Divider(),
            pn.Row(
                None, None, None
            ),
            pn.Row(
                None,
            ),
            title='Evaluate Synchro-PASEF Method',
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
            self.im: [self.update_acquisition_parameters, 'value'],
            self.ramp_steps: [self.update_acquisition_parameters, 'value'],
            self.data.upload_button: [self.update_layout_according_to_library, 'clicks'],
            self.evaluate_button: [self.run_evaluation, 'clicks'],
        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1],
                onlychanged=True
            )

        return self.layout

    def update_acquisition_parameters(self, event):
        global method_conf
        convertion_dict = {
            self.im.name: "IM",
            self.ramp_steps.name: "ramp_steps",
        }
        method_conf["graphs"][convertion_dict[event.obj.name]] = event.new

    def update_layout_according_to_library(self, *args):
        if 'IMlength' in self.data.library.columns:
            self.layout[1][0][0] = pn.Row(self.im, self.ramp_steps, sizing_mode='stretch_width')
        else:
            self.layout[1][0][0] = None

    def run_evaluation(self, *args):
        self.evaluate_progress.active = True

        self.layout[3][0] = pn.Column(
            pn.pane.Markdown(
                '### Acquisition scheme highlighting m/z position and lengths of scans',
                 align='center'
            ),
            pn.pane.Matplotlib(
            plots.plot_acquisition_scheme(
                self.load_method.df_method,
                os.path.join(
                    self.data.path_save_folder.value,
                    'final_method',
                    "acquisition_scheme.png"
                ),
                gui=True
            ),
            tight=True
        )
        )

        self.layout[3][1] = pn.Column(
            pn.pane.Markdown(
                '### Acquisition scheme plotted on top of precursor cloud',
                 align='center'
            ),
            pn.pane.Matplotlib(
            plots.plot_method_and_precursors(
                self.data.xi,
                self.data.yi,
                self.data.zi,
                method_conf["graphs"],
                self.load_method.df_method,
                os.path.join(
                    self.data.path_save_folder.value,
                    'final_method',
                    "acquisition_scheme_and_density_plot.png"
                ),
                gui=True
            ),
            tight=True
        )
        )

        plots.plot_acquisition_scheme(
            self.load_method.df_method,
            os.path.join(
                self.data.path_save_folder.value,
                'final_method',
                "acquisition_scheme.pdf"
            ),
        )

        plots.plot_method_and_precursors(
            self.data.xi,
            self.data.yi,
            self.data.zi,
            method_conf["graphs"],
            self.load_method.df_method,
            os.path.join(
                self.data.path_save_folder.value,
                'final_method',
                "acquisition_scheme_and_density_plot.pdf"
            ),
        )

        if ('IMlength' in self.data.library.columns) & (self.ramp_steps.value != 0):
            self.layout[3][2] = pn.Column(
            pn.pane.Markdown(
                '### Histogram showing ion mobility length of precursors and slicing cutoff per scan',
                align='center'
            ),
            pn.pane.Matplotlib(
                plots.histogram_precursor_slicing(
                    300,
                    self.data.library,
                    method_conf["graphs"],
                    self.load_method.df_method,
                    os.path.join(
                        self.data.path_save_folder.value,
                        'final_method',
                        "histogram_slicing.png"
                    ),
                    gui=True
                ),
                tight=True
            )
            )
            plots.histogram_precursor_slicing(
                300,
                self.data.library,
                method_conf["graphs"],
                self.load_method.df_method,
                os.path.join(
                    self.data.path_save_folder.value,
                    'final_method',
                    "histogram_slicing.pdf"
                )
            )
        else:
            self.layout[3][2] = None

        if "IMlength" not in self.data.library.columns:
            df_evaluated, _ = method_evaluator.calculate_coverage_total_per_scan_per_charge_state(self.load_method.df_method, self.data.library)
        else:
            df_coverage, df_temp = method_evaluator.calculate_coverage_total_per_scan_per_charge_state(self.load_method.df_method, self.data.library)
            df_slicing = method_evaluator.calculate_slicing_and_coverage_in_total(self.data.library, df_temp)
            df_evaluated = pd.concat([df_coverage, df_slicing])
        self.layout[4] = pn.widgets.Tabulator(
            df_evaluated,
            margin=(0, 0, 20, 100),
            # align='center',
            layout='fit_data_table', 
            width = 1200,
        )
        df_evaluated.to_csv(
            os.path.join(
                self.data.path_save_folder.value,
                'final_method',
                'df_evaluated.csv'
            )
        )

        self.evaluate_progress.active = False


class MakeGifCard(object):
    """
    MakeGifCard class includes properties and methods that provide users the ability to create a GIF of 
    a Synchro-PASEF method.

    Attributes:
    layout: Placeholder for the layout of the widget card.
    data: Provided data for the card.
    load_method: Details of the method loaded.
    make_gif_descr (Markdown): A pane that displays the description.
    gif_button (Button): A button widget that allows users to create a GIF of Synchro-PASEF method.
    upload_progress (Progress: An indicator that displays a progress bar during GIF upload.
    gif_duration (FloatInput): Widget that allows user to input the time per frame in ms.
    im_steps (IntInput): Widget that allows user to input the ion mobility steps.
    scans_plotted_at_once (IntInput): Widget that allows user to input the TOF triggers highlighted at once.

    Methods:
    init: Constructor for the LoadLibraryCard class.
    create_layout(): This function creates the panel layout for this widget card and sets up 
        dependencies for widget changes.
    run_gif(*args): Function that creates a GIF of a Synchro-PASEF method and displays it in the user 
        interface.
    """
    def __init__(self, data, load_method, description):
        # super().__init__(name="Make GIF")
        self.data = data
        self.load_method = load_method
        self.make_gif_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start',
        )
        self.gif_button = pn.widgets.Button(
            name='Create Synchro-PASEF method GIF',
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
            value=0.002,
            step=0.001,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.im_steps = pn.widgets.IntInput(
            name='Ion mobility steps',
            start=0,
            end=200,
            value=50,
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )
        self.scans_plotted_at_once = pn.widgets.IntInput(
            name='TOF triggers highlighted at once',
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
            title='Make a GIF of a Synchro-PASEF Method',
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

        files = glob.glob(save_at + "/*")
        for f in files:
            os.remove(f)

        boxes, df_temp_reset_index = plots.generate_boxes(self.load_method.df_method, self.im_steps.value)
        plots.generate_gif_single_windows(
            self.data.xi,
            self.data.yi,
            self.data.zi,
            method_conf["graphs"],
            boxes,
            range(1, len(df_temp_reset_index)+1),
            save_at + "/",
            facecolor="#FF0098",
            window_color="white_grey",
            scans_plotted_at_once=self.scans_plotted_at_once.value,
        )

        images = []
        for file_name in sorted(os.listdir(save_at)):
            if file_name.endswith('.png'):
                file_path = os.path.join(save_at, file_name)
                images.append(iio.imread(file_path))
        iio.imwrite(os.path.join(save_at, "SynchroPasef_single_windows.gif"), images, duration=self.gif_duration.value, loop=0)

        self.layout[3] = pn.Column(
            pn.pane.Markdown('### Synchro-PASEF Method GIF', align='center'),
            pn.pane.GIF(
                os.path.join(save_at, "SynchroPasef_single_windows.gif"),
                alt_text="Sorry, please try using a different browser",
                width=640, height=360
            ),
            sizing_mode='stretch_width',
            align='center',
        )

        self.gif_progress.active = False


class SynchroCardWidget(object):
    """
    SynchroCardWidget class contains descriptions and methods for loading the library, defining the scan 
    area, generating methods, loading methods, evaluating methods, and visualizing the synchro scan with 
    the help of the Card classes defined for each task.

    Attributes:
    layout: None by default, will be updated later.
    load_library_description (str): Descriptive text shown to the user during the library loading stage.
    scan_area_description (str): Descriptive text shown to the user while defining the scan area.
    generate_method_description (str): Descriptive text shown to the user while generating the methods.
    load_method_description (str): Descriptive text shown to the user when loading the methods.
    evaluate_method_description (str): Descriptive text shown to the user during the method evaluation 
        stage.
    make_gif_description (str): Descriptive text shown to the user while making the GIF.
    data: Instance of LoadLibraryCard to carry out the task of loading the library.
    scan_area: Instance of DefineScanAreaCard to carry out the task of defining the scan area.
    generate_method: Instance of GenerateMethodCard to carry out the task of generating the method.
    load_method: Instance of LoadMethodCard to carry out the task of loading the method.
    evaluate_method: Instance of EvaluateMethodCard to carry out the task of evaluating the method.
    make_gif_method: Instance of MakeGifCard to carry out the task of GIF creation.

    Methods:
    init: Constructor for the LoadLibraryCard class.
    create_layout: Creates and returns the panel layout for the class with the layout of each of the 
        cards as column elements.

    Note:
    The create_layout method takes an optional parameter 'tab_list'.
    """
    def __init__(self):
        self.layout = None
        self.load_library_description = "#### Load the library for the indicated analysis software to check the distribution of the precursors in the m/z-ion mobility plane."
        self.scan_area_description = "#### Define the scan area dependent on the precursor density."
        self.generate_method_description = "#### Generated a synchro-PASEF method with individual method design."
        self.load_method_description = "#### Load a synchro-PASEF method from a txt file or directly from the raw file."
        self.evaluate_method_description = "#### Evaluate a synchro-PASEF method based on a pre-acquired library. Optionally specify the ion mobility limits and TIMS ramp time of the method above."
        self.make_gif_description = "#### Visualize the synchro scan movement of your method with a GIF."

        self.data = LoadLibraryCard(self.load_library_description)
        self.scan_area = DefineScanAreaCard(
            self.data,
            self.scan_area_description
            )
        self.generate_method = GenerateMethodCard(
            self.data,
            self.scan_area,
            self.generate_method_description
            )
        self.load_method = LoadMethodCard(
            self.data,
            self.generate_method,
            self.load_method_description
            )
        self.evaluate_method = EvaluateMethodCard(
            self.data,
            self.generate_method,
            self.load_method,
            self.evaluate_method_description,
            )
        self.make_gif_method = MakeGifCard(
            self.data,
            self.load_method,
            self.make_gif_description
            )

    def create_layout(
        self,
        tab_list=None
    ):
        return pn.Column(
            self.data.create_layout(),
            self.scan_area.create_layout(),
            self.generate_method.create_layout(),
            self.load_method.create_layout(),
            self.evaluate_method.create_layout(),
            self.make_gif_method.create_layout(),
            sizing_mode='stretch_width',
        )
