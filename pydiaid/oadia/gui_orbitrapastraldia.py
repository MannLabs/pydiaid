import pandas as pd
import os
import platform
import json

# visualization
import panel as pn

# local imports
import pydiaid.oadia.method_generator as method_creator
import pydiaid.oadia.method_evaluation as method_evaluator
import pydiaid.loader as loader
import pydiaid.synchropasef.method_creator as syP_method_creator
import pydiaid.synchropasef.method_evaluator as syP_method_evaluator
import pydiaid.synchropasef.plots as syP_plots

# paths
BASE_PATH = os.path.dirname(__file__)
STYLE_PATH = os.path.join(BASE_PATH, "style")
DOCS_PATH = os.path.join(BASE_PATH, "docs")
IMG_PATH = os.path.join(BASE_PATH, "img")
DEFAULT_FILE = os.path.join(
    BASE_PATH,
    'static',
    "default_parameters.json"
)

with open(DEFAULT_FILE, "r") as infile:
    method_conf = json.load(infile)

# Platform-specific defaults
if platform.system() == 'Windows':
    method_path_placeholder = 'D:\pydiaid\pydiaid\oadia\static\MassListTable_targeted.csv'
    library_path_placeholder = 'D:\pydiaid\pydiaid\oadia\static\AlphaPept_results.csv'
    save_folder_placeholder = 'D:\pydiaid\pydiaid\static'
else:
    method_path_placeholder = '/Users/patriciaskowronek/Documents/pydiaid/pydiaid/oadia/static/MassListTable_targeted.csv'
    library_path_placeholder = '/Users/patriciaskowronek/Documents/pydiaid/pydiaid/oadia/static/AlphaPept_results.csv'
    save_folder_placeholder = '/Users/pydiaid/static'

class BaseWidget(object):
    """
    BaseWidget class that initializes a base widget with given name and contains methods to track
    updates and trigger dependencies.
    """
    def __init__(self, name):
        self.name = name
        self.__update_event = pn.widgets.IntInput(value=0)
        self.depends = pn.depends(self.__update_event.param.value)
        self.active_depends = pn.depends(self.__update_event.param.value, watch=True)
        
    def trigger_dependancy(self):
        self.__update_event.value += 1

class LoadLibraryCard(BaseWidget):
    """
    LoadLibraryCard class extends BaseWidget with properties and methods for loading and analyzing 
    proteomics libraries for Orbitrap Astral DIA analysis.

    This class provides graphical user interface components to input parameters for loading a proteomics 
    library, including file path, analysis software and post-translational modifications (PTM), and 
    provides components for visualizing precursor distributions as a histogram.

    Attributes:
        description (str): Textual description of card shown in the UI.
        library: Expected to be a DataFrame to hold the loaded library data.
        layout (Card): Panel Card that holds the widgets.
        path_library (TextInput): Widget for the user to specify the library path.
        path_save_folder (TextInput): Widget for the user to specify the save folder.
        ptm (LiteralInput): Widget for the user to specify the PTMs.
        analysis_software (Select): Widget for the user to select the analysis software.
        plot_mz (EditableRangeSlider): Widget controlling the m/z-range visualization in plot [Da].
        upload_button (Button): Button to upload and analyze selected library.
        upload_progress (Progress): Indicator to show progress during library upload.

    Methods:
        init: Constructor for the LoadLibraryCard class.
        create_layout: Creates the panel layout for the class.
        update_parameters: Updates the parameters in method configuration for input data.
        update_parameters_plotting: Updates the parameters for plotting configuration.
        upload_data: Handles library upload, creates folders, and generates visualizations.
    """

    def __init__(self, description):
        super().__init__(name="Data")
        self.library = None
        self.layout = None
        
        # Initialize widgets
        self.path_library = pn.widgets.TextInput(
            name='Specify the path to the library:',
            placeholder=library_path_placeholder,
            value=os.path.join(os.path.dirname(__file__), "static", "AlphaPept_results.csv"),
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )
        
        self.path_save_folder = pn.widgets.TextInput(
            name='Save the output to the following folder:',
            placeholder=save_folder_placeholder,
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )
        
        self.analysis_software = pn.widgets.Select(
            name='Analysis software',
            options=['AlphaPept', 'MaxQuant', 'MSFragger', 'Spectronaut single-run', 'Spectronaut library', "DIANN single-run", "DIANN library", 'AlphaPeptDeep library', 'AlphaDIA'],
            min_width=200,
            margin=(15, 15, 15, 15)
        )
        
        self.ptm = pn.widgets.LiteralInput(
            name='Specify the PTM:',
            # value = 'None',
            placeholder="['Phospho']",
            sizing_mode='stretch_width',
            margin=(15, 15, 15, 15)
        )
        
        self.plot_mz = pn.widgets.EditableRangeSlider(
            name='Plot m/z-range [Da]',
            start=100,
            end=1700,
            value=tuple(method_conf['graphs']['plot_mz']),
            step=50,
            margin=(15, 15, 0, 15)
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

        self.load_library_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start'
        )

    def create_layout(self):
        """Creates and returns the widget layout."""

        self.layout = pn.Card(
            self.load_library_descr,
            pn.Row(
                pn.Column(
                    self.path_library,
                    self.path_save_folder,
                    pn.Row(self.analysis_software, self.ptm),
                    self.plot_mz,
                    margin=(10, 30, 10, 10),
                    sizing_mode='stretch_width'
                ),
                pn.Column(
                    pn.Spacer(),
                    self.upload_button,
                    self.upload_progress,
                    pn.Spacer(),
                    width=400,
                    align='center'
                ),
                sizing_mode='stretch_width'
            ),
            pn.layout.Divider(),
            pn.Row(None, None),
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

        # Set up event handlers
        self.upload_button.on_click(self.upload_data)
        return self.layout
    
    def update_parameters_plotting(self, event):
        global method_conf
        convertion_dict = {
            self.plot_mz.name: "plot_mz",
        }
        method_conf['graphs'][convertion_dict[event.obj.name]] = event.new

    def upload_data(self, event):
        """Handles library upload and initial analysis."""
        self.upload_progress.active = True
        
        # Load library
        self.library = loader.load_library(
            self.path_library.value,
            method_conf["input"]["analysis_software"],
            method_conf["input"]["PTM"]
        )

        # Create output directories
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
        syP_method_creator.create_folder(folder_paths)
        self.layout[3][0] = pn.Column(
            pn.pane.Markdown(
                '### Histogram: Precursor distribution in m/z',
                 align='center'
            ),
            pn.pane.Matplotlib(
            syP_plots.plot_precursor_distribution_as_histogram(
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
    
        dict_charge_of_precursor = syP_method_evaluator.calculate_percentage_multiple_charged_precursors(self.library)
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
        self.layout[3][1] = pn.Column(
            pn.pane.Markdown(
                '### Percentage of multiple charged precursors',
                 align='center'
            ),
            pn.widgets.Tabulator(
                mult_charged_precursor_info,
                layout='fit_data_table', 
                width=400
            ),
            # margin=(0, 50),
            sizing_mode='stretch_width',
            # align='center'
        )
        self.trigger_dependancy()
        self.upload_progress.active = False


class SpecifyMethodParametersCard(BaseWidget):
    """
    SpecifyMethodParametersCard class for configuring DIA method parameters.
    Extends BaseWidget with properties and methods for specifying DIA windows.

    Attributes:
        data: Reference to loaded library data.
        layout (Card): Panel Card that holds the widgets.
        mz_range (EditableRangeSlider): Controls the m/z range for window generation [Da].
        num_bins (IntInput): Sets number of DIA scans per method.
        window_type (Select): Selects between fixed or dynamic window widths.
        output_format (Select): Determines output format (all, targeted, center_mass, mz_ranges).
        min_width (FloatInput): Sets minimum allowed window width [Da].
        max_width (IntInput): Sets maximum allowed window width [Da].
        calculate_button (Button): Triggers parameter calculation and validation.
        specify_parameter_descr (Markdown): Displays description and instructions.

    Methods:
        create_layout: Creates and organizes the parameter input interface.
        update_parameters: Updates method configuration based on widget changes.
        calculate_method: Validates parameters and calculates precursor coverage.
    """
    def __init__(self, data, description):
        super().__init__(name="Method")
        self.data = data
        self.layout = None

        # Initialize widgets
        self.mz_range = pn.widgets.EditableRangeSlider(
            name='m/z range [Da]',
            start=100,
            end=1700,
            value=(380, 980),
            step=50,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width'
        )

        self.num_bins = pn.widgets.IntInput(
            name='Number of scans',
            start=1,
            end=3000,
            value=60,
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width'
        )

        self.window_type = pn.widgets.Select(
            name='Window type',
            options=['fixed', 'dynamic'],
            value='dynamic',
            margin=(15, 15, 0, 15),
            min_width=200
        )

        self.output_format = pn.widgets.Select(
            name='Output format',
            options=['all', 'targeted', 'center_mass', 'mz_ranges'],
            value='all',
            margin=(15, 15, 0, 15),
            min_width=200
        )

        self.rt_time = pn.widgets.FloatInput(
            name='RT (min) (only for targeted output format)',
            start=1,
            end=300,
            value=16.5,
            step=0.1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width'
        )

        self.window = pn.widgets.IntInput(
            name='Window (only for targeted output format)',
            start=1,
            end=300,
            value=33,
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width'
        )

        self.min_width = pn.widgets.FloatInput(
            name='Lower limit (m/z width)',
            start=0.2,
            end=50,
            value=2,
            step=0.01,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width'
        )

        self.max_width = pn.widgets.IntInput(
            name='Upper limit (m/z width)',
            start=1,
            end=300,
            value=50,
            step=1,
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width'
        )

        self.forbidden_zone = pn.widgets.Checkbox(
            name='Adjust window boarders to forbidden zones',
            value=False,  # default value
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )

        self.phospho_method = pn.widgets.Checkbox(
            name='Calculate forbidden zones for phosphopeptides',
            value=False,  # default value
            margin=(15, 15, 0, 15),
            sizing_mode='stretch_width',
        )

        self.calculate_button = pn.widgets.Button(
            name='Calculate Method Parameters',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )

        self.specify_parameter_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start'
        )

    def create_layout(self):
        """Creates and returns the widget layout."""
        parameter_widgets = pn.Column(
            pn.Row(pn.Column(self.mz_range), pn.Column(self.num_bins)),
            pn.Row(self.window_type, self.output_format, sizing_mode='stretch_width'),
            pn.Row(pn.Column(self.min_width), pn.Column(self.max_width)),
            pn.Row(pn.Column(self.rt_time), pn.Column(self.window)),
            pn.Row(pn.Column(self.forbidden_zone), pn.Column(self.phospho_method)),
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
                    align='center'
                ),
                sizing_mode='stretch_width'
            ),
            pn.layout.Divider(),
            pn.Row(None, None),
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
            self.mz_range: [self.update_parameters, 'value'],
            self.num_bins: [self.update_parameters, 'value'],
            self.window_type: [self.update_parameters, 'value'],
            self.output_format: [self.update_parameters, 'value'],
            self.min_width: [self.update_parameters, 'value'],
            self.max_width: [self.update_parameters, 'value'],
            self.forbidden_zone: [self.update_parameters, 'value'],
            self.phospho_method: [self.update_parameters, 'value'],
            self.rt_time: [self.update_parameters, 'value'],
            self.window: [self.update_parameters, 'value'],
            self.calculate_button: [self.calculate_method, 'clicks'],
        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1],
                onlychanged=True
            )

        return self.layout

    def update_parameters(self, event):
        """Updates method configuration parameters based on widget changes."""
        global method_conf
        convertion_dict = {
            self.mz_range.name: "mz_range",
            self.num_bins.name: "num_bins",
            self.window_type.name: "window_type",
            self.output_format.name: "output_format",
            self.min_width.name: "min_width",
            self.max_width.name: "max_width",
            self.forbidden_zone.name: "forbidden_zone",
            self.phospho_method.name: "forbidden_zones_for_phospho",
            self.rt_time: "rt_time",
            self.window: "window"
        }
        method_conf['method_parameters'][convertion_dict[event.obj.name]] = event.new

    def calculate_method(self, event):
        """Handles method parameter calculation and result display."""

        # Verify library is loaded
        if self.data.library is None:
            raise ValueError("Please load a library first")

        # Calculate and display statistics
        dict_precursors_within_scan_area = method_evaluator.calculate_precursor_within_scan_area(
            self.data.library,
            method_conf['method_parameters']["mz_range"],
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


class CreateMethodCard(BaseWidget):
    """
    CreateMethodCard class for generating optimized DIA methods.
    Extends BaseWidget with properties and methods for creating DIA windows based on 
    specified parameters.

    Attributes:
        data: Reference to loaded library data.
        param: Reference to method parameters.
        create_button (Button): Triggers method generation.
        path_method (TextInput): Shows path where method will be saved.
        create_method_descr (Markdown): Displays description and instructions.

    Methods:
        create_layout: Creates the method generation interface.
        create_method: Generates DIA windows based on parameters and library distribution.
    """
    def __init__(self, data, param, description):
        super().__init__(name="Method_creation")
        self.data = data
        self.param = param
        
        self.create_button = pn.widgets.Button(
            name='Create Method',
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
            align='start'
        )

    def create_layout(self):
        self.results_placeholder = pn.Row()

        self.layout = pn.Card(
            self.create_method_descr,
            pn.Row(
                pn.Column(
                    None
                ),
                pn.Column(
                    pn.Spacer(),
                    self.create_button,
                    pn.Spacer(),
                    width=400,
                    align='center'
                ),
                sizing_mode='stretch_width'
            ),
            pn.layout.Divider(),
            self.results_placeholder, 
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
            self.create_button: [self.create_method, 'clicks'],
        }
        for k in dependances.keys():
            k.param.watch(
                dependances[k][0],
                dependances[k][1],
                onlychanged=True
            )
        return self.layout

    def create_method(self, event):
        try:
            self.results_placeholder.clear()
            # Verify library is loaded
            if self.data.library is None:
                raise ValueError("Please load a library first")
            
            # Create method
            if self.param.output_format.value == "all":
                final_method_path = os.path.join(
                    method_conf["input"]["save_at"],
                    'final_method',
                    f"OA_DIA_method_{self.param.mz_range.value[0]}_{self.param.mz_range.value[1]}_{self.param.num_bins.value}_{self.param.window_type.value}_center_mass.csv"
                )
            else:
                final_method_path = os.path.join(
                    method_conf["input"]["save_at"],
                    'final_method',
                    f"OA_DIA_method_{self.param.mz_range.value[0]}_{self.param.mz_range.value[1]}_{self.param.num_bins.value}_{self.param.window_type.value}_{self.param.output_format.value}.csv"
                )


            final_method_folder = os.path.join(
                method_conf["input"]["save_at"],
                'final_method'
            )

            df_window, bins = method_creator.create_method(
                mz_values=self.data.library["mz"],
                mz_range=self.param.mz_range.value,
                num_bins=self.param.num_bins.value,
                window_type=self.param.window_type.value,
                output_format=self.param.output_format.value,
                folder_path=final_method_folder,
                min_width=self.param.min_width.value,
                max_width=self.param.max_width.value,
                adjusted_for_forbidden_zones = self.param.forbidden_zone.value,
                phospho_method = self.param.phospho_method.value
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

        except Exception as e:
            self.layout[3] = pn.pane.Alert(
                f"Error creating method: {str(e)}",
                alert_type='danger'
            )

        finally:
            self.trigger_dependancy()


class EvaluateMethodCard(BaseWidget):
    """
    EvaluateMethodCard class for analyzing and visualizing DIA methods.
    Extends BaseWidget with properties and methods for comprehensive method evaluation.

    Attributes:
        data: Reference to loaded library data.
        method_creation: Reference to created method details.
        evaluate_button (Button): Triggers method evaluation.
        evaluate_method_descr (Markdown): Displays evaluation instructions.
        evaluate_progress (Progress): Shows evaluation progress.

    Methods:
        create_layout: Creates the evaluation interface.
        evaluate_method: Analyzes isolation window distribution and generates visualizations.
        update_parameters: Updates evaluation parameters.
    """
    def __init__(self, data, method_creation, description):
        super().__init__(name="Method_evaluation")
        self.data = data
        self.method_creation = method_creation
        
        # Initialize widgets
        self.evaluate_button = pn.widgets.Button(
            name='Load and Evaluate Method',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )

        self.evaluate_method_descr = pn.pane.Markdown(
            description,
            margin=(5, 0, 2, 15),
            align='start'
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
                    align='center'
                ),
                sizing_mode='stretch_width'
            ),
            pn.layout.Divider(),
            pn.Row(None),
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
            self.evaluate_button: [self.evaluate_method, 'clicks'],
            self.method_creation.path_method: [self.update_parameters, 'value'],
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
            self.method_creation.path_method.name: "dia_method_only_used_for_evaluate",
        }
        method_conf['input'][convertion_dict[event.obj.name]] = event.new
    

    def evaluate_method(self, event):

        try:
            # Verify library is loaded
            if self.data.library is None:
                raise ValueError("Please load a library first")

            # Load and evaluate method
            _, bins, df_method = method_evaluator.load_method_file(
                method_conf["input"]["dia_method_only_used_for_evaluate"]
            )

            # Get statistics
            stats = method_evaluator.analyze_bins(bins, self.data.library["mz"])
        
            df_window = method_evaluator.parse_stats_text(stats['tabular_stats'])

            # Determine window type
            window_type = "fixed" if stats['max_width']-stats['min_width']<2 else "dynamic"

            # Update layout with results in three columns
            self.layout[3] = pn.Row(
                # Left column: Method details
                pn.Column(
                    pn.pane.Markdown(
                        '### Method Table',
                        align='center'
                    ),
                    pn.widgets.Tabulator(
                        df_method,
                        layout='fit_data_table',
                        height=400
                    )
                ),  

                # Middle column: Method summary
                pn.Column(
                    pn.pane.Markdown(
                        '### Method Statistics',
                        align='center'
                    ),
                    pn.widgets.Tabulator(
                        pd.DataFrame({
                            'Metric': [
                                'Number of scans',
                                'Min width (Da)',
                                'Max width (Da)',
                                'Average width (Da)',
                                'Min items per scan',
                                'Max items per scan',
                                'Average items per scan',
                                'Window type'
                            ],
                            'Value': [
                                stats['num_bins'],
                                f"{stats['min_width']:.2f}",
                                f"{stats['max_width']:.2f}",
                                f"{stats['avg_width']:.2f}",
                                stats['min_count'],
                                stats['max_count'],
                                f"{stats['avg_count']:.2f}",
                                window_type
                            ]
                        }),
                        layout='fit_data_table',
                        width=400
                    ),
                ),
                # Right column: Distribution plot
                pn.Column(
                    pn.pane.Markdown(
                        '### Precursor per Window',
                        align='center'
                    ),
                    pn.pane.Matplotlib(
                        method_evaluator.plot_precursors_per_scan(
                            window_type,
                            df_window,
                            os.path.join(
                                self.data.path_save_folder.value,
                                'final_method',
                                'precursor_per_scan_evaluation.png'
                            ),
                            gui=True
                        ),
                        tight=True,
                    )
                )
            )

            # Save evaluation results
            stats_df = pd.DataFrame({
                'Metric': [
                    'Number of bins',
                    'Min width',
                    'Max width',
                    'Average width',
                    'Min items per bin',
                    'Max items per bin',
                    'Average items per bin',
                    'Window type'
                ],
                'Value': [
                    stats['num_bins'],
                    stats['min_width'],
                    stats['max_width'],
                    stats['avg_width'],
                    stats['min_count'],
                    stats['max_count'],
                    stats['avg_count'],
                    window_type
                ]
            })

            stats_df.to_csv(
                os.path.join(
                    self.data.path_save_folder.value,
                    'final_method',
                    'method_evaluation.csv'
                ),
                index=False
            )
            method_evaluator.plot_precursors_per_scan(
                window_type,
                df_window,
                os.path.join(
                    self.data.path_save_folder.value,
                    'final_method',
                    'precursor_per_scan_evaluation.pdf'
                ),
            ),

        except Exception as e:
            self.layout[3] = pn.pane.Alert(
                f"Error evaluating method: {str(e)}",
                alert_type='danger'
            )

        finally:
            self.trigger_dependancy()


class OrbitrapAstralCardWidget:
    """
    Main GUI class for Orbitrap Astral DIA Windows method development.
    Coordinates the complete workflow from library loading to method evaluation.

    This class integrates multiple specialized cards to provide a comprehensive 
    interface for:
    - Loading and analyzing proteomics libraries
    - Specifying DIA window parameters
    - Generating optimized methods
    - Evaluating method performance

    Attributes:
        load_library_description (str): Description for library loading interface.
        specify_parameter_descr (str): Description for parameter specification.
        create_method_descr (str): Description for method creation.
        evaluate_descr (str): Description for method evaluation.
        data (LoadLibraryCard): Handles library loading and analysis.
        param (SpecifyMethodParametersCard): Handles parameter specification.
        create (CreateMethodCard): Handles method generation.
        evaluate (EvaluateMethodCard): Handles method evaluation.

    Methods:
        create_layout: Assembles the complete interface layout.
    """
    def __init__(self):
        self.load_library_description = """
        #### Load a proteomics library to analyze precursor distribution in the m/z dimension for a selected analysis software.
        """
        
        self.specify_parameter_descr = """
        #### Specify DIA window parameters. Recommended m/z-range is 380-980. Choose the number of DIA scans based on the sample concentration and chromatographic peak width. Select window type (fixed/dynamic) and output format.
        """
        
        self.create_method_descr = """
        #### Generate DIA windows optimized for the precursor distribution. Window widths smaller than the lower limit will be automatically merged, which may result in fewer scans than initially specified.
        """

        self.evaluate_descr = """
        #### Load and evaluate a DIA method using the pre-acquired library.
        """
        
        self.data = LoadLibraryCard(self.load_library_description)
        self.param = SpecifyMethodParametersCard(
            self.data,
            self.specify_parameter_descr
        )
        self.create = CreateMethodCard(
            self.data,
            self.param,
            self.create_method_descr
        )

        self.evaluate = EvaluateMethodCard(
            self.data,
            self.create,
            self.evaluate_descr
        )

    def create_layout(self):
        """Creates the main application layout."""
        return pn.Column(
            self.data.create_layout(),
            self.param.create_layout(),
            self.create.create_layout(),
            self.evaluate.create_layout(),
        sizing_mode='stretch_width'
        )