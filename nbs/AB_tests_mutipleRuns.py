# Author: Patricia Skowronek, Max Planck Institute of Biochemistry

# Standard library imports
import math
import re
from difflib import SequenceMatcher
from typing import List, Dict, Any, Tuple, Callable, Union

# Data manipulation and analysis
import numpy as np
import pandas as pd

## Visualization libraries
### Matplotlib and related
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib_venn import venn2, venn3

### Seaborn
import seaborn as sns

### Plotly
import plotly.express as px
from plotly.express.colors import sample_colorscale
import plotly.graph_objects as go
from plotly.subplots import make_subplots

### Other visualization tools
import datashader as ds
from upsetplot import plot as upset_plot

# Matplotlib configuration
plt.rcParams.update({'font.size': 16.5})
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

# Citation for UpSet plot
# Alexander Lex, Nils Gehlenborg, Hendrik Strobelt, Romain Vuillemot, Hanspeter Pfister,
# UpSet: Visualization of Intersecting Sets, IEEE Transactions on Visualization and Computer Graphics (InfoVis '14),
# vol. 20, no. 12, pp. 1983–1992, 2014. doi: doi.org/10.1109/TVCG.2014.2346248

def process_proteomics_data(
    file_paths: List[str],
    min_data_completeness: float,
    output_directory: str,
    dataset_labels: List[str]
) -> Tuple[List[str], Dict[str, List[Any]], str, List[str], List[pd.DataFrame], List[pd.DataFrame]]:
    analysis_results: Dict[str, List[Any]] = {
        "totalIDs": [],
        "CV20Percent": [],
        "CV10Percent": [],
        "CV5Percent": [],
        "IDs_per_run": [],
        "CV": [],
        "DataCompleteness": []
    }
    intensity_only_dataframes: List[pd.DataFrame] = []
    full_dataframes: List[pd.DataFrame] = []

    for index, file_path in enumerate(file_paths):
        df: pd.DataFrame = pd.read_csv(file_path, sep='\t')
        data_source: str = identify_data_source(df)
        
        if data_source:
            column_types, analysis_results, intensity_only_dataframes, full_dataframes = analyze_dataframe(
                df=df,
                analysis_results=analysis_results,
                min_data_completeness=min_data_completeness,
                output_directory=output_directory,
                data_source=data_source,
                dataset_label=dataset_labels[index],
                intensity_only_dataframes=intensity_only_dataframes,
                full_dataframes=full_dataframes,
                file_path=file_path
            )

    return column_types, analysis_results, data_source, file_paths, intensity_only_dataframes, full_dataframes

def identify_data_source(df: pd.DataFrame) -> str:
    """
    Identify the data source based on the DataFrame columns.

    Args:
    df (pd.DataFrame): Input DataFrame.

    Returns:
    str: Identified data source.
    """
    data_source_indicators: Dict[str, bool] = {
        "maxquant_proteins": ("Only identified by site" in df.columns),
        "maxquant_peptides": ("A Count" in df.columns),
        "alphadia_proteins": ("pg" in df.columns and "base_width_mobility" not in df.columns),
        "alphadia_peptides": ("mod_seq_hash" in df.columns and "base_width_mobility" not in df.columns),
        "alphadia_precursors": ("mod_seq_charge_hash" in df.columns),
        "diann_proteins": ("First.Protein.Description" in df.columns and "Stripped.Sequence" not in df.columns),
        "diann_precursors": ("First.Protein.Description" in df.columns and "Stripped.Sequence" in df.columns),
        "spectronaut_proteins_1": ("PG.ProteinGroups" in df.columns),
        "spectronaut_proteins_2": ("PG.ProteinAccessions" in df.columns),
        "spectronaut_precursors": ("EG.PrecursorId" in df.columns)
    }
    
    for data_source, condition in data_source_indicators.items():
        if condition:
            return data_source
    
    return ""

def analyze_dataframe(
    df: pd.DataFrame,
    data_source: str,
    analysis_results: Dict[str, List[Any]],
    min_data_completeness: float,
    output_directory: str,
    dataset_label: str,
    intensity_only_dataframes: List[pd.DataFrame],
    full_dataframes: List[pd.DataFrame],
    file_path: str
) -> Tuple[List[str], Dict[str, List[Any]], List[pd.DataFrame], List[pd.DataFrame]]:
    """
    Analyze a single DataFrame and update analysis results.

    Args:
    df (pd.DataFrame): Input DataFrame to analyze.
    data_source (str): Identified data source.
    analysis_results (Dict[str, List[Any]]): Current analysis results to update.
    min_data_completeness (float): Minimum data completeness threshold.
    output_directory (str): Directory to save output files.
    dataset_label (str): Label for the current dataset.
    intensity_only_dataframes (List[pd.DataFrame]): List of intensity-only DataFrames to update.
    full_dataframes (List[pd.DataFrame]): List of full DataFrames to update.
    file_path (str): Path of the current file being processed.

    Returns:
    Tuple containing:
    - List of column types
    - Updated dictionary of analysis results
    - Updated list of intensity-only DataFrames
    - Updated list of full DataFrames
    """
    config: Dict[str, Union[str, List[str]]] = get_data_source_config(data_source)
    id_column: str = config["id_column"]
    column_types: List[str] = config["column_types"]
    
    for column_type in column_types:
        df_filtered, intensity_columns = preprocess_columns(df, column_type, data_source)
        
        analysis_results, intensity_only_dataframes, full_dataframes = analyze_and_filter_data(
            df_filtered, analysis_results, column_type, min_data_completeness,
            intensity_columns, id_column, output_directory, data_source,
            dataset_label, intensity_only_dataframes, full_dataframes, file_path
        )
    
    return column_types, analysis_results, intensity_only_dataframes, full_dataframes

def get_data_source_config(data_source: str) -> Dict[str, Union[str, List[str]]]:
    """
    Get configuration for a specific data source.

    Args:
    data_source (str): Identified data source.

    Returns:
    Dict[str, Union[str, List[str]]]: Configuration dictionary for the data source.
    """
    data_source_config: Dict[str, Dict[str, Union[str, List[str]]]] = {
        "maxquant_proteins": {"id_column": "Protein IDs", "column_types": ["Intensity", "iBAQ", "LFQ intensity"]},
        "maxquant_peptides": {"id_column": "Sequence", "column_types": ["Intensity", "LFQ intensity"]},
        "alphadia_proteins": {"id_column": "pg", "column_types": ["Intensity"]},
        "alphadia_peptides": {"id_column": "mod_seq_hash", "column_types": ["Intensity"]},
        "alphadia_precursors": {"id_column": "mod_seq_charge_hash", "column_types": ["Intensity"]},
        "diann_proteins": {"id_column": "Protein.Group", "column_types": ["Intensity"]},
        "diann_precursors": {"id_column": "Precursor.Id", "column_types": ["Intensity"]},
        "spectronaut_proteins_1": {"id_column": "PG.ProteinGroups", "column_types": ["Intensity"]},
        "spectronaut_proteins_2": {"id_column": "PG.ProteinAccessions", "column_types": ["Intensity"]},
        "spectronaut_precursors": {"id_column": "EG.PrecursorId", "column_types": ["Intensity"]}
    }
    return data_source_config[data_source]

def preprocess_columns(df: pd.DataFrame, column_type: str, data_source: str) -> Tuple[pd.DataFrame, List[str]]:
    """
    Preprocess columns based on the data source and column type.

    Args:
    df (pd.DataFrame): Input DataFrame.
    column_type (str): Type of column to process.
    data_source (str): Identified data source.

    Returns:
    Tuple containing:
    - Processed DataFrame
    - List of intensity column names
    """
    return filter_and_rename_columns(df, data_source, column_type)

def select_and_rename_columns(df: pd.DataFrame, selection_criteria: Callable[[str], bool], prefix: str = "") -> Tuple[pd.DataFrame, List[str]]:
    """
    Select and rename columns based on given criteria.

    Args:
    df (pd.DataFrame): Input DataFrame.
    selection_criteria (Callable[[str], bool]): Function to select columns.
    prefix (str, optional): Prefix for new column names. Defaults to "".

    Returns:
    Tuple containing:
    - DataFrame with selected and renamed columns
    - List of new column names
    """
    selected_columns: List[str] = [col for col in df.columns if selection_criteria(col)]
    print(f"Selected columns: {selected_columns}")
    new_column_names: List[str] = [f"{prefix}{i+1}" for i in range(len(selected_columns))]
    df[new_column_names] = df[selected_columns]
    df[new_column_names] = df[new_column_names].replace(0, np.nan)
    print(f"New column names: {new_column_names}")
    return df, new_column_names

def filter_and_rename_columns(df: pd.DataFrame, data_source: str, column_type: str) -> Tuple[pd.DataFrame, List[str]]:
    """
    Filter and rename columns based on the data source and column type.

    Args:
    df (pd.DataFrame): Input DataFrame.
    data_source (str): Identified data source.
    column_type (str): Type of column to process.

    Returns:
    Tuple containing:
    - DataFrame with filtered and renamed columns
    - List of new column names
    """
    if data_source.startswith('maxquant_'):
        selection_criteria = lambda col: column_type in col and col not in ['iBAQ peptides', 'iBAQ', 'Intensity']
        return select_and_rename_columns(df, selection_criteria)
    elif data_source.startswith('alphadia_'):
        selection_criteria = lambda col: not any(x in col for x in ['pg', 'mod_seq_charge_hash', 'mod_seq_hash', 'precursor'])
        return select_and_rename_columns(df, selection_criteria, "Intensity_")
    elif data_source.startswith('diann_'):
        selection_criteria = lambda col: ".d" in col
        return select_and_rename_columns(df, selection_criteria, "Intensity_")
    elif data_source.startswith('spectronaut_'):
        selection_criteria = lambda col: ("[" in col) and ("PG.Quantity" not in col) if "EG.PrecursorId" in df.columns else "[" in col
        df, new_columns = select_and_rename_columns(df, selection_criteria, "Intensity_")
        for col in new_columns:
            df[col] = pd.to_numeric(df[col].replace('Filtered', np.nan), errors='coerce')
        return df, new_columns
    else:
        raise ValueError(f"Unknown data source: {data_source}")

def prepare_intensity_dataframe(df_filtered: pd.DataFrame, intensity_columns: List[str], id_column: str, dataset_label: str) -> pd.DataFrame:
    """
    Prepare intensity DataFrame with unique identifiers.

    Args:
    df_filtered (pd.DataFrame): Filtered input DataFrame.
    intensity_columns (List[str]): List of intensity column names.
    id_column (str): Name of the ID column.
    dataset_label (str): Label for the current dataset.

    Returns:
    pd.DataFrame: Prepared intensity DataFrame.
    """
    intensity_df: pd.DataFrame = df_filtered[intensity_columns].add_prefix(dataset_label + "_")
    intensity_df["unique_id"] = df_filtered[id_column]
    return intensity_df

def analyze_and_filter_data(
    df_filtered: pd.DataFrame,
    analysis_results: Dict[str, List[Any]],
    column_type: str,
    min_data_completeness: float,
    intensity_columns: List[str],
    id_column: str,
    output_directory: str,
    data_source: str,
    dataset_label: str,
    intensity_only_dataframes: List[pd.DataFrame],
    full_dataframes: List[pd.DataFrame],
    file_path: str
) -> Tuple[Dict[str, List[Any]], List[pd.DataFrame], List[pd.DataFrame]]:
    """
    Analyze and filter data, update analysis results, and generate plots.

    Args:
    df_filtered (pd.DataFrame): Filtered input DataFrame.
    analysis_results (Dict[str, List[Any]]): Current analysis results to update.
    column_type (str): Type of column being processed.
    min_data_completeness (float): Minimum data completeness threshold.
    intensity_columns (List[str]): List of intensity column names.
    id_column (str): Name of the ID column.
    output_directory (str): Directory to save output files.
    data_source (str): Identified data source.
    dataset_label (str): Label for the current dataset.
    intensity_only_dataframes (List[pd.DataFrame]): List of intensity-only DataFrames to update.
    full_dataframes (List[pd.DataFrame]): List of full DataFrames to update.
    file_path (str): Path of the current file being processed.

    Returns:
    Tuple containing:
    - Updated dictionary of analysis results
    - Updated list of intensity-only DataFrames
    - Updated list of full DataFrames
    """
    analysis_results, df_filtered = calculate_statistics(
        df_filtered, analysis_results, column_type, min_data_completeness,
        intensity_columns, id_column
    )
    
    try:
        generate_correlation_plots(df_filtered, output_directory, data_source, column_type, dataset_label, intensity_columns)
    except Exception as e:
        print(f"Error processing {file_path}: {str(e)}")
        print("Please ensure experiment names end with 1, 2, 3, 4, ...")

    intensity_df: pd.DataFrame = prepare_intensity_dataframe(df_filtered, intensity_columns, id_column, dataset_label)
    
    intensity_only_dataframes.append(intensity_df)
    full_dataframes.append(df_filtered)

    return analysis_results, intensity_only_dataframes, full_dataframes

def calculate_statistics(
    df_filtered: pd.DataFrame,
    analysis_results: Dict[str, List[Any]],
    column_type: str,
    min_data_completeness: float,
    intensity_columns: List[str],
    id_column: str
) -> Tuple[Dict[str, List[Any]], pd.DataFrame]:
    """
    Calculate statistics for the filtered DataFrame and update analysis results.

    Args:
    df_filtered (pd.DataFrame): Filtered input DataFrame.
    analysis_results (Dict[str, List[Any]]): Current analysis results to update.
    column_type (str): Type of column being processed.
    min_data_completeness (float): Minimum data completeness threshold.
    intensity_columns (List[str]): List of intensity column names.
    id_column (str): Name of the ID column.

    Returns:
    Tuple containing:
    - Updated dictionary of analysis results
    - Filtered and processed DataFrame
    """
    non_numeric: pd.DataFrame = df_filtered[intensity_columns].applymap(lambda x: not pd.api.types.is_numeric_dtype(pd.Series([x])))
    if non_numeric.any().any():
        print("Warning: Non-numeric values found in columns:", intensity_columns)
        print(df_filtered[intensity_columns][non_numeric.any(axis=1)])
    
    df_filtered[f"CV_{column_type}"] = df_filtered[intensity_columns].apply(lambda x: calculate_cv(x), axis=1)
    
    df_filtered[f"data_completeness_{column_type}"] = df_filtered[intensity_columns].notna().sum(axis=1) / len(intensity_columns)
    
    required_columns: List[str] = [id_column, f"CV_{column_type}", f"data_completeness_{column_type}"] + intensity_columns
    df_filtered = df_filtered[required_columns][df_filtered[f"data_completeness_{column_type}"] >= min_data_completeness]
    df_filtered_nonzero: pd.DataFrame = df_filtered[df_filtered[f"data_completeness_{column_type}"] > 0]
    
    analysis_results = update_analysis_results(df_filtered, column_type, intensity_columns, analysis_results)
    
    return analysis_results, df_filtered

def update_analysis_results(df_filtered: pd.DataFrame, column_type: str, intensity_columns: List[str], analysis_results: Dict[str, List[Any]]) -> Dict[str, List[Any]]:
    """
    Update analysis results based on the filtered DataFrame.

    Args:
    df_filtered (pd.DataFrame): Filtered input DataFrame.
    column_type (str): Type of column being processed.
    intensity_columns (List[str]): List of intensity column names.
    analysis_results (Dict[str, List[Any]]): Current analysis results to update.

    Returns:
    Dict[str, List[Any]]: Updated dictionary of analysis results.
    """
    analysis_results["totalIDs"].append(len(df_filtered))
    for cv_threshold in [0.2, 0.1, 0.05]:
        key = f"CV{int(cv_threshold*100)}Percent"
        analysis_results[key].append(len(df_filtered[df_filtered[f"CV_{column_type}"] <= cv_threshold]))

    ids_per_run: List[int] = [len(df_filtered[column].dropna()) for column in intensity_columns]
    analysis_results["IDs_per_run"].append(ids_per_run)
    analysis_results["CV"].append(list(df_filtered[f"CV_{column_type}"]))

    data_completeness_values: List[float] = sorted(set(df_filtered[f"data_completeness_{column_type}"]))
    ids_above_threshold: List[int] = [len(df_filtered[df_filtered[f"data_completeness_{column_type}"] >= value]) for value in data_completeness_values]
    analysis_results["DataCompleteness"].append([data_completeness_values, ids_above_threshold])

    return analysis_results

def calculate_cv(x: pd.Series) -> float:
    """
    Calculate the coefficient of variation for a given series.

    Args:
    x (pd.Series): Input series of numeric values.

    Returns:
    float: Calculated coefficient of variation or NaN if all values are NaN.
    """
    x_numeric: pd.Series = pd.to_numeric(x, errors='coerce')
    if x_numeric.isna().all():
        return np.nan
    return np.nanstd(x_numeric, ddof=1) / np.nanmean(x_numeric)

def generate_correlation_plots(
    df_filtered: pd.DataFrame,
    output_directory: str,
    data_source: str,
    column_type: str,
    dataset_label: str,
    intensity_columns: List[str]
) -> None:
    """
    Generate correlation plots for the filtered DataFrame.

    Args:
    df_filtered (pd.DataFrame): Filtered input DataFrame.
    output_directory (str): Directory to save output files.
    data_source (str): Identified data source.
    column_type (str): Type of column being processed.
    dataset_label (str): Label for the current dataset.
    intensity_columns (List[str]): List of intensity column names.
    """
    common_prefix: str = find_common_prefix(intensity_columns)
    
    plot_sample_correlations(
        df_filtered,
        output_directory,
        data_source,
        column_type,
        dataset_label,
        data_columns=common_prefix + str(list(range(len(intensity_columns) + 1))),
        font_size=12
    )

def find_common_prefix(column_names: List[str]) -> str:
    """
    Find the common prefix among a list of column names.

    Args:
    column_names (List[str]): List of column names.

    Returns:
    str: Common prefix found among the column names.
    """
    if not column_names:
        return ""
    
    reference_name: str = column_names[0]
    for i in range(1, len(column_names)):
        current_name: str = column_names[i]
        match: SequenceMatcher.Match = SequenceMatcher(None, reference_name, current_name).find_longest_match(0, len(reference_name), 0, len(current_name))
        reference_name = reference_name[match.a: match.a + match.size]
    
    return reference_name

def create_correlation_plot(df_subset: pd.DataFrame, plot_type: str, correlation_function: Callable, binning: int, font_size: int) -> go.Figure:
    """ 
    Create a correlation plot based on the specified plot type.

    Args:
    df_subset (pd.DataFrame): Subset of the DataFrame to plot.
    plot_type (str): Type of plot to create ('scatter' or 'heatmap').
    correlation_function (Callable): Function to calculate correlation.
    binning (int): Number of bins for the plot.
    font_size (int): Font size for the plot.

    Returns:
    go.Figure: Created correlation plot.
    """
    if plot_type == "scatter":
        fig = create_scatter_correlation_plot(df_subset, correlation_function, binning)
    elif plot_type == "heatmap":
        fig = create_heatmap_correlation_plot(df_subset, correlation_function)
    else:
        raise ValueError("Invalid plot_type. Choose either 'scatter' or 'heatmap'.")

    fig.update_layout(template="simple_white", coloraxis2=dict(showscale=False, colorscale=["black", "black"]),
                      coloraxis1=dict(showscale=False), font_size=font_size)
    return fig

def create_scatter_correlation_plot(df_subset: pd.DataFrame, correlation_function: Callable, binning: int) -> go.Figure:
    """
    Create a scatter correlation plot.

    Args:
    df_subset (pd.DataFrame): Subset of the DataFrame to plot.
    correlation_function (Callable): Function to calculate correlation.
    binning (int): Number of bins for the plot.

    Returns:
    go.Figure: Created scatter correlation plot.
    """
    fig = make_subplots(rows=len(df_subset.columns), cols=len(df_subset.columns), start_cell='bottom-left',
                        shared_yaxes=True, shared_xaxes=True, horizontal_spacing=0.03, vertical_spacing=0.03)
    value_range = (np.floor(np.nanmin(df_subset)), np.ceil(np.nanmax(df_subset))+1/binning)
    plot_width = plot_height = int((value_range[1]-value_range[0]-1/binning)*binning+1)

    for i, col_i in enumerate(df_subset.columns):
        for j, col_j in enumerate(df_subset.columns):
            canvas = ds.Canvas(plot_width=plot_width, plot_height=plot_height, x_range=value_range, y_range=value_range)
            df_pair = df_subset[[col_i, col_j]].dropna() if i!=j else pd.DataFrame(df_subset[col_i].dropna())
            agg = canvas.points(df_pair, x=col_i, y=col_j)
            agg.values = agg.values.astype(float)
            agg.values[agg.values == 0] = np.nan

            fig.add_trace(go.Heatmap(z=agg, coloraxis="coloraxis1" if i!=j else "coloraxis2"), row=j+1, col=i+1)

            if j == 0:
                fig.update_xaxes(title_text=col_i, row=j+1, col=i+1, tickvals=list(range(0,plot_width,binning)),
                                 ticktext=np.round(agg[col_j].values[0:plot_width:binning]))
            if i == 0:
                fig.update_yaxes(title_text=col_j, row=j+1, col=i+1, tickvals=list(range(0,plot_height,binning)),
                                 ticktext=np.round(agg[col_i].values[0:plot_height:binning]))
            if i!=j:
                fig.add_annotation(dict(text=str(np.round(np.min(correlation_function(df_subset[[col_i,col_j]].dropna())),4)),
                                        x=binning, y=plot_height, showarrow=False), row=j+1, col=i+1)
    return fig

def create_heatmap_correlation_plot(df_subset: pd.DataFrame, correlation_function: Callable) -> go.Figure:
    """
    Create a heatmap correlation plot.

    Args:
    df_subset (pd.DataFrame): Subset of the DataFrame to plot.
    correlation_function (Callable): Function to calculate correlation.

    Returns:
    go.Figure: Created heatmap correlation plot.
    """
    correlation_matrix = np.ones((len(df_subset.columns), len(df_subset.columns)))
    for i, col_i in enumerate(df_subset.columns):
        for j, col_j in enumerate(df_subset.columns):
            if i!=j:
                # print(col_i,col_j)
                # print(df_subset[[col_i,col_j]].head())
                # print(correlation_function(df_subset[[col_i,col_j]].dropna()))
                correlation_matrix[i,j] = np.round(np.min(correlation_function(df_subset[[col_i,col_j]].dropna())),4)
    fig = go.Figure(data=go.Heatmap(z=correlation_matrix))
    fig.update_xaxes(tickvals=list(range(len(df_subset.columns))), ticktext=list(df_subset.columns))
    fig.update_yaxes(tickvals=list(range(len(df_subset.columns))), ticktext=list(df_subset.columns))
    return fig


def plot_sample_correlations(
    df: pd.DataFrame,
    output_directory: str,
    data_source: str,
    column_type: str,
    dataset_label: str = "",
    data_columns: str = "Intensity (.*)",
    correlation_function: Callable[[pd.DataFrame], np.ndarray] = lambda x: np.corrcoef(x.T),
    plot_type: str = "scatter",
    log_transform: bool = True,
    binning: int = 10,
    font_size: int = 10,
) -> go.Figure:
    """
    Plot sample correlations and save the resulting figure.

    Args:
    df (pd.DataFrame): Input DataFrame.
    output_directory (str): Directory to save output files.
    data_source (str): Identified data source.
    column_type (str): Type of column being processed.
    dataset_label (str, optional): Label for the current dataset. Defaults to "".
    data_columns (str, optional): Regex pattern to select data columns. Defaults to "Intensity (.*)".
    correlation_function (Callable, optional): Function to calculate correlation. Defaults to Pearson correlation.
    plot_type (str, optional): Type of plot to create ('scatter' or 'heatmap'). Defaults to "scatter".
    log_transform (bool, optional): Whether to apply log transformation. Defaults to True.
    binning (int, optional): Number of bins for the plot. Defaults to 10.
    font_size (int, optional): Font size for the plot. Defaults to 10.

    Returns:
    go.Figure: Created correlation plot.
    """
    df_subset = df[[col for col in df.columns if re.match(data_columns, col)]].copy()
    if log_transform:
        df_subset = df_subset.apply(np.log10)
    df_subset = df_subset.replace([np.inf, -np.inf], np.nan)
    df_subset.columns = [re.findall(data_columns, col)[0] for col in df_subset.columns]

    fig = create_correlation_plot(df_subset, plot_type, correlation_function, binning, font_size)
    
    if plot_type == "scatter":
        fig.update_layout(width=len(df_subset.columns)*200+100, height=len(df_subset.columns)*200+50, title=dataset_label)
        fig.write_image(f"{output_directory}/correlation_plots_pearson_{dataset_label}_{data_source}_{column_type}.pdf")
    else:
        fig.update_layout(width=500, height=500, title="Pearson Correlation")
        fig.write_image(f"{output_directory}/correlation_plots_heatmap_{data_source}_{column_type}.pdf")
    
    fig.show()
    return fig






# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------









def make_plots(
    column_indicators,
    dict_results,
    labels,
    save_output_here, 
    file_type,
    file_names,
    list_dfs_only_with_intensity_columns,
    list_dfs
    ):

    for column_indicator in column_indicators:

        # bar plots with CVs
        df_results = pd.DataFrame(dict_results)
        index_column_indicator = column_indicators.index(column_indicator)

        list_filtered_for_column_indicator = list()
        for file_num in range(len(file_names)):
            list_filtered_for_column_indicator.append(index_column_indicator + len(column_indicators) * file_num)

        df_results_filtered_for_column_indicator = df_results.iloc[list_filtered_for_column_indicator]
        print(df_results_filtered_for_column_indicator)

        # print([x for x in df_results_filtered_for_column_indicator["CV"].iloc[0] if not math.isnan(x)])

        print(np.median([x for x in df_results_filtered_for_column_indicator["CV"].iloc[0] if not math.isnan(x)]),
              np.median([x for x in df_results_filtered_for_column_indicator["CV"].iloc[1] if not math.isnan(x)]),
            #   np.median(df_results_filtered_for_column_indicator["CV"].iloc[2]),
             )

        make_bar_bar_plot(
            df_results_filtered_for_column_indicator,
            labels,
            save_output_here, 
            file_type,
            column_indicator,
        )

        print("The following plot requires an equal number of replicates in all analysis output files")
        bar_missing_value(
            df_results_filtered_for_column_indicator,
            labels,
            save_output_here,
            file_type,
            column_indicator    
        )

        # bar plot with error bar per run and average ID rates
        plot_IDs_per_run(
            df_results_filtered_for_column_indicator,
            labels,
            save_output_here,
            file_type,
            column_indicator,
            dotsize = 10,
            capsize = 22
        )

        # Violin plots of CVs
        generate_violin_plots(
            df_results_filtered_for_column_indicator,
            labels,
            save_output_here,
            file_type,
            column_indicator
        )

        # correlation heatmap
        dfs_intensity_for_column_indicator = pd.DataFrame()
        dfs_intensity_for_column_indicator["unique_id"] = []
        for df_temp in list_dfs_only_with_intensity_columns:
            column_list = [col for col in df_temp.columns if column_indicator in col]
            column_list.append("unique_id")
            df_selected = df_temp[column_list]
            df_selected['unique_id'] = df_selected['unique_id'].astype(str)

            dfs_intensity_for_column_indicator = dfs_intensity_for_column_indicator.merge(
                df_selected, 
                on="unique_id", 
                how='outer'
            )

        # function form Julia P. Schessner, Eugenia Voytik, Isabell Bludau,  https://doi.org/10.1002/pmic.202100103
        cross_corr2 = plot_sample_correlations(
            dfs_intensity_for_column_indicator,
            save_output_here,
            file_type,
            column_indicator,
            data_columns="(.*)"+column_indicator+"(.*)", 
            plot_type="heatmap",  # Changed from mode to plot_type
            font_size=12
        )

    #     # abundance range illustrated as rank plot
    #     rank_plot(
    #         list_dfs, 
    #         labels,
    #         column_indicator,
    #         save_output_here,
    #         file_type,
    #     )

    #     # abundance range illustrated as rank plot
    #     rank_plot_static(
    #         list_dfs, 
    #         labels,
    #         column_indicator,
    #         save_output_here,
    #         file_type,
    #     )

    # #     # CV values illustrated as rank plot
    #     list_scatter_plot_index = rank_plot_static(
    #         list_dfs, 
    #         labels,
    #         column_indicator,
    #         save_output_here,
    #         file_type,
    #         x_value="index",
    #         y_value= "CV_"+column_indicator,
    #     )

    #     scatter_density_plot(
    #         list_scatter_plot_index,
    #         labels,
    #         column_indicator,
    #         save_output_here,
    #         file_type,
    #         x_value="index",
    #         y_value= "CV_"+column_indicator,
    # #         kde=True
    #     )

    #     # missing value illustrated as rank plot
    #     list_scatter_plot_missing_index = rank_plot_static(
    #         list_dfs, 
    #         labels,
    #         column_indicator,
    #         save_output_here,
    #         file_type,
    #         x_value="index",
    #         y_value= "data_completeness_"+column_indicator,
    #     )
    #     scatter_density_plot(
    #         list_scatter_plot_missing_index,
    #         labels,
    #         column_indicator,
    #         save_output_here,
    #         file_type,
    #         x_value="index",
    #         y_value= "data_completeness_"+column_indicator,
    #     )

    #     # CV illustrated along log2_median
    #     list_scatter_plot_log2 = rank_plot_static(
    #         list_dfs, 
    #         labels,
    #         column_indicator,
    #         save_output_here,
    #         file_type,
    #         x_value="log2_median",
    #         y_value= "CV_"+column_indicator,
    #     )
    #     scatter_density_plot(
    #         list_scatter_plot_log2,
    #         labels,
    #         column_indicator,
    #         save_output_here,
    #         file_type,
    #         x_value="log2_median",
    #         y_value= "CV_"+column_indicator,
    # #         kde=True
    #     )

    #     # missing value illustrated along log2_median
    #     list_scatter_plot_missing_log2 = rank_plot_static(
    #         list_dfs, 
    #         labels,
    #         column_indicator,
    #         save_output_here,
    #         file_type,
    #         x_value="log2_median",
    #         y_value= "data_completeness_"+column_indicator,
    #     )
    #     scatter_density_plot(
    #         list_scatter_plot_missing_log2,
    #         labels,
    #         column_indicator,
    #         save_output_here,
    #         file_type,
    #         x_value="log2_median",
    #         y_value= "data_completeness_"+column_indicator,
    #     )

    #     # Venn diagram or upset plot to show shared IDs
    #     make_venn_or_upset_plot(
    #         list_dfs_only_with_intensity_columns,
    #         column_indicator,
    #         save_output_here
    #     )

def bar_plot(
    df,
    names,
    save_name,
    ylabel,
    xlabel,
    titel,
    label_box
):

    # colors = ['lightgrey', '#7AC7C9','#4EA7BB','#267FA5', 'lightgrey']
    colors = plt.cm.viridis(np.linspace(0,1,len(label_box)+1))

    r = list(range(0,len(df)))

    barWidth = 0.5

    # Create brown bars
    ctr=0
    for m in list(df.columns):
        bars =list(df[m].values)

        # Create brown bars
        plt.bar(r, bars, label=label_box[ctr], color=colors[ctr], edgecolor='white', width=barWidth,alpha=0.4)
        ctr+=1

    plt.xticks(r, names, fontweight='bold')##rotation='vertical')
    #plt.legend(loc="lower right")
    plt.legend(loc="lower right", bbox_to_anchor=(1.55, 0.75), ncol=1)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.xticks(rotation=90)
    plt.title(titel)
    plt.savefig(save_name+".pdf", bbox_inches='tight', pad_inches=0,dpi=300)
    plt.savefig(save_name+".png", bbox_inches='tight', pad_inches=0,dpi=300)
    plt.show()
    plt.clf()

def find_y_label(file_type):
    dict_ylabels = {
        "mq_proteins": "# protein groups",
        'mq_peptides': "# peptides",
        'alphadia_proteins': "# protein groups",
        "alphadia_peptides": "# peptides",
        "alphadia_precursors": "# precursors",
        'diann_proteins': "# protein groups",
        'diann_precursors': "# precursors",
        'spectronaut_proteins': "# protein groups",
        'spectronaut_precursors': "# precursors",
    }
    return dict_ylabels[file_type]

def make_bar_bar_plot(
    df_results_filtered_for_column_indicator,
    labels,
    save_output_here,
    file_type,
    column_indicator,
):

    save_name = save_output_here + "/bar_plot_CVs_"+ file_type + "_" + column_indicator
    ylabel = find_y_label(file_type) # Name of the y axis for the bar plot
    xlabel = 'Type' # Name or topic of all the columns
    label_box = ['Total No.', '<= CV 20%', '<= CV 10%', '<= CV 5%']


    bar_plot(
        df_results_filtered_for_column_indicator[['totalIDs', 'CV20Percent', 'CV10Percent', 'CV5Percent']],
        labels,
        save_name,
        ylabel,
        xlabel,
        column_indicator,
        label_box
    )

def plot_IDs_per_run(
    df_results_filtered_for_column_indicator,
    labels,
    save_output_here,
    file_type,
    column_indicator,
    dotsize = 15,
    capsize = 25
):

    mean_list_pep = list()
    error_pep = list()
    list_dfs = list()
    x_pos_pep = range(len(labels))

    for x in x_pos_pep:
        # calculate mean and std for each sample
        pep_temp = df_results_filtered_for_column_indicator["IDs_per_run"].iloc[x]
        mean_list_pep.append(np.mean(pep_temp))

        error_pep.append(np.std(pep_temp))

        #create joined data frame
        jitter = pd.DataFrame(pep_temp)
        jitter["Type"] = labels[x]
        list_dfs.append(jitter)
    print(mean_list_pep)
    jitter = pd.concat(list_dfs)
    jitter.columns =['value', 'Type']

    # Build the plot
    fig, ax = plt.subplots()
    ax.bar(x_pos_pep, mean_list_pep, yerr=error_pep, align='center', ecolor='black', capsize=capsize, color = '#3777ac', label="average")
    ax = sns.swarmplot(x='Type', y='value', data=jitter, color = "grey", alpha = 0.8, size = dotsize)
    ax.set_ylabel(find_y_label(file_type))
    ax.set_xticks(x_pos_pep)
    ax.set_xticklabels(labels)
    plt.title(column_indicator)
    plt.legend(loc="lower right", bbox_to_anchor=(1.4, 0.85), ncol=1)

    # Save the figure and show
    plt.savefig(save_output_here + "/bar_plot_IDs_per_run_"+ file_type + "_" + column_indicator + ".png", bbox_inches='tight', pad_inches=0,dpi=300)
    plt.savefig(save_output_here + "/bar_plot_IDs_per_run_"+ file_type + "_" + column_indicator + ".pdf", bbox_inches='tight', pad_inches=0,dpi=300)

    plt.show()

def generate_violin_plots(
    df_results_filtered_for_column_indicator,
    labels,
    save_output_here,
    file_type,
    column_indicator
):

    list_dfs = list()
    x_pos_pep = range(len(labels))

    for x in x_pos_pep:
        #create joined data frame
        df_temp = pd.DataFrame(df_results_filtered_for_column_indicator["CV"].iloc[x])
        df_temp["Type"] = labels[x]
        list_dfs.append(df_temp)

    df = pd.concat(list_dfs)
    df.columns =['value', 'Type']

    # Build the plot
    fig, ax = plt.subplots()
    sns.violinplot(data=df.drop_duplicates(), x='Type', y='value',palette="Blues")
    ax.set_ylabel("CV (" + find_y_label(file_type)[2:] + ")")
    plt.title(column_indicator)

    # Save the figure and show
    plt.savefig(save_output_here + "/violin_plot_IDs_per_run_"+ file_type + "_" + column_indicator+".pdf", bbox_inches='tight', pad_inches=0,dpi=300)
    plt.savefig(save_output_here + "/violin_plot_IDs_per_run_"+ file_type + "_" + column_indicator+".png", bbox_inches='tight', pad_inches=0,dpi=300)

    plt.show()


def bar_missing_value(
    df_results_filtered_for_column_indicator,
    labels,
    save_output_here,
    file_type,
    column_indicator
):
    list_dfs_temp = list()
    x_pos_pep = range(len(labels))

    # print(df_results_filtered_for_column_indicator)
    for x in x_pos_pep:
        #create joined data frame
        list_filtered = df_results_filtered_for_column_indicator["DataCompleteness"].iloc[x]
        df_temp = pd.DataFrame(np.array([list_filtered[1]]),
                               columns=list_filtered[0])
        df_temp["Type"] = labels[x]
        list_dfs_temp.append(df_temp)

    df = pd.concat(list_dfs_temp).reset_index(drop=True).fillna(0)
    print(df)
    columns_selected = sorted([column for column in df.columns if type(column) != str])

    save_name = save_output_here + "/bar_plot_DataCompleteness_"+ file_type + "_" + column_indicator
    ylabel = find_y_label(file_type) # Name of the y axis for the bar plot
    xlabel = 'Type' # Name or topic of all the columns
    label_box = columns_selected*len(df_results_filtered_for_column_indicator["IDs_per_run"].iloc[0])

    bar_plot(
        df[columns_selected],
        labels,
        save_name,
        ylabel,
        xlabel,
        column_indicator,
        label_box
    )

def rank_plot(
    list_dfs,
    labels,
    column_indicator,
    save_output_here,
    file_type,
    x_value="index",
    y_value= "log10_median",
):
    colors = sample_colorscale('blues', len(labels)+1)
    fig1 = px.scatter()
    num = 0
    for dfs_intensity in list_dfs:
        if True in dfs_intensity.columns.str.contains(column_indicator+"_"):
            label = labels[num]
            selected_columns = [col for col in dfs_intensity.columns if column_indicator in col]
            dfs_intensity['median'] = dfs_intensity[selected_columns].median(axis=1)
            dfs_intensity["log10_median"] = np.log10(dfs_intensity["median"])
            dfs_intensity_sorted = dfs_intensity.sort_values(by="log10_median", ascending=False)
            dfs_intensity_index = dfs_intensity_sorted.reset_index()
            dfs_intensity_index["index"] = dfs_intensity_index.index

            fig1.add_scatter(x=dfs_intensity_index[x_value], y=dfs_intensity_index[y_value], mode="markers", marker=dict(color=colors[num+1]), name=label)

            dict_labels = {
                "index": "Rank", "log10_median": "Median Intensities (log10)",
                "CV_"+column_indicator: "CV", "data_completeness_"+column_indicator: "missing value ratio"
            }
            x_axis_label = dict_labels[x_value]
            y_axis_label = dict_labels[y_value]
            fig1.update_layout(
                title=dict(
                    text=x_axis_label+" versus "+y_axis_label,
                    font=dict(
                        size=16,
                    ),
                    x=0.5,
                    xanchor='center',
                    yanchor='top',
                    y=0.92
                ),
                xaxis=dict(
                    title=x_axis_label,
                    titlefont_size=14,
                    #tickmode='auto',
                    tickfont_size=14,
                    #autorange=False,
                ),
                yaxis=dict(
                    title=y_axis_label,
                    # type="log"
                ),
                legend_title_text='Legend:',
                #hovermode="x",
                template='plotly_white',
                # width=1000,
                height=450
            )

            if (x_value=="index") & (y_value== "log10_median"):
                # caclulate the position of the quartiles
                print(label+":", "quartile, position on the y-axis:")
                df_quartile = dfs_intensity["log10_median"].quantile([0.25,0.5,0.75])
                for quartile in [0.25, 0.5, 0.75]:
                    print(quartile, df_quartile[quartile])
        #             fig1.add_hline(y=df_quartile[quartile], line_color=colors[num+1], annotation_text=str(quartile), annotation_position="bottom right",annotation_font_color=colors[num+1], opacity=0.25, legend='legend2')


                # calculate the number of precursors per 25% of the dynamic range
                diff_max_min = (dfs_intensity["log10_median"].max() - dfs_intensity["log10_median"].min()) / 4
                percentile_ranges = [[dfs_intensity["log10_median"].min(), dfs_intensity["log10_median"].min()+diff_max_min],
                [dfs_intensity["log10_median"].min()+diff_max_min, dfs_intensity["log10_median"].min()+diff_max_min*2],
                [dfs_intensity["log10_median"].min()+diff_max_min*2, dfs_intensity["log10_median"].min()+diff_max_min*3],
                [dfs_intensity["log10_median"].min()+diff_max_min*3, dfs_intensity["log10_median"].max()]]

                print("number of precursors per 25% of the dynamic range")
                num2 = 0
                for percentile in percentile_ranges:
                    precursors_within_range = len(dfs_intensity[(dfs_intensity["log10_median"]>=percentile[0]) & (dfs_intensity["log10_median"]<percentile[1])])
                    if num != len(percentile_ranges)-1:
                        print(str(percentile) + "(25%):")
                        print("count:", precursors_within_range, ", ratio from total IDs:", round(precursors_within_range/len(dfs_intensity["log10_median"]),2))
                    else:
                        print(str(percentile) + "(25%):")
                        print("count:", precursors_within_range, ", ratio from total IDs:", round(precursors_within_range/len(dfs_intensity["log10_median"]),2))
                    num2+=1
            num += 1
        else:
            next
    fig1.show()
    fig1.write_image(save_output_here + "/plot_"+x_axis_label+"_versus_"+y_axis_label+ "_" + file_type + "_" + column_indicator + ".pdf")
    return fig1


def rank_plot_static(
    list_dfs,
    labels,
    column_indicator,
    save_output_here,
    file_type,
    x_value="index",
    y_value= "log10_median",
):
    colors = plt.cm.viridis(np.linspace(0,1,len(labels)+1))
    # colors = ["#E48E25", "#4E3468", "#E48E25",]
    num = 0
    list_scatter_plot = list()
    # fig, ax = plt.subplots()
    for dfs_intensity in list_dfs:
        if True in dfs_intensity.columns.str.contains(column_indicator+"_"):
            label = labels[num]
            selected_columns = [col for col in dfs_intensity.columns if column_indicator in col]
            dfs_intensity['median'] = dfs_intensity[selected_columns].median(axis=1)
            dfs_intensity["log10_median"] = np.log10(dfs_intensity["median"])
            dfs_intensity_sorted = dfs_intensity.sort_values(by="log10_median", ascending=False)
            dfs_intensity_index = dfs_intensity_sorted.reset_index()
            dfs_intensity_index["index"] = dfs_intensity_index.index

            if "log10" in y_value:
                plt.plot(x_value, y_value, data = dfs_intensity_index, color = colors[num+1], linewidth = 2, label = label) #marker=".",linewidth = 0, label = label)
            else:
                dfs_intensity["log2_median"] = np.log2(dfs_intensity["median"])
                dfs_intensity_sorted = dfs_intensity.sort_values(by="log2_median", ascending=False)
                dfs_intensity_index = dfs_intensity_sorted.reset_index()
                dfs_intensity_index["index"] = dfs_intensity_index.index
                df_temp = pd.DataFrame()
                df_temp[x_value] = dfs_intensity_index[x_value]
                df_temp[y_value] = dfs_intensity_index[y_value]
                df_temp["Type: "] = label
                list_scatter_plot.append(df_temp)

            num += 1
        else:
            next
    dict_labels = {
        "index": "Rank", "log10_median": "Median Intensities (log10)", "log2_median": "Median Intensities (log2)",
        "CV_"+column_indicator: "CV", "data_completeness_"+column_indicator: "missing value ratio"
    }
    x_axis_label = dict_labels[x_value]
    y_axis_label = dict_labels[y_value]

    if "log10" in y_value:
        plt.xlabel(x_axis_label)
        plt.ylabel(y_axis_label)
        plt.title(x_axis_label+" versus "+y_axis_label)
        plt.legend()
        plt.savefig(save_output_here + "/plot_"+x_axis_label+"_versus_"+y_axis_label+ "_" + file_type + "_" + column_indicator + ".pdf",bbox_inches='tight', pad_inches=0, dpi=300)
        plt.savefig(save_output_here + "/plot_"+x_axis_label+"_versus_"+y_axis_label+ "_" + file_type + "_" + column_indicator + ".png",bbox_inches='tight', pad_inches=0, dpi=300)
        plt.show()
        plt.clf()
    else:
        return list_scatter_plot

def scatter_density_plot(
    list_scatter_plot,
    labels,
    column_indicator,
    save_output_here,
    file_type,
    x_value="index",
    y_value= "log10_median",
    kde=False,
):
    dict_labels = {
        "index": "Rank", "log10_median": "Median Intensities (log10)", "log2_median": "Median Intensities (log2)",
        "CV_"+column_indicator: "CV", "data_completeness_"+column_indicator: "missing value ratio"
    }
    x_axis_label = dict_labels[x_value]
    y_axis_label = dict_labels[y_value]

    df_scatter_plot = pd.concat(list_scatter_plot).reset_index(drop=True)
    # Generate a scatter plot comparing sample A and B
    print("Please wait - it may take several minutes")
    # delete kind="kde" if you don't want to wait
    colors = sns.color_palette('viridis', len(labels)+1)
    fig = sns.jointplot(data=df_scatter_plot, x=x_value, y=y_value, hue="Type: ", palette="viridis", s=5, ec="face")#, kind="kde")
    plt.xlabel(x_axis_label)
    plt.ylabel(y_axis_label)
    fig.savefig(save_output_here + "/plot_"+x_axis_label+"_versus_"+y_axis_label+ "_" + file_type + "_" + column_indicator + "_wo_kde.pdf",bbox_inches='tight', pad_inches=0, dpi=300)
    plt.show()
    if kde==True:
        fig = sns.jointplot(data=df_scatter_plot, x=x_value, y=y_value, hue="Type: ", palette="viridis", kind="kde")
        plt.xlabel(x_axis_label)
        plt.ylabel(y_axis_label)
        fig.savefig(save_output_here + "/plot_"+x_axis_label+"_versus_"+y_axis_label+ "_" + file_type + "_" + column_indicator + "_w_kde.pdf",bbox_inches='tight', pad_inches=0, dpi=300)
        plt.show()


def make_venn_or_upset_plot(list_dfs_only_with_intensity_columns, column_indicator, save_output_here):

    set_all = set()
    dict_set = dict()

    for df_temp in list_dfs_only_with_intensity_columns:
        if True in df_temp.columns.str.contains(column_indicator):
            label = df_temp.columns[0].split(column_indicator)[0]
            set_temp = set(df_temp["unique_id"])
            set_all = set_all.union(set_temp)
            dict_set[label] = set_temp
        else:
            next

    if len(dict_set) == 2:
        print("Make Venn diagram")

        make_venn2_diagram(dict_set, save_output_here, column_indicator)

    elif len(dict_set) == 3:
        print("Make Venn diagram")

        make_venn3_diagram(dict_set, save_output_here, column_indicator)

    else:
        print ("More than three files: Make upset plot")

        make_upset_plot(dict_set, set_all, save_output_here, column_indicator)

def make_venn3_diagram(
    dict_set,
    save_output_here,
    column_indicator
):
    viridis = cm.get_cmap('viridis', 3)
    colors = viridis.colors

    venn3(
        [value for key, value in dict_set.items()],
        set_labels = [key for key, value in dict_set.items()],
        set_colors=colors, alpha = 0.7
    )
    plt.title(column_indicator)
    # Save the figure and show
    plt.savefig(save_output_here + "/venn_diagram_"+ column_indicator + ".png", bbox_inches='tight', pad_inches=0,dpi=300)
    plt.savefig(save_output_here + "/venn_diagram_"+ column_indicator + ".pdf", bbox_inches='tight', pad_inches=0,dpi=300)
    plt.show()

def make_venn2_diagram(
    dict_set,
    save_output_here,
    column_indicator
):
    viridis = cm.get_cmap('viridis', 2)
    colors = viridis.colors

    venn2(
        [value for key, value in dict_set.items()],
        set_labels = [key for key, value in dict_set.items()],
        set_colors=colors, alpha = 0.7
    )
    plt.title(column_indicator)
    # Save the figure and show
    plt.savefig(save_output_here + "/venn_diagram_"+ column_indicator + ".png", bbox_inches='tight', pad_inches=0,dpi=300)
    plt.savefig(save_output_here + "/venn_diagram_"+ column_indicator + ".pdf", bbox_inches='tight', pad_inches=0,dpi=300)
    plt.show()

def make_upset_plot(dict_set, set_all, save_output_here, column_indicator):
    dict_df = dict()
    for key, value in dict_set.items():
        list_temp = [True if e in value else False for e in set_all]
        dict_df[key] = list_temp

    df = pd.DataFrame(dict_df)
    df_up = df.groupby(list(dict_set.keys())).size()
    upset_plot(df_up, orientation='horizontal')
    plt.suptitle(column_indicator)
    # Save the figure and show
    plt.savefig(save_output_here + "/upset_plot_"+ column_indicator + ".png", bbox_inches='tight', pad_inches=0,dpi=300)
    plt.savefig(save_output_here + "/upset_plot_"+ column_indicator + ".pdf", bbox_inches='tight', pad_inches=0,dpi=300)
    plt.show()



