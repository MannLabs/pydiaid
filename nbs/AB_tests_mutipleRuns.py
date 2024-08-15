# Author: Patricia Skowronek, Max Planck Institute of Biochemistry

# Standard library imports
import math
import re
import os
from difflib import SequenceMatcher
from typing import List, Dict, Any, Tuple, Callable, Union, Set, Optional

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
    dataset_labels: List[str],
    filter_columns: List[str] = None
) -> Tuple[List[str], Dict[str, List[Any]], str, List[str], List[pd.DataFrame], List[pd.DataFrame]]:
    """
    Process proteomics data from multiple files and perform analysis.

    This function reads proteomics data from the provided file paths, identifies the data source,
    analyzes the data, and compiles the results. It can optionally filter the data based on
    specified columns.

    Args:
        file_paths (List[str]): A list of file paths to the proteomics data files.
        min_data_completeness (float): The minimum data completeness threshold for filtering.
        output_directory (str): The directory where output files will be saved.
        dataset_labels (List[str]): Labels for each dataset corresponding to the file paths.
        filter_columns (List[str], optional): Columns to use for filtering the data. Defaults to None.

    Returns:
        Tuple[List[str], Dict[str, List[Any]], str, List[str], List[pd.DataFrame], List[pd.DataFrame]]:
            - List[str]: Types of columns identified in the data.
            - Dict[str, List[Any]]: Analysis results containing various metrics.
            - str: Identified data source type.
            - List[str]: List of processed file paths.
            - List[pd.DataFrame]: List of DataFrames containing only intensity data.
            - List[pd.DataFrame]: List of full DataFrames with all processed data.

    The analysis_results dictionary contains the following keys:
        - "totalIDs": Total number of identifications.
        - "CV20Percent": Number of identifications with CV <= 20%.
        - "CV10Percent": Number of identifications with CV <= 10%.
        - "CV5Percent": Number of identifications with CV <= 5%.
        - "IDs_per_run": Number of identifications per run.
        - "CV": Coefficient of variation values.
        - "DataCompleteness": Data completeness metrics.

    Each file is processed individually, and the results are aggregated. If filter_columns
    are provided, the analysis is performed for each filter condition separately.
    """
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
            if filter_columns:
                for filter_condition in filter_columns:
                    column_types, analysis_results, intensity_only_dataframes, full_dataframes = analyze_dataframe(
                        df=df,
                        analysis_results=analysis_results,
                        min_data_completeness=min_data_completeness,
                        output_directory=output_directory,
                        data_source=data_source,
                        dataset_label=f"{dataset_labels[index]}_{filter_condition}",
                        intensity_only_dataframes=intensity_only_dataframes,
                        full_dataframes=full_dataframes,
                        file_path=file_path,
                        filter_columns=filter_condition
                    )
            else:
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
    file_path: str,
    filter_columns: str = None  
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
    filter_columns (str): str to filter column names for.

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
        df_filtered, intensity_columns = preprocess_columns(df, column_type, data_source, filter_columns) 
        
        analysis_results, intensity_only_dataframes, full_dataframes = analyze_and_filter_data(
            df_filtered, analysis_results, column_type, min_data_completeness,
            intensity_columns, id_column, output_directory, data_source,
            dataset_label, intensity_only_dataframes, full_dataframes, file_path,
            filter_columns 
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

def preprocess_columns(df: pd.DataFrame, column_type: str, data_source: str, filter_columns: str = None) -> Tuple[pd.DataFrame, List[str]]:
    """
    Preprocess columns based on the data source and column type.

    Args:
    df (pd.DataFrame): Input DataFrame.
    column_type (str): Type of column to process.
    data_source (str): Identified data source.
    filter_columns (str): str to filter column names for.

    Returns:
    Tuple containing:
    - Processed DataFrame
    - List of intensity column names
    """
    return filter_and_rename_columns(df, data_source, column_type, filter_columns)

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

def filter_and_rename_columns(df: pd.DataFrame, data_source: str, column_type: str, filter_columns: str = None) -> Tuple[pd.DataFrame, List[str]]:
    """
    Filter and rename columns based on the data source and column type.

    Args:
    df (pd.DataFrame): Input DataFrame.
    data_source (str): Identified data source.
    column_type (str): Type of column to process.
    filter_columns (str): str to filter column names for.

    Returns:
    Tuple containing:
    - DataFrame with filtered and renamed columns
    - List of new column names
    """
    if data_source.startswith('maxquant_'):
        selection_criteria = lambda col: column_type in col and col not in ['iBAQ peptides', 'iBAQ', 'Intensity'] and (filter_columns is None or filter_columns in col)
        return select_and_rename_columns(df, selection_criteria)
    elif data_source.startswith('alphadia_'):
        selection_criteria = lambda col: not any(x in col for x in ['pg', 'mod_seq_charge_hash', 'mod_seq_hash', 'precursor']) and (filter_columns is None or filter_columns in col)
        return select_and_rename_columns(df, selection_criteria, "Intensity_")
    elif data_source.startswith('diann_'):
        selection_criteria = lambda col: ".d" in col and (filter_columns is None or filter_columns in col)
        return select_and_rename_columns(df, selection_criteria, "Intensity_")
    elif data_source.startswith('spectronaut_'):
        selection_criteria = lambda col: ("[" in col) and ("PG.Quantity" not in col) and (filter_columns is None or filter_columns in col) if "EG.PrecursorId" in df.columns else "[" in col and (filter_columns is None or filter_columns in col)
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
    file_path: str,
    filter_columns: str = None 
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
    filter_columns (str): str to filter column names for.

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
    intensity_columns: List[str],
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
        font_size=12,
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




def generate_all_plots(
    data_column_identifiers: List[str],
    results_dictionary: Dict[str, Any],
    labels: List[str],
    output_directory: str,
    data_file_type: str,
    file_names: List[str],
    intensity_data_frames: List[pd.DataFrame],
    all_data_frames: List[pd.DataFrame],
    filter_columns: List[str] = None
) -> None:
    """
    Generate all plots for the given data.

    Args:
        data_column_identifiers (List[str]): List of column identifiers.
        results_dictionary (Dict[str, Any]): Dictionary containing results data.
        labels (List[str]): List of labels for the data.
        output_directory (str): Directory to save output files.
        data_file_type (str): Type of data file being analyzed.
        file_names (List[str]): List of file names.
        intensity_data_frames (List[pd.DataFrame]): List of DataFrames with intensity data.
        all_data_frames (List[pd.DataFrame]): List of all DataFrames.
        filter_columns (List[str], optional): List of columns to filter by. Defaults to None.

    Returns:
        None
    """
    os.makedirs(output_directory, exist_ok=True)

    for column_identifier in data_column_identifiers:
        labels2: List[str] = [f"{label}_{filter_col}" for label in labels for filter_col in (filter_columns or [])] or labels

        results_df: pd.DataFrame = pd.DataFrame(results_dictionary)
        print(labels2, results_df)
        column_identifier_index: int = data_column_identifiers.index(column_identifier)

        filtered_indices: List[int] = _get_filtered_indices(
            column_identifier_index, len(data_column_identifiers), len(file_names), filter_columns
        )

        filtered_results_df: pd.DataFrame = results_df.iloc[filtered_indices]

        # Calculate and print median CV values
        print_median_cv_values(filtered_results_df, labels2)

        generate_and_save_individual_plots(
            filtered_results_df, labels2, output_directory, data_file_type, column_identifier, 
            intensity_data_frames, all_data_frames
        )

def print_median_cv_values(filtered_results_df: pd.DataFrame, labels: List[str]) -> None:
    """
    Print median CV values for each sample in the filtered results DataFrame.

    Args:
        filtered_results_df (pd.DataFrame): DataFrame containing filtered results.
        labels (List[str]): List of labels for each sample.

    Returns:
        None
    """
    print("Median CV values:")
    for i, label in enumerate(labels):
        if i < len(filtered_results_df["CV"]):
            cv_values = filtered_results_df["CV"].iloc[i]
            if isinstance(cv_values, list):
                median_cv = np.median([x for x in cv_values if not np.isnan(x)])
                print(f"Sample {label}: {median_cv}")
            else:
                print(f"Sample {label}: CV data is not in the expected list format")
        else:
            print(f"Sample {label}: No CV data available")
    print()  # Add a blank line for better readability

def generate_and_save_individual_plots(
    filtered_results_df: pd.DataFrame,
    labels: List[str],
    output_directory: str,
    data_file_type: str,
    column_identifier: str,
    intensity_data_frames: List[pd.DataFrame],
    all_data_frames: List[pd.DataFrame]
) -> None:
    """
    Generate and save individual plots for the given data.

    Args:
        filtered_results_df (pd.DataFrame): DataFrame containing filtered results.
        labels (List[str]): List of labels for the data.
        output_directory (str): Directory to save output files.
        data_file_type (str): Type of data file being analyzed.
        column_identifier (str): Identifier for the current column.
        intensity_data_frames (List[pd.DataFrame]): List of DataFrames with intensity data.
        all_data_frames (List[pd.DataFrame]): List of all DataFrames.

    Returns:
        None
    """
    create_cv_bar_plot(filtered_results_df, labels, output_directory, data_file_type, column_identifier)
    plot_identifications_per_run(filtered_results_df, labels, output_directory, data_file_type, column_identifier)
    plot_missing_value_bars(filtered_results_df, labels, output_directory, data_file_type, column_identifier)
    generate_violin_plots(filtered_results_df, labels, output_directory, data_file_type, column_identifier)
    
    if intensity_data_frames:
        combined_df = merge_intensity_dataframes(intensity_data_frames)
        plot_sample_correlations(
            combined_df,
            output_directory,
            data_file_type,
            column_identifier,
            data_columns="(.*)"+column_identifier+"(.*)", 
            plot_type="heatmap",
            font_size=12
        )
    else:
        print("No data available for correlation plot.")

    create_interactive_rank_plot(all_data_frames, labels, column_identifier, output_directory, data_file_type)
    create_static_rank_plot(all_data_frames, labels, column_identifier, output_directory, data_file_type)
    create_venn_or_upset_plot(intensity_data_frames, column_identifier, output_directory)
    # create_scatter_density_plot(all_data_frames, labels, column_identifier, output_directory, data_file_type, x_value="index", y_value="CV")
    # create_scatter_density_plot(all_data_frames, labels, column_identifier, output_directory, data_file_type, x_value="index", y_value="missing_value_ratio")
    create_scatter_density_plot(all_data_frames, labels, column_identifier, output_directory, data_file_type, x_value="log2_median", y_value="CV")
    # create_scatter_density_plot(all_data_frames, labels, column_identifier, output_directory, data_file_type, x_value="log2_median", y_value="missing_value_ratio")

def merge_intensity_dataframes(intensity_data_frames: List[pd.DataFrame]) -> pd.DataFrame:
    """
    Merge multiple intensity dataframes based on the 'unique_id' column.

    This function takes a list of pandas DataFrames, each containing intensity data
    and a 'unique_id' column. It merges these dataframes on the 'unique_id' column,
    keeping all unique IDs even if they don't appear in all dataframes.

    Args:
        intensity_data_frames (List[pd.DataFrame]): A list of pandas DataFrames
            containing intensity data. Each DataFrame must have a 'unique_id' column.

    Returns:
        pd.DataFrame: A merged DataFrame containing all intensity columns from all
            input dataframes, with rows aligned based on the 'unique_id'.

    Note:
        If a unique_id is not present in one of the dataframes, the corresponding
        intensity values will be NaN in the merged DataFrame.
    """
    if not intensity_data_frames:
        return pd.DataFrame()

    # Convert 'unique_id' to string in all dataframes
    for df in intensity_data_frames:
        df['unique_id'] = df['unique_id'].astype(str)

    # Start with the first dataframe
    merged_df: pd.DataFrame = intensity_data_frames[0]

    # Merge with the rest of the dataframes
    for df in intensity_data_frames[1:]:
        merged_df = pd.merge(merged_df, df, on='unique_id', how='outer')

    return merged_df

def _get_filtered_indices(column_identifier_index: int, num_identifiers: int, num_files: int, filter_columns: List[str]) -> List[int]:
    """
    Get filtered indices based on the given parameters.

    Args:
        column_identifier_index (int): Index of the column identifier.
        num_identifiers (int): Number of identifiers.
        num_files (int): Number of files.
        filter_columns (List[str]): List of columns to filter by.

    Returns:
        List[int]: List of filtered indices.
    """
    if filter_columns:
        indices = [column_identifier_index + i 
                   for i in range(len(filter_columns) * num_files)]
        return indices
    indices = [column_identifier_index + num_identifiers * file_num 
               for file_num in range(num_files)]
    return indices

def get_viridis_colors(n_colors: int) -> List[Tuple[float, float, float, float]]:
    """
    Get a list of colors from the Viridis colormap.

    Args:
        n_colors (int): Number of colors to generate.

    Returns:
        List[Tuple[float, float, float, float]]: List of RGBA color tuples.
    """
    viridis = cm.get_cmap('viridis')
    return [viridis(i / (n_colors - 1)) for i in range(n_colors)]

def create_stacked_bar_plot(
    df: pd.DataFrame,
    names: List[str],
    save_name: str,
    ylabel: str,
    xlabel: str,
    title: str,
    label_box: List[str]
) -> None:
    """
    Create a stacked bar plot.

    Args:
        df (pd.DataFrame): DataFrame containing the data to plot.
        names (List[str]): List of names for x-axis categories.
        save_name (str): Name to use when saving the plot.
        ylabel (str): Label for y-axis.
        xlabel (str): Label for x-axis.
        title (str): Title of the plot.
        label_box (List[str]): List of labels for each segment in the stacked bar.

    Returns:
        None
    """
    n_colors = len(label_box)
    plot_colors: List[Tuple[float, float, float, float]] = get_viridis_colors(n_colors)
    
    fig, ax = plt.subplots(figsize=(12, 6))

    bottom: np.ndarray = np.zeros(len(df))

    for i, column in enumerate(df.columns[::-1]):
        ax.bar(names, df[column], bottom=bottom, label=label_box[-(i+1)], 
               color=plot_colors[-(i+1)], edgecolor='white', width=0.5)
        bottom += df[column]

    ax.set_xlabel(xlabel, fontweight='bold')
    ax.set_ylabel(ylabel, fontweight='bold')
    ax.set_title(title, fontweight='bold')
    ax.set_xticklabels(names, rotation=0, ha='center')
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1), ncol=1)

    # Rotate x-axis labels
    ax.set_xticklabels(names, rotation=90, ha='center', va='top')

    plt.tight_layout()
    save_plot(save_name)

def save_plot(save_name: str) -> None:
    """
    Save the current plot as both PDF and PNG files.

    Args:
        save_name (str): Base name to use when saving the plot.

    Returns:
        None
    """
    plt.savefig(f"{save_name}.pdf", dpi=300)
    plt.savefig(f"{save_name}.png", dpi=300)
    plt.show()
    plt.close()

def get_y_label(data_file_type: str) -> str:
    """
    Get the appropriate y-axis label based on the data file type.

    Args:
        data_file_type (str): Type of data file being analyzed.

    Returns:
        str: Appropriate y-axis label.
    """
    y_labels: Dict[str, str] = {
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
    return y_labels.get(data_file_type, "Count")

def create_cv_bar_plot(
    filtered_results_df: pd.DataFrame,
    labels: List[str],
    output_directory: str,
    data_file_type: str,
    column_identifier: str
) -> None:
    """
    Create a CV bar plot.

    Args:
        filtered_results_df (pd.DataFrame): DataFrame containing filtered results.
        labels (List[str]): List of labels for the data.
        output_directory (str): Directory to save output files.
        data_file_type (str): Type of data file being analyzed.
        column_identifier (str): Identifier for the current column.

    Returns:
        None
    """
    save_name: str = f"{output_directory}/bar_plot_CVs_{data_file_type}_{column_identifier}"
    ylabel: str = get_y_label(data_file_type)
    xlabel: str = 'Type'
    title: str = column_identifier
    label_box: List[str] = ['Total No.', '<= CV 20%', '<= CV 10%', '<= CV 5%']
    
    required_columns: List[str] = ['totalIDs', 'CV20Percent', 'CV10Percent', 'CV5Percent']
    if not all(col in filtered_results_df.columns for col in required_columns):
        print(f"Error: Missing required columns. Available columns: {filtered_results_df.columns}")
        return

    print(filtered_results_df[required_columns])
    plot_data: pd.DataFrame = pd.DataFrame({
        'Total No.': filtered_results_df['totalIDs'] - filtered_results_df['CV20Percent'],
        '<= CV 20%': filtered_results_df['CV20Percent'] - filtered_results_df['CV10Percent'],
        '<= CV 10%': filtered_results_df['CV10Percent'] - filtered_results_df['CV5Percent'],
        '<= CV 5%': filtered_results_df['CV5Percent']
    })

    create_stacked_bar_plot(
        df=plot_data,
        names=labels,
        save_name=save_name,
        ylabel=ylabel,
        xlabel=xlabel,
        title=title,
        label_box=label_box
    )

def plot_identifications_per_run(
    filtered_results_df: pd.DataFrame,
    labels: List[str],
    output_directory: str,
    data_file_type: str,
    column_identifier: str,
    dot_size: int = 15,
    cap_size: int = 25
) -> None:
    """
    Plot identifications per run.

    Args:
        filtered_results_df (pd.DataFrame): DataFrame containing filtered results.
        labels (List[str]): List of labels for the data.
        output_directory (str): Directory to save output files.
        data_file_type (str): Type of data file being analyzed.
        column_identifier (str): Identifier for the current column.
        dot_size (int, optional): Size of dots in the swarm plot. Defaults to 15.
        cap_size (int, optional): Size of error bar caps. Defaults to 25.

    Returns:
        None
    """
    mean_list: List[float] = []
    error_list: List[float] = []
    jitter_data: List[Tuple[float, str]] = []

    for i, label in enumerate(labels):
        if 'IDs_per_run' in filtered_results_df.columns:
            ids_per_run = filtered_results_df['IDs_per_run'].iloc[i]
            if isinstance(ids_per_run, list):
                mean_list.append(np.mean(ids_per_run))
                error_list.append(np.std(ids_per_run))
                jitter_data.extend([(id_value, label) for id_value in ids_per_run])
            else:
                print(f"Warning: IDs_per_run for {label} is not a list. Skipping this entry.")
        else:
            print(f"Warning: 'IDs_per_run' column not found in DataFrame. Available columns: {filtered_results_df.columns}")
            return

    jitter_df: pd.DataFrame = pd.DataFrame(jitter_data, columns=['value', 'Type'])

    print("Mean IDs:", mean_list)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(labels, mean_list, yerr=error_list, align='center', ecolor='black', capsize=cap_size, color='#3777ac', label="average")
    sns.swarmplot(x='Type', y='value', data=jitter_df, color="grey", alpha=0.8, size=dot_size, ax=ax)

    ax.set_ylabel(get_y_label(data_file_type))
    ax.set_title(column_identifier)
    plt.legend(loc="lower right", bbox_to_anchor=(1.4, 0.85), ncol=1)
    plt.xticks(rotation=90, ha='right')

    plt.tight_layout()
    save_name = f"{output_directory}/bar_plot_IDs_per_run_{data_file_type}_{column_identifier}"
    save_plot(save_name)

def generate_violin_plots(
    filtered_results_df: pd.DataFrame,
    labels: List[str],
    output_directory: str,
    data_file_type: str,
    column_identifier: str
) -> None:
    """
    Generate violin plots for CV values.

    Args:
        filtered_results_df (pd.DataFrame): DataFrame containing filtered results.
        labels (List[str]): List of labels for the data.
        output_directory (str): Directory to save output files.
        data_file_type (str): Type of data file being analyzed.
        column_identifier (str): Identifier for the current column.

    Returns:
        None
    """
    cv_data: List[Tuple[float, str]] = []

    if 'CV' not in filtered_results_df.columns:
        print(f"Warning: 'CV' column not found in DataFrame. Available columns: {filtered_results_df.columns}")
        return

    for i, label in enumerate(labels):
        cv_values = filtered_results_df['CV'].iloc[i]
        if isinstance(cv_values, list):
            cv_data.extend([(cv, label) for cv in cv_values])
        else:
            print(f"Warning: CV values for {label} is not a list. Skipping this entry.")

    if not cv_data:
        print("No valid CV data found. Unable to generate violin plot.")
        return

    cv_df: pd.DataFrame = pd.DataFrame(cv_data, columns=['CV', 'Type'])

    fig, ax = plt.subplots(figsize=(12, 6))
    sns.violinplot(data=cv_df, x='Type', y='CV', palette="Blues", ax=ax)

    ax.set_ylabel(f"CV ({get_y_label(data_file_type)[2:]})")
    ax.set_title(column_identifier)
    plt.xticks(rotation=90, ha='right')

    plt.tight_layout()
    save_name = f"{output_directory}/violin_plot_CVs_{data_file_type}_{column_identifier}"
    save_plot(save_name)

def is_consistent_pattern(levels: List[Tuple[float, ...]]) -> bool:
    """
    Check if the data completeness levels follow a consistent pattern.
    
    Args:
        levels (List[Tuple[float, ...]]): List of data completeness level tuples.
    
    Returns:
        bool: True if the pattern is consistent, False otherwise.
    """
    # Check if all sequences end with 1.0
    if not all(seq[-1] == 1.0 for seq in levels):
        return False
    
    # Check if all sequences have a consistent step size
    step_sizes = set()
    for seq in levels:
        if len(seq) > 1:
            step = round(seq[1] - seq[0], 4)  # Round to 4 decimal places to avoid floating point issues
            step_sizes.add(step)
            if not all(round(seq[i+1] - seq[i], 4) == step for i in range(len(seq)-1)):
                return False
    
    # Check if all sequences use the same step size
    return len(step_sizes) == 1

def plot_missing_value_bars(
    filtered_results_df: pd.DataFrame,
    labels: List[str],
    output_directory: str,
    data_file_type: str,
    column_identifier: str
) -> None:
    """
    Plot stacked bar chart for missing value ratios if all datasets have an equal number of entries.

    Args:
        filtered_results_df (pd.DataFrame): DataFrame containing filtered results.
        labels (List[str]): List of labels for each dataset.
        output_directory (str): Directory to save output files.
        data_file_type (str): Type of data file being analyzed.
        column_identifier (str): Identifier for the current column.

    Returns:
        None
    """
    print("Checking data completeness levels...")
    
    # Extract data completeness levels
    data_completeness_levels = [tuple(data_completeness[0]) for data_completeness in filtered_results_df["DataCompleteness"]]
    print("Data completeness levels:", data_completeness_levels)
    
    # Check if the pattern is consistent
    consistent = is_consistent_pattern(data_completeness_levels)
    
    if not consistent:
        print("Error: Inconsistent data completeness levels across datasets.")
        for label, levels in zip(labels, data_completeness_levels):
            print(f"{label}: {levels}")
        print("Skipping the missing value bar plot.")
        return

    print("Creating stacked bar plot for missing values")
    save_name: str = f"{output_directory}/bar_plot_DataCompleteness_{data_file_type}_{column_identifier}"
    ylabel: str = get_y_label(data_file_type)
    xlabel: str = 'Type'
    title: str = f'Data Completeness - {column_identifier}'
    label_box: List[str] = ['Total No.', '<= CV 20%', '<= CV 10%', '<= CV 5%']

    temp_data_frames: List[pd.DataFrame] = []

    for label, data_completeness in zip(labels, filtered_results_df["DataCompleteness"]):
        df_temp = pd.DataFrame({
            "DataCompleteness": data_completeness[0],
            "Count": data_completeness[1],
            "Type": label
        })
        temp_data_frames.append(df_temp)

    df = pd.concat(temp_data_frames).reset_index(drop=True)
    
    # Pivot the dataframe to get it in the right format for stacking
    df_pivot = df.pivot(index='Type', columns='DataCompleteness', values='Count')
    print(df_pivot)
    
    # Sort columns in descending order (1.0 to 0.25)
    df_pivot = df_pivot.sort_index(axis=1, ascending=False)
    
    # Calculate the differences for stacking
    df_stacked = df_pivot.copy()
    for col in df_pivot.columns[1:]:
        df_stacked[col] = df_pivot[col] - df_pivot[df_pivot.columns[df_pivot.columns.get_loc(col) - 1]]

    create_stacked_bar_plot(
        df=df_stacked.sort_index(axis=1, ascending=True),
        names=labels,
        save_name=save_name,
        ylabel=ylabel,
        xlabel=xlabel,
        title=title,
        label_box=df_stacked.columns
    )

def create_interactive_rank_plot(
    all_data_frames: List[pd.DataFrame],
    labels: List[str],
    column_identifier: str,
    output_directory: str,
    data_file_type: str,
    x_value: str = "index",
    y_value: str = "log10_median"
) -> px.scatter:
    """
    Create an interactive rank plot.

    Args:
        all_data_frames (List[pd.DataFrame]): List of all DataFrames.
        labels (List[str]): List of labels for the data.
        column_identifier (str): Identifier for the current column.
        output_directory (str): Directory to save output files.
        data_file_type (str): Type of data file being analyzed.
        x_value (str, optional): Column name for x-axis values. Defaults to "index".
        y_value (str, optional): Column name for y-axis values. Defaults to "log10_median".

    Returns:
        px.scatter: The created interactive scatter plot.
    """
    plot_colors: List[str] = px.colors.sample_colorscale('blues', len(labels) + 1)
    fig: px.scatter = px.scatter()

    for num, (intensity_df, label) in enumerate(zip(all_data_frames, labels)):
        if column_identifier + "_" in intensity_df.columns.str.cat():
            processed_df: pd.DataFrame = prepare_dataframe_for_plotting(intensity_df, column_identifier, x_value, y_value)
            fig.add_scatter(x=processed_df[x_value], y=processed_df[y_value], 
                            mode="markers", marker=dict(color=plot_colors[num+1]), name=label)

            if x_value == "index" and y_value == "log10_median":
                print_quartile_info(processed_df, label)
                print_dynamic_range_info(processed_df)

    update_layout(fig, x_value, y_value, column_identifier)
    fig.show()
    fig.write_image(f"{output_directory}/plot_{get_axis_label(x_value, column_identifier)}_{get_axis_label(y_value, column_identifier)}_{data_file_type}_{column_identifier}.pdf")
    return fig

def print_quartile_info(df: pd.DataFrame, label: str) -> None:
    """
    Print quartile information for the given DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing the data.
        label (str): Label for the current dataset.

    Returns:
        None
    """
    print(f"{label}:", "quartile, position on the y-axis:")
    df_quartile: pd.Series = df["log10_median"].quantile([0.25, 0.5, 0.75])
    for quartile in [0.25, 0.5, 0.75]:
        print(quartile, df_quartile[quartile])

def print_dynamic_range_info(df: pd.DataFrame) -> None:
    """
    Print dynamic range information for the given DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing the data.

    Returns:
        None
    """
    diff_max_min: float = (df["log10_median"].max() - df["log10_median"].min()) / 4
    percentile_ranges: List[List[float]] = [
        [df["log10_median"].min() + i * diff_max_min, df["log10_median"].min() + (i + 1) * diff_max_min]
        for i in range(4)
    ]

    print("number of precursors per 25% of the dynamic range")
    for percentile in percentile_ranges:
        precursors_within_range: int = len(df[(df["log10_median"] >= percentile[0]) & (df["log10_median"] < percentile[1])])
        print(f"{percentile} (25%):")
        print(f"count: {precursors_within_range}, ratio from total IDs: {precursors_within_range / len(df['log10_median']):.2f}")

def update_layout(fig: px.scatter, x_value: str, y_value: str, column_identifier: str) -> None:
    """
    Update the layout of the given scatter plot.

    Args:
        fig (px.scatter): The scatter plot to update.
        x_value (str): Column name for x-axis values.
        y_value (str): Column name for y-axis values.
        column_identifier (str): Identifier for the current column.

    Returns:
        None
    """
    x_axis_label: str = get_axis_label(x_value, column_identifier)
    y_axis_label: str = get_axis_label(y_value, column_identifier)
    fig.update_layout(
        title=dict(text=f"{x_axis_label} versus {y_axis_label}", font=dict(size=16), x=0.5, xanchor='center', yanchor='top', y=0.92),
        xaxis=dict(title=x_axis_label, titlefont_size=14, tickfont_size=14),
        yaxis=dict(title=y_axis_label),
        legend_title_text='Legend:',
        template='plotly_white',
        height=450
    )

def create_static_rank_plot(
    all_data_frames: List[pd.DataFrame],
    labels: List[str],
    column_identifier: str,
    output_directory: str,
    data_file_type: str,
    x_value: str = "index",
    y_value: str = "log10_median"
) -> Optional[List[pd.DataFrame]]:
    """
    Create a static rank plot.

    Args:
        all_data_frames (List[pd.DataFrame]): List of all DataFrames.
        labels (List[str]): List of labels for the data.
        column_identifier (str): Identifier for the current column.
        output_directory (str): Directory to save output files.
        data_file_type (str): Type of data file being analyzed.
        x_value (str, optional): Column name for x-axis values. Defaults to "index".
        y_value (str, optional): Column name for y-axis values. Defaults to "log10_median".

    Returns:
        Optional[List[pd.DataFrame]]: List of DataFrames for scatter plots if y_value is not "log10", otherwise None.
    """
    plot_colors: np.ndarray = plt.cm.viridis(np.linspace(0, 1, len(labels) + 1))
    scatter_plot_data: List[pd.DataFrame] = []

    for num, (intensity_df, label) in enumerate(zip(all_data_frames, labels)):
        if column_identifier + "_" in intensity_df.columns.str.cat():
            processed_df: pd.DataFrame = prepare_dataframe_for_plotting(intensity_df, column_identifier, x_value, y_value)
            
            if "log10" in y_value:
                plt.plot(processed_df[x_value], processed_df[y_value], color=plot_colors[num + 1], linewidth=2, label=label)
            else:
                temp_df: pd.DataFrame = pd.DataFrame({
                    x_value: processed_df[x_value],
                    y_value: processed_df[y_value],
                    "Type: ": label
                })
                scatter_plot_data.append(temp_df)

    x_axis_label: str = get_axis_label(x_value, column_identifier)
    y_axis_label: str = get_axis_label(y_value, column_identifier)

    if "log10" in y_value:
        plt.xlabel(x_axis_label)
        plt.ylabel(y_axis_label)
        plt.title(f"{x_axis_label} versus {y_axis_label}")
        plt.legend()
        save_name = f"{output_directory}/plot_{x_axis_label}_versus_{y_axis_label}_{data_file_type}_{column_identifier}"
        save_plot(save_name)
        plt.clf()
        return None
    else:
        return scatter_plot_data

def prepare_dataframe_for_plotting(df: pd.DataFrame, column_identifier: str, x_value: str, y_value: str) -> pd.DataFrame:
    """
    Prepare a DataFrame for plotting by calculating various metrics.

    Args:
        df (pd.DataFrame): Input DataFrame.
        column_identifier (str): Identifier for the current column.
        x_value (str): Column name for x-axis values.
        y_value (str): Column name for y-axis values.

    Returns:
        pd.DataFrame: Processed DataFrame ready for plotting.
    """
    selected_columns: List[str] = [col for col in df.columns if column_identifier in col]
    df['median'] = df[selected_columns].median(axis=1)
    df["log10_median"] = np.log10(df["median"])
    df["log2_median"] = np.log2(df["median"])
    df["CV"] = df[f"CV_{column_identifier}"]
    df["missing_value_ratio"] = df[f"data_completeness_{column_identifier}"]
    sorted_df: pd.DataFrame = df.sort_values(by="log10_median", ascending=False).reset_index(drop=True)
    sorted_df["index"] = sorted_df.index
    return sorted_df

def create_scatter_density_plot(
    all_data_frames: List[pd.DataFrame],
    labels: List[str],
    column_identifier: str,
    output_directory: str,
    data_file_type: str,
    x_value: str = "index",
    y_value: str = "log10_median",
    use_kernel_density_estimation: bool = False
) -> None:
    """
    Create a scatter density plot.

    Args:
        all_data_frames (List[pd.DataFrame]): List of all DataFrames.
        labels (List[str]): List of labels for the data.
        column_identifier (str): Identifier for the current column.
        output_directory (str): Directory to save output files.
        data_file_type (str): Type of data file being analyzed.
        x_value (str, optional): Column name for x-axis values. Defaults to "index".
        y_value (str, optional): Column name for y-axis values. Defaults to "log10_median".
        use_kernel_density_estimation (bool, optional): Whether to use kernel density estimation. Defaults to False.

    Returns:
        None
    """
    x_axis_label: str = get_axis_label(x_value, column_identifier)
    y_axis_label: str = get_axis_label(y_value, column_identifier)

    scatter_plot_data: List[pd.DataFrame] = []
    for df, label in zip(all_data_frames, labels):
        if column_identifier + "_" in df.columns.str.cat():
            processed_df: pd.DataFrame = prepare_dataframe_for_plotting(df, column_identifier, x_value, y_value)
            temp_df: pd.DataFrame = pd.DataFrame({
                x_value: processed_df[x_value],
                y_value: processed_df[y_value],
                "Type": label
            })
            scatter_plot_data.append(temp_df)

    combined_scatter_data: pd.DataFrame = pd.concat(scatter_plot_data).reset_index(drop=True)
    
    print("Please wait - it may take several minutes")
    
    joint_plot: sns.JointGrid = create_joint_plot(combined_scatter_data, x_value, y_value, use_kernel_density_estimation)
    
    plt.xlabel(x_axis_label)
    plt.ylabel(y_axis_label)
    save_name = f"{output_directory}/plot_{x_axis_label}_versus_{y_axis_label}_{data_file_type}_{column_identifier}_{'w' if use_kernel_density_estimation else 'wo'}_kde"
    joint_plot.savefig(f"{save_name}.pdf", bbox_inches='tight', pad_inches=0, dpi=300)
    plt.show()

def create_joint_plot(df: pd.DataFrame, x_value: str, y_value: str, use_kernel_density_estimation: bool) -> sns.JointGrid:
    """
    Create a joint plot using seaborn.

    Args:
        df (pd.DataFrame): DataFrame containing the data to plot.
        x_value (str): Column name for x-axis values.
        y_value (str): Column name for y-axis values.
        use_kernel_density_estimation (bool): Whether to use kernel density estimation.

    Returns:
        sns.JointGrid: The created joint plot.
    """
    if use_kernel_density_estimation:
        return sns.jointplot(data=df, x=x_value, y=y_value, hue="Type", kind="kde")
    else:
        return sns.jointplot(data=df, x=x_value, y=y_value, hue="Type", s=5, edgecolor="face")

def get_axis_label(value: str, column_identifier: str) -> str:
    """
    Get the appropriate axis label based on the value and column identifier.

    Args:
        value (str): The value to get the label for.
        column_identifier (str): Identifier for the current column.

    Returns:
        str: The appropriate axis label.
    """
    axis_labels: Dict[str, str] = {
        "index": "Rank", 
        "log10_median": "Median Intensities (log10)", 
        "log2_median": "Median Intensities (log2)",
        f"CV_{column_identifier}": "CV", 
        f"data_completeness_{column_identifier}": "Missing Value Ratio"
    }
    return axis_labels.get(value, value)

def create_venn_or_upset_plot(
    intensity_data_frames: List[pd.DataFrame],
    column_identifier: str,
    output_directory: str
) -> None:
    """
    Create either a Venn diagram or an UpSet plot based on the number of datasets.

    Args:
        intensity_data_frames (List[pd.DataFrame]): List of DataFrames with intensity data.
        column_identifier (str): Identifier for the current column.
        output_directory (str): Directory to save output files.

    Returns:
        None
    """
    all_identifications: Set[str]
    set_dictionary: Dict[str, Set[str]]
    all_identifications, set_dictionary = process_dataframes_for_set_analysis(intensity_data_frames, column_identifier)

    if len(set_dictionary) == 2:
        print("Creating Venn diagram")
        create_venn2_diagram(set_dictionary, output_directory, column_identifier)
    elif len(set_dictionary) == 3:
        print("Creating Venn diagram")
        create_venn3_diagram(set_dictionary, output_directory, column_identifier)
    else:
        print("More than three files: Creating UpSet plot")
        create_upset_plot(set_dictionary, all_identifications, output_directory, column_identifier)

def process_dataframes_for_set_analysis(
    data_frames: List[pd.DataFrame],
    column_identifier: str
) -> Tuple[Set[str], Dict[str, Set[str]]]:
    """
    Process DataFrames to create sets of unique IDs for set analysis.

    Args:
        data_frames (List[pd.DataFrame]): List of DataFrames to process.
        column_identifier (str): Identifier for the current column.

    Returns:
        Tuple[Set[str], Dict[str, Set[str]]]: A tuple containing a set of all unique IDs and a dictionary mapping labels to sets of unique IDs.
    """
    all_identifications: Set[str] = set()
    set_dictionary: Dict[str, Set[str]] = {}

    for df in data_frames:
        if column_identifier in df.columns.str.cat():
            label: str = df.columns[0].split(column_identifier)[0]
            identifications: Set[str] = set(df["unique_id"])
            all_identifications = all_identifications.union(identifications)
            set_dictionary[label] = identifications

    return all_identifications, set_dictionary

def create_venn2_diagram(
    set_dictionary: Dict[str, Set[str]],
    output_directory: str,
    column_identifier: str
) -> None:
    """
    Create a 2-set Venn diagram.

    Args:
        set_dictionary (Dict[str, Set[str]]): Dictionary mapping labels to sets of unique IDs.
        output_directory (str): Directory to save output files.
        column_identifier (str): Identifier for the current column.

    Returns:
        None
    """
    plot_colors: np.ndarray = cm.get_cmap('viridis', 2).colors

    venn2(
        [value for value in set_dictionary.values()],
        set_labels = list(set_dictionary.keys()),
        set_colors=plot_colors, alpha=0.7
    )
    plt.title(column_identifier)
    save_venn_diagram(output_directory, column_identifier)

def create_venn3_diagram(
    set_dictionary: Dict[str, Set[str]],
    output_directory: str,
    column_identifier: str
) -> None:
    """
    Create a 3-set Venn diagram.

    Args:
        set_dictionary (Dict[str, Set[str]]): Dictionary mapping labels to sets of unique IDs.
        output_directory (str): Directory to save output files.
        column_identifier (str): Identifier for the current column.

    Returns:
        None
    """
    plot_colors: np.ndarray = cm.get_cmap('viridis', 3).colors

    venn3(
        [value for value in set_dictionary.values()],
        set_labels = list(set_dictionary.keys()),
        set_colors=plot_colors, alpha=0.7
    )
    plt.title(column_identifier)
    save_venn_diagram(output_directory, column_identifier)

def save_venn_diagram(output_directory: str, column_identifier: str) -> None:
    """
    Save the Venn diagram to a file.

    Args:
        output_directory (str): Directory to save output files.
        column_identifier (str): Identifier for the current column.

    Returns:
        None
    """
    save_name = f"{output_directory}/venn_diagram_{column_identifier}"
    save_plot(save_name)

def create_upset_plot(
    set_dictionary: Dict[str, Set[str]],
    all_identifications: Set[str],
    output_directory: str,
    column_identifier: str
) -> None:
    """
    Create an UpSet plot.

    Args:
        set_dictionary (Dict[str, Set[str]]): Dictionary mapping labels to sets of unique IDs.
        all_identifications (Set[str]): Set of all unique IDs.
        output_directory (str): Directory to save output files.
        column_identifier (str): Identifier for the current column.

    Returns:
        None
    """
    df: pd.DataFrame = pd.DataFrame({key: [e in value for e in all_identifications] for key, value in set_dictionary.items()})
    upset_data: pd.Series = df.groupby(list(set_dictionary.keys())).size()
    
    upset_plot(upset_data, orientation='horizontal')
    plt.suptitle(column_identifier)
    save_name = f"{output_directory}/upset_plot_{column_identifier}"
    save_plot(save_name)