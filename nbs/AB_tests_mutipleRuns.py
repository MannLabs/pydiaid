# Author: Patricia Skowronek, Max Planck Institute of Biochemistry

import pandas as pd
import numpy as np

from matplotlib import rc
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16.5})
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

from difflib import SequenceMatcher

import re
from plotly.subplots import make_subplots
from upsetplot import plot as upset_plot #Alexander Lex, Nils Gehlenborg, Hendrik Strobelt, Romain Vuillemot, Hanspeter Pfister, UpSet: Visualization of Intersecting Sets, IEEE Transactions on Visualization and Computer Graphics (InfoVis ‘14), vol. 20, no. 12, pp. 1983–1992, 2014. doi: doi.org/10.1109/TVCG.2014.2346248
from matplotlib_venn import venn2, venn3
from matplotlib import cm
import datashader as ds
import plotly.graph_objects as go
import plotly.express as px
from plotly.express.colors import sample_colorscale

import warnings
warnings.filterwarnings('ignore')

def process_tables_depending_on_search_software(
    file_names,
    percentage_of_missing_values,
    save_output_here,
    labels,
    list_dfs_only_with_intensity_columns,
    list_dfs,
    dict_results,
):

    num = 0
    for file_name in file_names:
        
        df = pd.read_csv(file_name, sep='\t')  # .xls, .tsv, .txt
        
        # filter MaxQuant protein groups dataframe
        if "Only identified by site" in df.columns:
            file_type = "mq_proteins"
            unique_ID_column = "Protein IDs"
            df_filtered = df[(df["Only identified by site"]!="+") & (df["Reverse"]!="+") & (df["Potential contaminant"]!="+")]
            column_indicators = ["Intensity", "iBAQ",  "LFQ intensity"]
            
            for column_indicator in column_indicators:
                df_filtered, column_indicator, selected_columns = filter_MaxQuant_output_for_column_indicator(df_filtered, column_indicator)
                
                dict_results, list_dfs_only_with_intensity_columns, list_dfs = filter_tables_and_make_correlation_plots(
                    df_filtered, 
                    dict_results,
                    column_indicator, 
                    percentage_of_missing_values,
                    selected_columns,
                    unique_ID_column,
                    save_output_here,
                    file_type,
                    labels[num],
                    list_dfs_only_with_intensity_columns,
                    list_dfs,
                    file_name
                ) 
            num += 1
            
        # filter MaxQuant peptides dataframe
        elif "A Count" in df.columns:
            file_type = "mq_peptides"
            unique_ID_column = "Sequence"
            df_filtered = df[(df["Reverse"]!="+") & (df["Potential contaminant"]!="+")]
            column_indicators = ["Intensity", "LFQ intensity"]
            
            for column_indicator in column_indicators:
                df_filtered, column_indicator, selected_columns = filter_MaxQuant_output_for_column_indicator(df_filtered, column_indicator)
                
                dict_results, list_dfs_only_with_intensity_columns, list_dfs = filter_tables_and_make_correlation_plots(
                    df_filtered, 
                    dict_results,
                    column_indicator, 
                    percentage_of_missing_values,
                    selected_columns,
                    unique_ID_column,
                    save_output_here,
                    file_type,
                    labels[num],
                    list_dfs_only_with_intensity_columns,
                    list_dfs,
                    file_name
                )
            num += 1
            
        # filter AlphaDIA protein groups dataframe
        if ("pg" in df.columns) & ("base_width_mobility" not in df.columns):
            file_type = "alphadia_proteins"
            unique_ID_column = "pg"
            df_filtered = df.copy()
            column_indicators = ["Intensity"]

            for column_indicator in column_indicators:
                df_filtered, column_indicator, selected_columns = filter_alphadia_output_for_column_indicator(df)

                dict_results, list_dfs_only_with_intensity_columns, list_dfs = filter_tables_and_make_correlation_plots(
                    df_filtered, 
                    dict_results,
                    column_indicator, 
                    percentage_of_missing_values,
                    selected_columns,
                    unique_ID_column,
                    save_output_here,
                    file_type,
                    labels[num],
                    list_dfs_only_with_intensity_columns,
                    list_dfs,
                    file_name
                )
            num += 1
            
        # filter AlphaDIA precursor dataframe
        if "mod_seq_charge_hash" in df.columns:
            file_type = "alphadia_precursors"
            unique_ID_column = "mod_seq_charge_hash"
            df_filtered = df.copy()
            column_indicators = ["Intensity"]

            for column_indicator in column_indicators:
                df_filtered, column_indicator, selected_columns = filter_alphadia_output_for_column_indicator(df)

                dict_results, list_dfs_only_with_intensity_columns, list_dfs = filter_tables_and_make_correlation_plots(
                    df_filtered, 
                    dict_results,
                    column_indicator, 
                    percentage_of_missing_values,
                    selected_columns,
                    unique_ID_column,
                    save_output_here,
                    file_type,
                    labels[num],
                    list_dfs_only_with_intensity_columns,
                    list_dfs,
                    file_name
                )
            num += 1

        # filter DIA-NN protein dataframe
        elif ("First.Protein.Description" in df.columns) & ("Stripped.Sequence" not in df.columns):
            file_type = "diann_proteins"
            unique_ID_column = "Protein.Group"
            df_filtered = df.copy()
            column_indicators = ["Intensity"]
            
            for column_indicator in column_indicators:
                df_filtered, column_indicator, selected_columns = filter_DIANN_output_for_column_indicator(df)
                
                dict_results, list_dfs_only_with_intensity_columns, list_dfs = filter_tables_and_make_correlation_plots(
                    df_filtered, 
                    dict_results,
                    column_indicator, 
                    percentage_of_missing_values,
                    selected_columns,
                    unique_ID_column,
                    save_output_here,
                    file_type,
                    labels[num],
                    list_dfs_only_with_intensity_columns,
                    list_dfs,
                    file_name
                )
            num += 1
            
        # filter DIA-NN precursors dataframe
        elif ("First.Protein.Description" in df.columns) & ("Stripped.Sequence" in df.columns):
            file_type = "diann_precursors"
            unique_ID_column = "Precursor.Id"
            df_filtered = df.copy()
            column_indicators = ["Intensity"]
            
            for column_indicator in column_indicators:
                df_filtered, column_indicator, selected_columns = filter_DIANN_output_for_column_indicator(df)
                
                dict_results, list_dfs_only_with_intensity_columns, list_dfs = filter_tables_and_make_correlation_plots(
                    df_filtered, 
                    dict_results,
                    column_indicator, 
                    percentage_of_missing_values,
                    selected_columns,
                    unique_ID_column,
                    save_output_here,
                    file_type,
                    labels[num],
                    list_dfs_only_with_intensity_columns,
                    list_dfs,
                    file_name
                )
            num += 1
            
        # filter spectronaut protein dataframe
        elif "PG.ProteinAccessions" in df.columns:
            file_type = "spectronaut_proteins"
            unique_ID_column = "PG.ProteinAccessions"
            df_filtered = df.copy()
            column_indicators = ["Intensity"]
            
            for column_indicator in column_indicators:
                df_filtered, column_indicator, selected_columns = filter_spectronaut_output_for_column_indicator(df)
                dict_results, list_dfs_only_with_intensity_columns, list_dfs = filter_tables_and_make_correlation_plots(
                    df_filtered, 
                    dict_results,
                    column_indicator, 
                    percentage_of_missing_values,
                    selected_columns,
                    unique_ID_column,
                    save_output_here,
                    file_type,
                    labels[num],
                    list_dfs_only_with_intensity_columns,
                    list_dfs,
                    file_name
                )
            num += 1
        # filter spectronaut peptides dataframe
        elif "EG.PrecursorId" in df.columns:
            file_type = "spectronaut_precursors"
            unique_ID_column = "EG.PrecursorId"
            df_filtered = df.copy()
            column_indicators = ["Intensity"]
            
            for column_indicator in column_indicators:
                df_filtered, column_indicator, selected_columns = filter_spectronaut_output_for_column_indicator(df)
                df_filtered_w_nan = df_filtered.replace('Filtered',np.NaN)
                df_filtered_floats = df_filtered_w_nan[selected_columns].astype(float)
                df_filtered_floats["EG.PrecursorId"] = df_filtered_w_nan["EG.PrecursorId"]
                
                dict_results, list_dfs_only_with_intensity_columns, list_dfs = filter_tables_and_make_correlation_plots(
                    df_filtered_floats, 
                    dict_results,
                    column_indicator, 
                    percentage_of_missing_values,
                    selected_columns,
                    unique_ID_column,
                    save_output_here,
                    file_type,
                    labels[num],
                    list_dfs_only_with_intensity_columns,
                    list_dfs,
                    file_name
                )
            num += 1
            
    else:
        next

    return column_indicators, dict_results, file_type, file_names, list_dfs_only_with_intensity_columns, list_dfs


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

        print(np.median(df_results_filtered_for_column_indicator["CV"].iloc[0]),
              np.median(df_results_filtered_for_column_indicator["CV"].iloc[1]),
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
        dfs_intensity_for_column_indicator["unique_ID_column"] = []
        for df_temp in list_dfs_only_with_intensity_columns:
            column_list = [col for col in df_temp.columns if column_indicator in col]
            column_list.append("unique_ID_column")
            df_selected = df_temp[column_list]
            dfs_intensity_for_column_indicator = dfs_intensity_for_column_indicator.merge(df_selected, left_on="unique_ID_column", right_on="unique_ID_column", how='outer')

        # function form Julia P. Schessner, Eugenia Voytik, Isabell Bludau,  https://doi.org/10.1002/pmic.202100103
        cross_corr2 = plot_sample_correlations(
            dfs_intensity_for_column_indicator,
            save_output_here,
            file_type,
            column_indicator,
            data_columns="(.*)"+column_indicator+"(.*)", 
            mode="heatmap",
            font_size=12
        )

        # abundance range illustrated as rank plot
        rank_plot(
            list_dfs, 
            labels,
            column_indicator,
            save_output_here,
            file_type,
        )

        # abundance range illustrated as rank plot
        rank_plot_static(
            list_dfs, 
            labels,
            column_indicator,
            save_output_here,
            file_type,
        )

    #     # CV values illustrated as rank plot
        list_scatter_plot_index = rank_plot_static(
            list_dfs, 
            labels,
            column_indicator,
            save_output_here,
            file_type,
            x_value="index",
            y_value= "CV_"+column_indicator,
        )

        scatter_density_plot(
            list_scatter_plot_index,
            labels,
            column_indicator,
            save_output_here,
            file_type,
            x_value="index",
            y_value= "CV_"+column_indicator,
    #         kde=True
        )

        # missing value illustrated as rank plot
        list_scatter_plot_missing_index = rank_plot_static(
            list_dfs, 
            labels,
            column_indicator,
            save_output_here,
            file_type,
            x_value="index",
            y_value= "missing_value_"+column_indicator,
        )
        scatter_density_plot(
            list_scatter_plot_missing_index,
            labels,
            column_indicator,
            save_output_here,
            file_type,
            x_value="index",
            y_value= "missing_value_"+column_indicator,
        )

        # CV illustrated along log2_median
        list_scatter_plot_log2 = rank_plot_static(
            list_dfs, 
            labels,
            column_indicator,
            save_output_here,
            file_type,
            x_value="log2_median",
            y_value= "CV_"+column_indicator,
        )
        scatter_density_plot(
            list_scatter_plot_log2,
            labels,
            column_indicator,
            save_output_here,
            file_type,
            x_value="log2_median",
            y_value= "CV_"+column_indicator,
    #         kde=True
        )

        # missing value illustrated along log2_median
        list_scatter_plot_missing_log2 = rank_plot_static(
            list_dfs, 
            labels,
            column_indicator,
            save_output_here,
            file_type,
            x_value="log2_median",
            y_value= "missing_value_"+column_indicator,
        )
        scatter_density_plot(
            list_scatter_plot_missing_log2,
            labels,
            column_indicator,
            save_output_here,
            file_type,
            x_value="log2_median",
            y_value= "missing_value_"+column_indicator,
        )

        # Venn diagram or upset plot to show shared IDs
        make_venn_or_upset_plot(
            list_dfs_only_with_intensity_columns,
            column_indicator,
            save_output_here
        )


def calc_cv(x):
    cv = np.nanstd(x, ddof=1) / np.nanmean(x)
    return cv

def filter_MaxQuant_output_for_column_indicator(df_filtered, column_indicator):
    selected_columns = [col for col in df_filtered.columns if column_indicator in col]

    if column_indicator == "iBAQ":
        selected_columns.remove('iBAQ peptides')
        selected_columns.remove('iBAQ')
    elif column_indicator == "Intensity":
        selected_columns.remove('Intensity')
    else:
        next
    for column in selected_columns:
        df_filtered[column][df_filtered[column] == 0] = np.nan

    return df_filtered, column_indicator, selected_columns

def filter_dataframe_for_results(df_prefiltered, column_indicator, selected_columns, dict_results):

    dict_results["totalIDs"].append(len(df_prefiltered))
    dict_results["CV20"].append(len(df_prefiltered[df_prefiltered["CV_"+column_indicator]<=0.2]))
    dict_results["CV10"].append(len(df_prefiltered[df_prefiltered["CV_"+column_indicator]<=0.1]))
    dict_results["CV5"].append(len(df_prefiltered[df_prefiltered["CV_"+column_indicator]<=0.05]))

    list_IDs_per_run = list()
    for column in selected_columns:
        list_IDs_per_run.append(len(df_prefiltered[column].dropna()))
    dict_results["IDs_per_run"].append(list_IDs_per_run)
    dict_results["CV"].append(list(df_prefiltered["CV_"+column_indicator]))

    missing_values = sorted(list(set(df_prefiltered["missing_value_"+column_indicator])))
    missing_values_total_list = list()
    for value in missing_values:
        missing_values_total_list.append(len(df_prefiltered[df_prefiltered["missing_value_"+column_indicator]>=value]))
    #dict_results["missing_values"].append(missing_values)
    dict_results["missing_values"].append([missing_values, missing_values_total_list])

    return dict_results

def filter_results_according_to_column_indicator(
    df_filtered,
    dict_results,
    column_indicator,
    percentage_of_missing_values,
    selected_columns,
    unique_ID_column
):

    df_filtered["CV_"+column_indicator] = df_filtered[selected_columns].apply(lambda x: calc_cv(x), axis = 1)
    df_filtered["missing_value_"+column_indicator] = df_filtered[selected_columns].apply(lambda x: len(x.dropna())/len(selected_columns), axis = 1)
    needed_columns = [unique_ID_column, "CV_"+column_indicator,"missing_value_"+column_indicator]+selected_columns
    df_prefiltered = df_filtered[needed_columns][df_filtered["missing_value_"+column_indicator]>=percentage_of_missing_values]
    df_prefiltered_2 = df_prefiltered[df_prefiltered["missing_value_"+column_indicator]>0]
    dict_results = filter_dataframe_for_results(df_prefiltered_2, column_indicator, selected_columns, dict_results)
    dict_results.update(dict_results)

    return dict_results, df_prefiltered_2

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
        df_results_filtered_for_column_indicator[['totalIDs', 'CV20', 'CV10', 'CV5']],
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
        list_filtered = df_results_filtered_for_column_indicator["missing_values"].iloc[x]
        df_temp = pd.DataFrame(np.array([list_filtered[1]]),
                               columns=list_filtered[0])
        df_temp["Type"] = labels[x]
        list_dfs_temp.append(df_temp)

    df = pd.concat(list_dfs_temp).reset_index(drop=True)

    save_name = save_output_here + "/bar_plot_missing_values_"+ file_type + "_" + column_indicator
    ylabel = find_y_label(file_type) # Name of the y axis for the bar plot
    xlabel = 'Type' # Name or topic of all the columns
    label_box = df.columns[:-1]*len(df_results_filtered_for_column_indicator["IDs_per_run"].iloc[0])

    bar_plot(
        df[df.columns[:-1]],
        labels,
        save_name,
        ylabel,
        xlabel,
        column_indicator,
        label_box
    )


def find_common_string(selected_columns):
    # https://stackoverflow.com/questions/58585052/find-most-common-substring-in-a-list-of-strings
    string2 = selected_columns[0]
    for i in range(1, len(selected_columns)):
        string1 = string2
        string2 = selected_columns[i]
        match = SequenceMatcher(None, string1, string2).find_longest_match(0, len(string1), 0, len(string2))

    return string1[match.a: match.a + match.size]


# This funciton is taken form https://github.com/MannLabs/ProteomicsVisualization: Notebook Figure4
# A practical guide to interpreting and generating bottom-up proteomics data visualizations,
# Julia Patricia Schessner, Eugenia Voytik, Isabell Bludau,  https://doi.org/10.1002/pmic.202100103

def plot_sample_correlations(
    df: pd.DataFrame,
    save_output_here,
    file_type,
    column_indicator,
    label="",
    data_columns: str = "Intensity (.*)",
    correlation_function: callable = lambda x: np.corrcoef(x.T),
    mode: str = "scatter", log10: bool = True, binning: int = 10,
    font_size=10,

    ):
    """
    Generates either a grid of paired scatter plots, or a heatmap of pairwise correlation values.

    Parameters
    ----------
    df: pandas DataFrame
    data_columns: regular expression string, default = "Intensity (.*)"
        Regular expression matching to columns containing data and a capture group for sample labels.
    correlation_function: callable, default = lambda x: np.corrcoef(x.T)
        Callable function to calculate correlation between samples.
    mode: str, default = "scatter"
        One of 'scatter' and 'heatmap' to swtich modes.
    log10: bool = True
        If True, data is log10 transformed prior to analysis
    binning: int = 10
        Only relevant in scatter mode. Number of bins per full axis unit (e.g. 1e10 Intensity)
        used for datashader density rendering.

    Returns
    -------
    a plotly Figure
        either a subplotted, or a single heatmap depending on the mode
    """
    # pick and process data
    df_sub = df[[el for el in df.columns if re.match(data_columns, el)]].copy()
    if log10:
        df_sub = df_sub.apply(np.log10)
    df_sub = df_sub.replace([np.inf, -np.inf], np.nan)
    df_sub.columns = [re.findall(data_columns, el)[0] for el in df_sub.columns]

    if mode == "scatter":
        # setup subplots and axes
        fig = make_subplots(rows=len(df_sub.columns), cols=len(df_sub.columns), start_cell='bottom-left',
                            shared_yaxes=True, shared_xaxes=True, horizontal_spacing=0.03, vertical_spacing=0.03)
        i_range = (np.floor(np.nanmin(df_sub)), np.ceil(np.nanmax(df_sub))+1/binning)
        j_range = (np.floor(np.nanmin(df_sub)), np.ceil(np.nanmax(df_sub))+1/binning)
        i_width = int((i_range[1]-i_range[0]-1/binning)*binning+1)
        j_width = int((j_range[1]-j_range[0]-1/binning)*binning+1)

        # fill plots
        for i,ni in enumerate(df_sub.columns):
            for j,nj in enumerate(df_sub.columns):
                # apply datashader
                dc = ds.Canvas(plot_width=i_width, plot_height=j_width, x_range=i_range, y_range=j_range)
                df_ij = df_sub[[ni,nj]].dropna() if i!=j else pd.DataFrame(df_sub[ni].dropna())
                da = dc.points(df_ij, x=ni, y=nj)
                zero_mask = da.values == 0
                da.values = da.values.astype(float)
                da.values[zero_mask] = np.nan

                # add trace
                fig.add_trace(
                    go.Heatmap(z=da,coloraxis="coloraxis1" if i!=j else "coloraxis2"),
                    row=j+1, col=i+1
                )

                # add annotations
                if j == 0:
                    fig.update_xaxes(title_text=ni, row=j+1, col=i+1, tickvals=list(range(0,i_width,binning)),
                                     ticktext=np.round(da[nj].values[0:i_width:binning]))
                if i == 0:
                    fig.update_yaxes(title_text=nj, row=j+1, col=i+1, tickvals=list(range(0,j_width,binning)),
                                     ticktext=np.round(da[ni].values[0:j_width:binning]))
                if i!=j:
                    fig.add_annotation(dict(text=str(np.round(np.min(correlation_function(df_sub[[ni,nj]].dropna())),4)),
                                            x=binning, y=j_width, showarrow=False), row=j+1, col=i+1)

        # layout figure
        fig.update_layout(template="simple_white", coloraxis2=dict(showscale=False, colorscale=["black", "black"]),
                          width=i*200+100, height=j*200+50, title=label)
        fig.update_layout(coloraxis1=dict(showscale=False), font_size=font_size)
        fig.write_image(save_output_here + "/correlation_plots_pearson_"+ label + "_" +file_type + "_" + column_indicator + ".pdf")
        fig.show()

    elif mode=="heatmap":
        da = np.ones((len(df_sub.columns), len(df_sub.columns)))
        for i,ni in enumerate(df_sub.columns):
            for j,nj in enumerate(df_sub.columns):
                # filter data and store correlation values
                df_ij = df_sub[[ni,nj]].dropna() if i!=j else pd.DataFrame(df_sub[ni].dropna())
                if i!=j:
                    da[i,j] = np.round(np.min(correlation_function(df_sub[[ni,nj]].dropna())),4)
        # create figure and label axes
        fig = go.Figure(data=go.Heatmap(z=da))
        fig.update_xaxes(tickvals=list(range(0,i+1,1)),
                          ticktext=list(df_sub.columns))
        fig.update_yaxes(tickvals=list(range(0,j+1,1)),
                          ticktext=list(df_sub.columns))
        fig.update_layout(template="simple_white", width=i*50+100, height=j*50+100, title="pearson correlation")
        fig.update_layout(coloraxis1=dict(showscale=False), font_size=font_size, width=500, height=500)
        fig.write_image(save_output_here + "/correlation_plots_heatmap_"+ file_type + "_" + column_indicator + ".pdf")
        fig.show()
    else:
        raise ValueError
    return fig


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
                "CV_"+column_indicator: "CV", "missing_value_"+column_indicator: "missing value ratio"
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
        "CV_"+column_indicator: "CV", "missing_value_"+column_indicator: "missing value ratio"
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
        "CV_"+column_indicator: "CV", "missing_value_"+column_indicator: "missing value ratio"
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
            set_temp = set(df_temp["unique_ID_column"])
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

def filter_DIANN_output_for_column_indicator(df):
    column_indicator = "Intensity"
    selected_columns = [col for col in df.columns if ".d" in col]

    ind = 0
    new_selected_columns = []
    for intensity_column in selected_columns:
        str_new = column_indicator+"_"+str(ind+1)
        new_selected_columns.append(str_new)
        df[str_new] = df[intensity_column]
        ind += 1
    print("Column_names changed from:")
    print(selected_columns)
    print("to:")
    print(new_selected_columns)
    return df, column_indicator, new_selected_columns


def filter_alphadia_output_for_column_indicator(df):
    column_indicator = "Intensity"
    selected_columns = [col for col in df.columns if (not "pg" in col) & (not 'mod_seq_charge_hash' in col) & (not 'precursor' in col)]
    for column in selected_columns:
        df[column][df[column] == 0] = np.nan

    ind = 0
    new_selected_columns = []
    for intensity_column in selected_columns:
        str_new = column_indicator+"_"+str(ind+1)
        new_selected_columns.append(str_new)
        df[str_new] = df[intensity_column]
        ind += 1
    print("Column_names changed from:")
    print(selected_columns)
    print("to:")
    print(new_selected_columns)
    return df, column_indicator, new_selected_columns


def filter_spectronaut_output_for_column_indicator(df):
    column_indicator = "Intensity"
    if "EG.PrecursorId" in df.columns:
        selected_columns = [col for col in df.columns if ("[" in col) & ("PG.Quantity" not in col)]
    else:
        selected_columns = [col for col in df.columns if "[" in col]
    ind = 0
    new_selected_columns = []
    for intensity_column in selected_columns:
        str_new = column_indicator+"_"+str(ind+1)
        new_selected_columns.append(str_new)
        df[str_new] = df[intensity_column]
        ind += 1
    print("Column_names changed from:")
    print(selected_columns)
    print("to:")
    print(new_selected_columns)
    return df, column_indicator, new_selected_columns

def filter_tables_and_make_correlation_plots(
    df_filtered,
    dict_results,
    column_indicator,
    percentage_of_missing_values,
    selected_columns,
    unique_ID_column,
    save_output_here,
    file_type,
    label,
    list_dfs_only_with_intensity_columns,
    list_dfs,
    file_name

):
    # Here starts the general part, which is not specific for MaxQuant
    dict_results, df_prefiltered = filter_results_according_to_column_indicator(
        df_filtered,
        dict_results,
        column_indicator,
        percentage_of_missing_values,
        selected_columns,
        unique_ID_column
    )
    # pearson correlation plots
    error_files = []
    try:
        common_string = find_common_string(selected_columns)

        # function form Julia Patricia Schessner, Eugenia Voytik, Isabell Bludau,  https://doi.org/10.1002/pmic.202100103
        cross_corr = plot_sample_correlations(
            df_prefiltered,
            save_output_here,
            file_type,
            column_indicator,
            label,
            data_columns=common_string+str(list(range(len(selected_columns)+1))),
        #     font_size=12
        )
    except Exception as e:
        print(e)
        error_files.append(file_name)
        print("error processing: {}".format(file_name))
        print("Did you name the experiments ending with 1, 2, 3, 4, ...?")

    df_prefiltered_mod = df_prefiltered[selected_columns].add_prefix(label)
    df_prefiltered_mod["unique_ID_column"] = df_prefiltered[unique_ID_column]

    list_dfs_only_with_intensity_columns.append(df_prefiltered_mod)
    list_dfs.append(df_prefiltered)

    return dict_results, list_dfs_only_with_intensity_columns, list_dfs
