import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

from typing import List, Tuple, Dict


# importing components for visualization
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams.update({'font.size': 16.5})
mpl.rcParams['pdf.fonttype'] = 42


def calculate_precursor_within_scan_area(
    library: pd.DataFrame,
    mz: tuple,
) -> dict:
    """This function calculates the number of precursors within the scan area
    in percent.

    Parameters:
    library (pd.DataFrame): a pre-filtered data frame with unified column names
        containing all required precursor information.
    mz (tuple): lower and upper value of the m/z range.
    im (tuple): lower and upper value of the ion mobility range.

    Returns:
    dict: this dictionary has only one item. The key is the value description,
        and the value is the ratio of precursors [%] within the set scan area.
    """
    filtered_selection = library['mz'] > mz[0]
    filtered_selection &= library['mz'] < mz[1]

    prec_scan_area = (len(library[filtered_selection])/len(library['mz']))*100

    return {'precursors within m/z-range [%]': round(prec_scan_area, 2)}


def load_method_file(filename: str) -> Tuple[pd.DataFrame, List[List[float]]]:
    """Load DIA method file and extract window metrics and bins.
    
    Parameters:
        filename: Path to method file containing window information in any of these formats:
            - Thermo format: 'm/z' and 'Isolation Window (m/z)' columns
            - Center mass format: 'Center Mass (m/z)' and 'Window Width (m/z)' columns
            - Range format: 'm/z range' column with 'start-end' format
        
    Returns:
        Tuple[pd.DataFrame, List[List[float]]]: Tuple containing:
            - df_window DataFrame with columns: window_start, window_end, middle_of_window, window_width
            - bins as list of [start, end] ranges
    """
    df = pd.read_csv(filename)
    columns = df.columns.tolist()
    
    # Check format based on column names
    if 'm/z' in columns and 'Isolation Window (m/z)' in columns:
        # Thermo targeted format
        df_window = pd.DataFrame({
            'middle_of_window': df['m/z'],
            'window_width': df['Isolation Window (m/z)']
        })
        df_window['window_start'] = df_window['middle_of_window'] - df_window['window_width']/2
        df_window['window_end'] = df_window['middle_of_window'] + df_window['window_width']/2
        
    elif 'Center Mass (m/z)' in columns and 'Window Width (m/z)' in columns:
        # Center mass format
        df_window = pd.DataFrame({
            'middle_of_window': df['Center Mass (m/z)'],
            'window_width': df['Window Width (m/z)']
        })
        df_window['window_start'] = df_window['middle_of_window'] - df_window['window_width']/2
        df_window['window_end'] = df_window['middle_of_window'] + df_window['window_width']/2
        
    elif 'm/z range' in columns:
        # Range format
        # Parse the range strings, handling potential whitespace
        ranges = [range_str.strip().split('-') for range_str in df['m/z range']]
        df_window = pd.DataFrame({
            'window_start': [float(r[0].strip()) for r in ranges],
            'window_end': [float(r[1].strip()) for r in ranges]
        })
        df_window['middle_of_window'] = (df_window['window_start'] + df_window['window_end']) / 2
        df_window['window_width'] = df_window['window_end'] - df_window['window_start']
        
    else:
        raise ValueError(
            f"Unrecognized file format. File must contain either:\n"
            f"- 'm/z' and 'Isolation Window (m/z)' columns, or\n"
            f"- 'Center Mass (m/z)' and 'Window Width (m/z)' columns, or\n"
            f"- 'm/z range' column\n"
            f"Found columns: {columns}"
        )
    
    # Create bins list from window_start and window_end
    bins = [[start, end] for start, end in zip(df_window['window_start'], df_window['window_end'])]
    
    # Ensure consistent column order
    df_window = df_window[['window_start', 'window_end', 'middle_of_window', 'window_width']]
    
    return df_window, bins, df


def analyze_bins(
    bins: List[List[float]], 
    values: List[float],
    max_width: float = 50.0,
    min_width: float = 2.0,
) -> Dict:
    """Analyzes the distribution of items across bins and generates comprehensive
    statistics including a detailed tabular view.

    Parameters:
        bins (List[List[float]]): List of [start, end] bin ranges
        values (List[float]): Original values that were binned
        max_width (float): Maximum allowed width for validation

    Returns:
        Dict: Dictionary containing bin statistics and detailed analysis:
            - num_bins: Total number of bins
            - min_width: Minimum bin width
            - max_width: Maximum bin width
            - avg_width: Average bin width
            - min_count: Minimum items in a bin
            - max_count: Maximum items in a bin
            - avg_count: Average items per bin
            - total_width: Sum of all bin widths
            - coverage: Ratio of total bin width to data range
            - bins: Original bin ranges
            - tabular_stats: Detailed bin-by-bin analysis in tabular format
    """
    # Step 1: Calculate statistics for each individual bin
    bin_stats = []
    widths = []
    counts = []
    
    for i, bin_range in enumerate(bins):
        # Calculate basic bin metrics
        start, end = bin_range
        width = end - start
        count = sum(1 for v in values if start <= v <= end)
        width_ok = width <= max_width and width >= min_width
        
        # Store for summary statistics
        widths.append(width)
        counts.append(count)
        
        # Store detailed bin info
        bin_stats.append({
            'bin': i,
            'start': start,
            'end': end,
            'width': width,
            'items': count,
            'width_ok': width_ok
        })

    # Step 2: Create formatted table view
    table_output = "Bin statistics:\n"
    table_output += "Bin | " + "Range".ljust(25) + " | " + "Width".ljust(8) + " | " + "Items".ljust(7) + " | Width OK\n"
    table_output += "-" * 75 + "\n"
    
    for stat in bin_stats:
        range_str = f"{stat['start']:.1f}- {stat['end']:.1f}"
        table_output += f"{stat['bin']:3d} | {range_str:23s} | {stat['width']:6.2f} | {stat['items']:5d} | {stat['width_ok']}\n"

    # Step 3: Compile all statistics into final report
    stats = {
        'num_bins': len(bins),
        'min_width': min(widths),
        'max_width': max(widths),
        'avg_width': sum(widths) / len(widths),
        'min_count': min(counts),
        'max_count': max(counts),
        'avg_count': sum(counts) / len(counts),
        'total_width': sum(widths),
        'bins': bins,
        'bin_details': bin_stats,      # Detailed stats for each bin
        'tabular_stats': table_output  # Formatted table view
    }
    
    return stats


def parse_stats_text(text):
    """Convert stats table text into a pandas DataFrame"""
    
    # Split by newlines and get only data lines
    lines = text.strip().split('\n')
    data_lines = [line for line in lines if line and '|' in line and not line.startswith('-')]
    data_lines = data_lines[1:]  # Skip header
    
    # Initialize lists to store data
    bins = []
    starts = []
    ends = []
    widths = []
    items_list = []
    
    # Parse each line
    for line in data_lines:
        # Split by | and clean whitespace
        parts = [p.strip() for p in line.split('|')]
        
        # Parse bin number
        bins.append(int(parts[0]))
        
        # Parse range
        range_parts = parts[1].strip().split('-')
        starts.append(float(range_parts[0]))
        ends.append(float(range_parts[1]))
        
        # Parse width and items
        widths.append(float(parts[2]))
        items_list.append(int(parts[3]))
    
    # Create DataFrame
    df = pd.DataFrame({
        'bin': bins,
        'start': starts,
        'end': ends,
        'width': widths,
        'items': items_list
    })
    
    return df


def plot_precursors_per_scan(window_type, df_window, file_name, gui=False):
    """
    Function to plot precursors per scan with bar color indications for scan width.
    
    Arguments:
    df_window: DataFrame containing columns 'start', 'end', 'width', and 'items'
    window_type: string indicating window type ('dynamic' or 'fixed')
    file_name: string for output file path
    gui: boolean to determine if figure should be returned
    
    Returns: matplotlib figure if gui=True, else None
    """
    # Calculate middle points for each window
    df_window['middle_of_window'] = (df_window['start'] + df_window['end']) / 2
    
    # Normalize your data to 0-1 range
    norm = plt.Normalize(min(df_window['width']), max(df_window['width']))
    
    # Setting color map
    cmap = plt.cm.viridis
    colors = cmap(norm(df_window['width']))
    
    fig, ax = plt.subplots()
    
    if window_type == "dynamic":
        bars = ax.bar(df_window['middle_of_window'], 
                     df_window['items'],  # Using pre-calculated items
                     width=df_window['width'], 
                     color=colors,
                     edgecolor='black', 
                     linewidth=0.25)
        # Create colorbar
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, orientation='vertical')
        cb.set_label('scan width')
        
    elif window_type == "fixed":
        bars = ax.bar(df_window['middle_of_window'], 
                     df_window['items'],  # Using pre-calculated items
                     width=df_window['width'], 
                     color=[0.127568, 0.566949, 0.550556, 1.],
                     edgecolor='black', 
                     linewidth=0.25)
    
    plt.xlabel('m/z range [Da]')
    plt.ylabel('# precursors')
    
    plt.savefig(file_name, bbox_inches='tight', pad_inches=0, dpi=300)
    
    if gui:
        return fig
    else:
        plt.clf()




