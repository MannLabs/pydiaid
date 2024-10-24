## import public packages
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# from matplotlib import rc

import pandas as pd
# from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib.ticker import NullFormatter
# from tqdm import tqdm

from typing import List, Tuple, Dict

## import os for setting of directory
import os

# importing components for visualization
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams.update({'font.size': 16.5})
mpl.rcParams['pdf.fonttype'] = 42


def create_method(
    mz_values: List[float],
    mz_range: Tuple[float, float],
    num_bins: int,
    window_type: str,
    output_format: str = "all",
    folder_path: str = ".",
    min_width: float = 2.0,
    max_width: float = 50.0,
    base_name: str = "OA_DIA_method"
) -> Tuple[pd.DataFrame, List[List[float]]]:
    """Create DIA method with specified window type and output format.
    
    Parameters:
        mz_values: List of m/z values to bin
        mz_range: Tuple of (min_mz, max_mz)
        num_bins: Number of desired bins
        window_type: Type of windowing ("fixed" or "variable")
        output_format: Type of output files to create ("all", "targeted", 
                      "center_mass", or "mz_ranges")
        folder_path: Directory to save output files
        min_width: Minimum allowed bin width (for variable windows)
        max_width: Maximum allowed bin width (for variable windows)
        base_name: Base name for output files
        
    Returns:
        Tuple[pd.DataFrame, List[List[float]]]: Window metrics DataFrame and bins list
        
    Raises:
        ValueError: If window_type or output_format is not supported
    """
    # Validate inputs
    """creates bins that have an approximately equal number of items in 
    each bin, while respecting the min/max width constraints and the total 
    number of bins constraint."""
    if window_type not in ["fixed", "variable"]:
        raise ValueError(f"Window type '{window_type}' not supported. Use 'fixed' or 'variable'")
        
    if output_format not in ["all", "targeted", "center_mass", "mz_ranges"]:
        raise ValueError(
            f"Output format '{output_format}' not supported. "
            "Use 'all', 'targeted', 'center_mass', or 'mz_ranges'"
        )
    
    # Create bins based on window type
    if window_type == "variable":
        bins = create_variable_bins(
            mz_values,
            mz_range,
            num_bins,
            min_width,
            max_width,
            merging_for_small_bins=True
        )
    else:  # fixed
        bins = create_fixed_bins(
            mz_values,
            mz_range,
            num_bins
        )
    
    # Calculate window metrics
    df_window = calculate_window_metrics(bins)
    
    # Generate base method name including parameters
    base_method_name = f"{base_name}_{mz_range[0]}_{mz_range[1]}_{num_bins}_{window_type}"
    
    # Create output files based on format
    if output_format == "all":
        create_all_formats(
            df_window,
            os.path.join(folder_path, base_method_name)
        )
    elif output_format == "targeted":
        create_targeted_mass_list(
            df_window,
            os.path.join(folder_path, f"{base_method_name}_targeted.csv")
        )
    elif output_format == "center_mass":
        create_center_mass_list(
            df_window,
            os.path.join(folder_path, f"{base_method_name}_center_mass.csv")
        )
    elif output_format == "mz_ranges":
        create_mz_range_list(
            df_window,
            os.path.join(folder_path, f"{base_method_name}_mz_ranges.csv")
        )
    
    return df_window, bins


def get_initial_equal_bins(
    values: List[float], 
    num_splits: int
) -> List[List[float]]:
    """Creates initial bins with equal number of items by splitting the data into
    roughly equal chunks.

    Parameters:
        values (List[float]): List of m/z values to be binned
        num_splits (int): Number of desired splits/bins

    Returns:
        List[List[float]]: List of [start, end] ranges for each bin
    """
    # If data can't be evenly split, pad with copies of last value
    # Example: if 10 items and 3 splits wanted -> add 2 copies of last value
    remainder = len(values) % num_splits
    if remainder:
        values = values + [values[-1]] * (num_splits - remainder)
    
    # Split data into roughly equal-sized chunks using numpy
    # Example: [1,2,3,4,5,6] with 3 splits -> [[1,2], [3,4], [5,6]]
    splits = np.array_split(np.array(values), num_splits)
    
    # Convert splits into [start, end] ranges
    # Example: [[1,2], [3,4], [5,6]] -> [[1,2], [2,4], [4,6]]
    ranges = []
    for i, split in enumerate(splits):
        if i == 0:
            # First range starts with first value in first split
            ranges.append([split[0], split[-1]])
        else:
            # Subsequent ranges start where previous ended
            ranges.append([ranges[-1][1], split[-1]])
    return ranges


def split_oversized_bin_fixed_width(
    bin_range: List[float], 
    max_width: float, 
    mz_range: Tuple[float, float], 
    remaining_bins: int
) -> Tuple[List[List[float]], float]:
    """Splits an oversized bin into fixed-width bins while considering the number
    of remaining target bins.

    Parameters:
        bin_range (List[float]): [start, end] range of the bin to split
        max_width (float): Maximum allowed width for any bin
        mz_range (Tuple[float, float]): Overall (min_mz, max_mz) range
        remaining_bins (int): Number of bins still available for use

    Returns:
        Tuple[List[List[float]], float]: List of new bins and the new start position
    """
    start, end = bin_range
    new_bins = []
    total_width = end - start
    
    # Special case: if only one bin left and width is close to max, keep as is
    if remaining_bins == 1 and total_width <= max_width * 1.2:  # 20% tolerance
        new_bins.append([start, end])
        return new_bins, end
    
    # Split into fixed-width segments
    current_start = start
    while current_start + max_width <= end and current_start + max_width <= mz_range[1]:
        # Create a new bin of exactly max_width
        new_bins.append([current_start, current_start + max_width])
        current_start += max_width
        
        # Stop if we've used all available bin slots
        if len(new_bins) >= remaining_bins:
            break
    
    return new_bins, current_start


def merge_small_bins(
    bins: List[List[float]], 
    min_width: float
) -> List[List[float]]:
    """Merges bins that are smaller than the minimum width with their neighbors.

    Parameters:
        bins (List[List[float]]): List of [start, end] bin ranges
        min_width (float): Minimum allowed width for any bin

    Returns:
        List[List[float]]: List of merged bin ranges
    """
    print("\nStarting bin merging process...")
    print(f"Minimum allowed width: {min_width}")
    
    # Process bins sequentially, merging small ones with their next neighbor
    merged = []
    i = 0
    while i < len(bins):
        current = bins[i]
        width = current[1] - current[0]
        
        print(f"\nChecking bin {i}: [{current[0]:.2f}, {current[1]:.2f}] (width: {width:.2f})")
        
        # If bin is too small and not the last bin, merge with next bin
        if width < min_width and i + 1 < len(bins):
            next_bin = bins[i+1]
            new_bin = [current[0], next_bin[1]]
            print(f"Width {width:.2f} is below minimum {min_width}")
            print(f"Merging with next bin [{next_bin[0]:.2f}, {next_bin[1]:.2f}]")
            print(f"Result: [{new_bin[0]:.2f}, {new_bin[1]:.2f}]")
            merged.append(new_bin)
            i += 2  # Skip next bin since we used it in merge
        else:
            print(f"Width {width:.2f} is acceptable, keeping bin as is")
            merged.append(current)
            i += 1
            
    print(f"\nFinal number of bins after merging: {len(merged)}")
    return merged


def create_variable_bins(
    values: List[float], 
    mz_range: Tuple[float, float],
    num_bins: int,
    min_width: float = 2.0,
    max_width: float = 50.0,
    merging_for_small_bins: bool = False
) -> List[List[float]]:
    """Creates variable bins for mass spectrometry data that balance even distribution
    of items with reasonable bin widths.

    Parameters:
        values (List[float]): List of m/z values to bin
        mz_range (Tuple[float, float]): (min_mz, max_mz) range to consider
        num_bins (int): Target number of bins
        min_width (float): Minimum allowed bin width in Da
        max_width (float): Maximum allowed bin width in Da
        merging_for_small_bins (bool): Defines if only max_width contrain or both, 
            min_width and max_width, is applied

    Returns:
        List[List[float]]: List of [start, end] ranges for each bin
    """
    # Filter values to specified range and sort them
    filtered_values = sorted([v for v in values if mz_range[0] <= v <= mz_range[1]])
    if not filtered_values:
        return []

    # Initialize tracking variables
    remaining_range = list(mz_range)  # Track what part of range still needs binning
    final_bins = []                   # Store our final bin list
    remaining_bins = num_bins         # Track how many bins we can still create
    
    # Main binning loop - continues until we use all bins or cover whole range
    while remaining_bins > 0 and remaining_range[0] < remaining_range[1]:
        # Get values in current remaining range (m/z)
        current_values = [v for v in filtered_values if remaining_range[0] <= v <= remaining_range[1]]
        if not current_values:
            break
        
        # Create initial equal-sized bins for remaining range
        bins = get_initial_equal_bins(current_values, remaining_bins)
        if not bins:
            break
        
        # Process each proposed bin
        for bin_range in bins:
            width = bin_range[1] - bin_range[0]
            
            if width > max_width:
                # Bin too wide - split it into fixed-width pieces
                new_bins, new_start = split_oversized_bin_fixed_width(
                    bin_range, max_width, mz_range, remaining_bins
                )
                final_bins.extend(new_bins)
                remaining_bins -= len(new_bins)
                remaining_range[0] = new_start  # Update where to start next iteration
                break  # Recalculate remaining range with new bins
            else:
                # Bin width acceptable - keep it
                final_bins.append(bin_range)
                remaining_bins -= 1
        
        # Safety check to prevent infinite loops
        if remaining_range[0] >= remaining_range[1]:
            break

    if merging_for_small_bins is True:
        # Final step: merge any bins that ended up too small
        final_bins = merge_small_bins(final_bins, min_width)
    return final_bins


def get_initial_equal_width_bins(
    values: List[float],
    num_splits: int
) -> List[List[float]]:
    """Creates initial bins with equal m/z width by splitting the total range into
    equal segments.

    Parameters:
        values (List[float]): List of m/z values to be binned
        num_splits (int): Number of desired splits/bins

    Returns:
        List[List[float]]: List of [start, end] ranges for each bin
    """
    # Find the total range to split
    min_val = min(values)
    max_val = max(values)
    total_range = max_val - min_val
    
    # Calculate the width of each bin
    bin_width = total_range / num_splits
    
    # Create the ranges
    ranges = []
    for i in range(num_splits):
        start = min_val + (i * bin_width)
        end = min_val + ((i + 1) * bin_width)
        
        # For the last bin, ensure we include the max value
        if i == num_splits - 1:
            end = max_val
            
        ranges.append([start, end])
    
    return ranges

def create_fixed_bins(
    values: List[float], 
    mz_range: Tuple[float, float],
    num_bins: int,
) -> List[List[float]]:
    """Creates variable bins for mass spectrometry data that balance even distribution
    of items with reasonable bin widths.

    Parameters:
        values (List[float]): List of m/z values to bin
        mz_range (Tuple[float, float]): (min_mz, max_mz) range to consider
        num_bins (int): Target number of bins
        min_width (float): Minimum allowed bin width in Da
        max_width (float): Maximum allowed bin width in Da

    Returns:
        List[List[float]]: List of [start, end] ranges for each bin
    """
    # Filter values to specified range and sort them
    filtered_values = sorted([v for v in values if mz_range[0] <= v <= mz_range[1]])
    if not filtered_values:
        return []
    
    return get_initial_equal_width_bins(
        filtered_values,
        num_bins
    )


def calculate_window_metrics(bins: List[List[float]]) -> pd.DataFrame:
    """Calculate window metrics from bin ranges.
    
    Parameters:
        bins (List[List[float]]): List of [start, end] m/z ranges
        
    Returns:
        pd.DataFrame: DataFrame containing:
            - window_start: Start m/z of each window
            - window_end: End m/z of each window
            - middle_of_window: Center m/z of each window
            - window_width: Width of each window
    """
    metrics = {
        'window_start': [],
        'window_end': [],
        'middle_of_window': [],
        'window_width': []
    }
    
    for start, end in bins:
        metrics['window_start'].append(start)
        metrics['window_end'].append(end)
        metrics['middle_of_window'].append((start + end) / 2)
        metrics['window_width'].append(end - start)
        
    return pd.DataFrame(metrics)


def create_targeted_mass_list(
    df_window: pd.DataFrame,
    filename: str,
    rt_time: str = '16.5',
    window: str = '33',
    charge: str = '2'
) -> None:
    """Creates targeted mass list format. This function uses your default values 
    for RT Time (min) and Window (min). If you want to use different 
    values, pass them as arguments: 
    create_dia_scheme(df_window, "OA_DIA_method_targeted.csv", '17', '34')
    
    Parameters:
        df_window: DataFrame with window metrics
        filename: Output filename
        rt_time: RT Time in minutes
        window: Window in minutes
        charge: Default charge state
    """
    method_df = pd.DataFrame({
        'Compound': range(1, len(df_window) + 1),
        'Formula': '',
        'Adduct': '',
        'm/z': df_window['middle_of_window'],
        'z': charge,
        'RT Time (min)': rt_time,
        'Window (min)': window,
        'Isolation Window (m/z)': df_window['window_width']
    })
    
    columns = [
        'Compound', 'Formula', 'Adduct', 'm/z', 'z',
        'RT Time (min)', 'Window (min)', 'Isolation Window (m/z)'
    ]
    method_df[columns].to_csv(filename, index=False)


def create_center_mass_list(
    df_window: pd.DataFrame,
    filename: str
) -> None:
    """Creates a center mass list format.
    
    Parameters:
        df_window: DataFrame with window metrics
        filename: Output filename
    """
    center_mass_df = pd.DataFrame({
        'Center Mass (m/z)': df_window['middle_of_window'],
        'Window Width (m/z)': df_window['window_width']
    })
    center_mass_df.to_csv(filename, index=False)


def create_mz_range_list(
    df_window: pd.DataFrame,
    filename: str
) -> None:
    """Creates an m/z range list format.
    
    Parameters:
        df_window: DataFrame with window metrics
        filename: Output filename
    """
    # Format each range with proper spacing
    ranges = [
        f" {start:.7f}-{end:.7f}" 
        for start, end in zip(df_window['window_start'], df_window['window_end'])
    ]
    
    range_df = pd.DataFrame({
        'm/z range': ranges
    })
    range_df.to_csv(filename, index=False)


def create_all_formats(
    df_window: pd.DataFrame,
    base_filename: str,
    rt_time: str = '16.5',
    window: str = '33',
    charge: str = '2'
) -> None:
    """Creates all three format types from bin ranges.
    
    Parameters:
        df_window: DataFrame with window metrics
        base_filename: Base name for output files (without extension)
        rt_time: RT Time in minutes
        window: Window in minutes
        charge: Default charge state
    """
    
    # Create all three formats
    create_targeted_mass_list(
        df_window, 
        f"{base_filename}_targeted.csv",
        rt_time,
        window,
        charge
    )
    create_center_mass_list(
        df_window,
        f"{base_filename}_center_mass.csv"
    )
    create_mz_range_list(
        df_window,
        f"{base_filename}_mz_ranges.csv"
    )



