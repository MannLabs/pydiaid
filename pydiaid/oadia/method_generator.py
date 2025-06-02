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


import numpy as np
from typing import List, Tuple, Dict
import math


def create_binning_with_constraints(
    mz_values: List[float],
    mz_range: Tuple[float, float],
    num_bins: int,
    max_width: float,
    max_iterations: int = 5
) -> List[List[float]]:
    """
    Creates bins for mass spectrometry data with constraints on bin width.
    
    Args:
        mz_values: List of m/z values to bin
        mz_range: (min_mz, max_mz) range to consider
        num_bins: Target number of bins
        max_width: Maximum allowed bin width
        max_iterations: Maximum number of refinement iterations
        
    Returns:
        List of [start, end] ranges for each bin
        
    Raises:
        ValueError: If constraints cannot be satisfied
    """
    # Check if basic constraints can be satisfied
    if (mz_range[1] - mz_range[0])/num_bins > max_width:
        raise ValueError(
            "Not all constraints can be fulfilled. Either choose a smaller m/z-range, "
            "more bins, or a wider maximal width."
        )
    
    # Initialize bins list
    bins_extended = list()
    
    # Create initial binning
    bins_temp, oversized_indices = create_initial_bins_and_check_for_overlength(
        mz_values,
        mz_range,
        num_bins,
        max_width
    )
    
    # Handle initial binning result
    if len(oversized_indices) == 0:
        bins_extended.append(bins_temp)
    else:
        bins_extended = reduce_oversized_bins(
            mz_range,
            max_width,
            bins_temp,
            oversized_indices,
            bins_extended
        )
    
    # Iterative refinement
    iteration_count = 0
    while iteration_count < max_iterations:
        print(iteration_count)
        bins_extended = flatten_and_sort_bins(bins_extended)
        
        # Get new ranges and bin numbers
        new_mz_ranges, new_num_bins = redistribute_bins_general(
            mz_values,
            mz_range,
            num_bins,
            bins_extended,
            max_width
        )
        
        # Process each range
        all_ranges_processed = True
        for i in range(len(new_mz_ranges)):
            bins_temp, oversized_indices = create_initial_bins_and_check_for_overlength(
                mz_values,
                new_mz_ranges[i],
                new_num_bins[i],
                max_width
            )
            
            if len(oversized_indices) == 0:
                bins_extended.append(bins_temp)
            else:
                all_ranges_processed = False
                bins_extended = reduce_oversized_bins(
                    new_mz_ranges[i],
                    max_width,
                    bins_temp,
                    oversized_indices,
                    bins_extended
                )
        
        # Check if we're done
        if all_ranges_processed:
            break
            
        iteration_count += 1
        
        if iteration_count == max_iterations:
            print("Warning: Reached maximum number of iterations")
    
    # Final cleanup and sorting
    return flatten_and_sort_bins(bins_extended)


def reduce_oversized_bins(
    mz_range,
    max_width,
    bins_temp, 
    oversized_indices,    
    bins_extended,
):

    categorized_bins = categorize_and_merge_oversized_bins(
        bins_temp, 
        oversized_indices, 
        mz_range
    )
    print("categoraized", categorized_bins)
    new_bins = transfer_oversized_bins_to_bins_with_maximal_widths(
        categorized_bins, 
        max_width
    )
    print("new_bins", new_bins)
    bins_extended.extend(new_bins) 
    print("current bins", bins_extended, len(bins_extended))

    return bins_extended


def redistribute_bins_general(
    mz_values, 
    mz_range,
    num_bins,
    bins_extended, 
    max_width    
):
    new_num_bins = num_bins - len(bins_extended)
    print("num bins for next iteration", new_num_bins)

    new_mz_ranges = find_gaps(bins_extended, mz_range)
    print("new mz ranges", new_mz_ranges)

    items_per_range, weights, bins_per_range_rounded = calculate_bins_per_range(
        mz_values,
        new_mz_ranges,
        new_num_bins
    )
    print("bins_per_range_rounded", bins_per_range_rounded)

    min_num_bins = list()
    for num in range(len(new_mz_ranges)):
        min_num_bins.append(math.ceil((new_mz_ranges[num][1]-new_mz_ranges[num][0])/max_width))
    print("min num bins", min_num_bins)

    new_bin_num = redistribute_bins(bins_per_range_rounded, min_num_bins)
    print("new num bins", new_bin_num)

    return new_mz_ranges, new_bin_num


def categorize_and_merge_oversized_bins(
    bins: List[List[float]], 
    oversized_indices: List[int], 
    mz_range: Tuple[float, float]
) -> Dict[str, List[Dict[str, any]]]:
    """
    Categorizes oversized bins by their position relative to the mz range midpoint
    and merges adjacent oversized bins.
    
    Parameters:
    bins (List[List[float]]): List of [start, end] ranges for each bin
    oversized_indices (List[int]): Indices of bins that exceed max_width
    mz_range (Tuple[float, float]): (min_mz, max_mz) range to consider
    
    Returns:
    Dict with keys 'left_region' and 'right_region', each containing a list of
    dictionaries with merged bin information
    """
    if not oversized_indices:
        return {'left_region': [], 'right_region': []}
    
    # Calculate midpoint of mz range
    mz_midpoint = (mz_range[0] + mz_range[1]) / 2
    
    # Initialize categories
    left_region = []
    right_region = []
    
    # Group adjacent oversized bins
    grouped_indices = []
    current_group = [oversized_indices[0]]
    
    for i in range(1, len(oversized_indices)):
        if oversized_indices[i] == oversized_indices[i-1] + 1:
            current_group.append(oversized_indices[i])
        else:
            grouped_indices.append(current_group)
            current_group = [oversized_indices[i]]
    grouped_indices.append(current_group)
    
    # Process each group of adjacent oversized bins
    for group in grouped_indices:
        start_idx = group[0]
        end_idx = group[-1]
        
        merged_bin = {
            'mz_range': [bins[start_idx][0], bins[end_idx][1]],
            'bin_width': bins[end_idx][1] - bins[start_idx][0],
        }
        
        # Calculate center of merged bin
        bin_center = (merged_bin['mz_range'][0] + merged_bin['mz_range'][1]) / 2
        
        # Categorize based on position relative to midpoint
        if bin_center < mz_midpoint:
            left_region.append(merged_bin)
        else:
            right_region.append(merged_bin)
    
    return {
        'left_region': left_region,
        'right_region': right_region
    }


def transfer_oversized_bins_to_bins_with_maximal_widths(
    categorized_bins: Dict[str, List[Dict[str, any]]], 
    max_width: float
) -> Dict[str, Dict[str, any]]:
    """
    Generates new windows from oversized regions that respect max_width.
    For left region, starts from lower m/z; for right region, starts from higher m/z.
    
    Parameters:
    categorized_bins (Dict): Output from categorize_and_merge_oversized_bins
    max_width (float): Maximum allowed window width
    
    Returns:
    Dict containing new windows, remaining ranges, and bin counts for each region
    """

    new_bins = list()
    # new_mz_range = list()
    
    # Process left region (start from lower m/z)
    for merged_bin in categorized_bins['left_region']:
        start_mz = merged_bin['mz_range'][0]
        end_mz = merged_bin['mz_range'][1]
        current_mz = start_mz
        
        while current_mz + max_width < end_mz:
            new_window = [current_mz, current_mz + max_width]
            new_bins.append(new_window)
            current_mz += max_width
    
    # Process right region (start from higher m/z)
    for merged_bin in categorized_bins['right_region']:
        start_mz = merged_bin['mz_range'][0]
        end_mz = merged_bin['mz_range'][1]
        current_mz = end_mz
        
        while current_mz - max_width > start_mz:
            new_window = [current_mz - max_width, current_mz]
            new_bins.append(new_window) # Insert at beginning to maintain order
            current_mz -= max_width

    return sorted(new_bins)
    
    
def find_gaps(intervals, mz_range):
    """
    Find gaps between intervals in a list of [start, end] intervals.

    Args:
        intervals: List of [start, end] intervals
        
    Returns:
        List of tuples (gap_start, gap_end) representing gaps between intervals
    """
    # Sort intervals by start time
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    sorted_intervals.insert(0, [mz_range[0], mz_range[0]])
    sorted_intervals.append([mz_range[1], mz_range[1]])
    gaps = []


    # Compare adjacent intervals
    for i in range(len(sorted_intervals) - 1):
        current_end = sorted_intervals[i][1]
        next_start = sorted_intervals[i + 1][0]
        
        # If there's a gap between current end and next start
        if current_end < next_start:
            gaps.append((current_end, next_start))

    return gaps


def calculate_bins_per_range(
    mz_values: List[float],
    mz_ranges: List[Tuple[float, float]],
    total_num_bins: int
) -> Tuple[List[int], List[float], List[int]]:
    """
    Calculate the optimal number of bins for each m/z range based on data density.
    
    Args:
        mz_values: List of m/z values
        mz_ranges: List of tuples defining (start, end) of each m/z range
        total_num_bins: Total number of bins to distribute
        
    Returns:
        Tuple containing:
        - List of items per range
        - List of weights for each range
        - List of number of bins per range (rounded to integers)
    """
    # Calculate number of items per m/z range
    items_per_range = []
    for range_start, range_end in mz_ranges:
        filtered_values = [v for v in mz_values if range_start <= v <= range_end]
        items_per_range.append(len(filtered_values))
    
    # Calculate weights based on item density
    total_items = sum(items_per_range)
    weights = [count / total_items for count in items_per_range]
    
    # Calculate number of bins per range
    bins_per_range = [weight * total_num_bins for weight in weights]
    
    # Round to nearest integer
    bins_per_range_rounded = np.round(bins_per_range).astype(int)
    
    # Adjust if rounding caused total to be off
    total_assigned = sum(bins_per_range_rounded)
    if total_assigned != total_num_bins:
        diff = total_num_bins - total_assigned
        # Add/subtract difference to/from range with highest weight
        max_weight_idx = weights.index(max(weights))
        bins_per_range_rounded[max_weight_idx] += diff
    
    return items_per_range, weights, bins_per_range_rounded.tolist()


def redistribute_bins(weighted_bins, min_bins):
    """
    Redistributes bins while maintaining minimum requirements and total count.
    
    Args:
        weighted_bins (list): Initial bin distribution based on weights
        min_bins (list): Minimum number of bins required for each range
        
    Returns:
        list: Adjusted bin distribution meeting minimum requirements
    """
    # First ensure minimum requirements
    result = []
    total_original = sum(weighted_bins)
    total_min = sum(min_bins)
    
    # If total minimum required is equal to total available, just return min_bins
    if total_min >= total_original:
        return min_bins
    
    # First set all to their minimum values
    remaining_bins = total_original - total_min
    result = min_bins.copy()
    
    # Then distribute remaining bins proportionally where possible
    weights = []
    for i in range(len(weighted_bins)):
        if weighted_bins[i] > min_bins[i]:
            weights.append(weighted_bins[i])
        else:
            weights.append(0)
            
    # Normalize weights
    total_weight = sum(weights)
    if total_weight > 0:
        weights = [w/total_weight if w > 0 else 0 for w in weights]
        
        # Distribute remaining bins
        for i in range(len(result)):
            if weights[i] > 0:
                extra = round(remaining_bins * weights[i])
                result[i] += extra
    
    return result


def create_initial_bins_and_check_for_overlength(
    mz_values, 
    mz_range,
    num_bins, 
    max_width
):
    """Creates 
    Parameters:
        values (List[float]): List of m/z values to bin
        mz_range (Tuple[float, float]): (min_mz, max_mz) range to consider
        num_bins (int): Target number of bins
        max_width (float): Maximum allowed bin width in Da
    Returns:
        .....
    """
    bins_temp = create_bins_with_varibale_length(
        mz_values, 
        mz_range,
        num_bins, 
    )
    print("bins", bins_temp)
    oversized_indices = identify_oversized_bins(bins_temp, max_width)
    print("oversized bin indices", oversized_indices)
    return bins_temp, oversized_indices


def create_bins_with_varibale_length(
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
        merging_for_small_bins (bool): Defines if only max_width contrain or both, 
            min_width and max_width, is applied
    Returns:
        List[List[float]]: List of [start, end] ranges for each bin
    """
    # Filter values to specified range and sort them
    filtered_values = sorted([v for v in values if mz_range[0] <= v <= mz_range[1]])
    if not filtered_values:
        return []
    # Create initial equal-sized bins for remaining range
    bins = get_initial_equal_bins(filtered_values, num_bins)
    # replace first and last value with m/z-range limits.
    bins[0][0] = mz_range[0]
    bins[-1][1] = mz_range[1] 
    return bins


def identify_oversized_bins(bins: List[List[float]], max_width: float) -> List[int]:
    """
    Identifies bins that exceed the specified maximum width.
    Parameters:
    bins (List[List[float]]): List of [start, end] ranges for each bin
    max_width (float): Maximum allowed bin width in Da
    Returns:
    List[int]: Indices of bins that exceed max_width
    """
    oversized_bins = []
    for i, bin_range in enumerate(bins):
        bin_width = bin_range[1] - bin_range[0]
        if bin_width > max_width:
            oversized_bins.append(i)
    return oversized_bins


def flatten_and_sort_bins(bins):
    """
    Flattens and sorts a list of bins where some elements might be nested lists of bins.
    Args:
        bins (list): List of [start, end] ranges, some elements might be lists of ranges
    Returns:
        list: Flattened and sorted list of [start, end] ranges
    """
    # Initialize flattened list
    flattened = []
    # Process each element
    for item in bins:
        if len(item) == 2 and isinstance(item[0], (int, float)):
            # This is a simple [start, end] range
            flattened.append(item)
        else:
            # This is a nested list of ranges
            flattened.extend(item)
    # Sort based on start value (first element of each sublist)
    return sorted(flattened, key=lambda x: x[0])


"""-----------------------------------------------------------------------------------------"""

def create_method(
    mz_values: List[float],
    mz_range: Tuple[float, float],
    num_bins: int,
    window_type: str,
    output_format: str = "all",
    folder_path: str = ".",
    min_width: float = 2.0,
    max_width: float = 50.0,
    adjusted_for_forbidden_zones: bool = False,
    phospho_method: bool = False,
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
    if window_type not in ["fixed", "dynamic"]:
        raise ValueError(f"Window type '{window_type}' not supported. Use 'fixed' or 'dynamic'")
        
    if output_format not in ["all", "targeted", "center_mass", "mz_ranges"]:
        raise ValueError(
            f"Output format '{output_format}' not supported. "
            "Use 'all', 'targeted', 'center_mass', or 'mz_ranges'"
        )
    
    # Create bins based on window type
    if window_type == "dynamic":
        bins = create_variable_bins(
            mz_values,
            mz_range,
            num_bins,
            min_width,
            max_width,
        )
    else:  # fixed
        bins = get_initial_equal_width_bins(
            mz_range,
            num_bins
        )

    if adjusted_for_forbidden_zones is True:
        final_bins = adjust_bin_boundaries(bins, phospho_enriched=phospho_method)
    else:
        final_bins = bins
    
    # Calculate window metrics
    df_window = calculate_window_metrics(final_bins)
    
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
    
    return df_window, final_bins


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


# def merge_small_bins(
#     bins: List[List[float]], 
#     min_width: float
# ) -> List[List[float]]:
#     """Merges bins that are smaller than the minimum width with their neighbors.

#     Parameters:
#         bins (List[List[float]]): List of [start, end] bin ranges
#         min_width (float): Minimum allowed width for any bin

#     Returns:
#         List[List[float]]: List of merged bin ranges
#     """
#     print("\nStarting bin merging process...")
#     print(f"Minimum allowed width: {min_width}")
    
#     # Process bins sequentially, merging small ones with their next neighbor
#     merged = []
#     i = 0
#     while i < len(bins):
#         current = bins[i]
#         width = current[1] - current[0]
        
#         print(f"\nChecking bin {i}: [{current[0]:.2f}, {current[1]:.2f}] (width: {width:.2f})")
        
#         # If bin is too small and not the last bin, merge with next bin
#         if width < min_width and i + 1 < len(bins):
#             next_bin = bins[i+1]
#             new_bin = [current[0], next_bin[1]]
#             print(f"Width {width:.2f} is below minimum {min_width}")
#             print(f"Merging with next bin [{next_bin[0]:.2f}, {next_bin[1]:.2f}]")
#             print(f"Result: [{new_bin[0]:.2f}, {new_bin[1]:.2f}]")
#             merged.append(new_bin)
#             i += 2  # Skip next bin since we used it in merge
#         else:
#             print(f"Width {width:.2f} is acceptable, keeping bin as is")
#             merged.append(current)
#             i += 1
            
#     print(f"\nFinal number of bins after merging: {len(merged)}")
#     return merged

def merge_small_bins(bins, min_width):
    """
    Merge bins that are smaller than min_width using a forward-looking approach.
    
    Parameters:
    bins (list): List of [start, end] bin edges
    min_width (float): Minimum allowed bin width
    
    Returns:
    list: New list of [start, end] bin edges after merging
    """
    result = []
    i = 0
    
    while i < len(bins):
        start = bins[i][0]
        end = bins[i][1]
        width = end - start
        
        # If bin is wide enough, add it and continue
        if width >= min_width:
            result.append([start, end])
            i += 1
            continue
        
        # If bin is too small, look ahead and accumulate width until we have enough
        total_width = width
        j = i + 1
        while j < len(bins) and total_width < min_width:
            total_width += bins[j][1] - bins[j][0]
            j += 1
            
        # Create new bins from the accumulated range
        accumulated_range = bins[j-1][1] - start
        n_new_bins = max(1, int(np.floor(accumulated_range / min_width)))
        new_width = accumulated_range / n_new_bins
        
        # Add the new evenly-spaced bins
        for k in range(n_new_bins):
            new_start = start + k * new_width
            new_end = start + (k + 1) * new_width
            result.append([new_start, new_end])
            
        i = j

    # Simple end case: if last bin is too small, merge with previous
    if len(result) > 1 and result[-1][1] - result[-1][0] < min_width:
        last_bin = result.pop()
        result[-1][1] = last_bin[1]

    print("Number of bins merged (reason: size below min_width):", len(bins)-len(result))

    return result


def create_variable_bins(
    values: List[float], 
    mz_range: Tuple[float, float],
    num_bins: int,
    min_width: float = 2.0,
    max_width: float = 50.0,
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

    bins = create_binning_with_constraints(
        values,
        mz_range,
        num_bins,
        max_width
    ) 

    # Final step: merge any bins that ended up too small
    final_bins = merge_small_bins(bins, min_width)
    return final_bins


def get_initial_equal_width_bins(
    mz_range: Tuple[float, float],
    num_splits: int
) -> List[List[float]]:
    """Creates initial bins with equal m/z width by splitting the total range into
    equal segments.

    Parameters:
        mz_range (Tuple[float, float]): (min_mz, max_mz) range to consider
        num_splits (int): Number of desired splits/bins

    Returns:
        List[List[float]]: List of [start, end] ranges for each bin
    """
    # Find the total range to split
    min_val = mz_range[0]
    max_val = mz_range[1]
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


# def create_fixed_bins(
#     values: List[float], 
#     mz_range: Tuple[float, float],
#     num_bins: int,
# ) -> List[List[float]]:
#     """Creates variable bins for mass spectrometry data that balance even distribution
#     of items with reasonable bin widths.

#     Parameters:
#         values (List[float]): List of m/z values to bin
#         mz_range (Tuple[float, float]): (min_mz, max_mz) range to consider
#         num_bins (int): Target number of bins
#         min_width (float): Minimum allowed bin width in Da
#         max_width (float): Maximum allowed bin width in Da

#     Returns:
#         List[List[float]]: List of [start, end] ranges for each bin
#     """
#     # Filter values to specified range and sort them
#     filtered_values = sorted([v for v in values if mz_range[0] <= v <= mz_range[1]])
#     if not filtered_values:
#         return []
    
#     return get_initial_equal_width_bins(
#         mz_range,
#         num_bins
#     )


import numpy as np

def calculate_forbidden_zones(nominal_mass, phospho_enriched):
    """
    Calculate forbidden zones for peptide m/z values based on mass defects.
    
    Parameters:
    -----------
    nominal_mass : float
        The nominal m/z value to calculate the forbidden zone for
    charge_state : int, optional
        Charge state of the peptides (default: 2)
    phospho_enriched : bool, optional
        Whether the sample is phospho-enriched (default: False)
    
    Returns:
    --------
    float
        The m/z value of the forbidden zone
    """
    optimal_mz_increment = 1.00045475
    # Adjust constant based on phospho-enrichment
    optimal_mz_constant = 0.18 if phospho_enriched else 0.25
    
    return (np.ceil(nominal_mass / optimal_mz_increment) * 
            optimal_mz_increment + optimal_mz_constant)

def find_closest_forbidden_zone(value, phospho_enriched):
        # Calculate forbidden zones for value ± 0.5
        lower_forbidden = calculate_forbidden_zones(value - 1, phospho_enriched)
        upper_forbidden = calculate_forbidden_zones(value, phospho_enriched)
        
        # Return the forbidden zone value that's closest to the original value
        return lower_forbidden if abs(value - lower_forbidden) < abs(value - upper_forbidden) else upper_forbidden

def adjust_bin_boundaries(bins, phospho_enriched=False):
    """
    Adjust bin boundaries to be closer to forbidden zones.
    
    Parameters:
    -----------
    bins : list of lists
        List of [start, end] bin boundaries
    phospho_enriched : bool, optional
        Whether the sample is phospho-enriched (default: False)
        
    Returns:
    --------
    list of lists
        Adjusted bin boundaries
    """
    
    adjusted_bins = []
    
    for i in range(len(bins)):
        start_value = bins[i][0]
        end_value = bins[i][1]
        
        # Adjust value (including first and lastbin)
        adjusted_start = find_closest_forbidden_zone(start_value, phospho_enriched)
        adjusted_end = find_closest_forbidden_zone(end_value, phospho_enriched)
            
        adjusted_bins.append([adjusted_start, adjusted_end])
    
    return adjusted_bins


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



