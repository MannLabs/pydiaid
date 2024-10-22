# for data manipulation:
import pandas as pd

def load_MS_method_from_txt_file(
    open_name,
):
    """
    Loads a text file of a mass spectrometry method as a pandas DataFrame,
    after skipping lines that have a '#' symbol or string 'type' and renaming the dataframe 
    column names.

    The function reads the file, identifies and removes any comment lines (lines that
    contain '#' or 'type'), and then loads the remained contents into a dataframe.
    The column names are set for alining the data formatting to the input file.

    Parameters:
    open_name (str): path to the file that is to be loaded.

    Returns:
    pd.DataFrame: returns a data frame with the file data with the original method column names,
    with the comment lines skipped.
    """
    a_file = open(open_name)
    file_contents = a_file.read()
    contents_split = file_contents.splitlines(True)
    indices = [i for i, s in enumerate(contents_split) if '#' in s or "type" in s]

    return pd.read_csv(
                open_name,
                skiprows=indices,
                names=[
                    "type",
                    "mobility pos.1 [1/K0]",
                    "mass pos.1 start [m/z]",
                    "mass pos.1 end [m/z]",
                    "mobility pos.2 [1/K0]",
                    "mass pos.2 start [m/z]"
                ]
            )
