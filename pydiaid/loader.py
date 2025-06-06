# for data manipulation
import pandas as pd
import os
import re

import warnings
warnings.filterwarnings("ignore")


def load_library(
    library_name: str,
    analysis_software: str,
    ptm_list: list,
    require_im: bool = True
) -> pd.DataFrame:
    """Loads an output file of a proteomics analysis software as a data frame,
        which should represent the diversity of possible peptides.

    Parameters:
    library_name (str): path, where the proteomics library is stored.
    analysis_software (str): an identifier for the analysis software used to
        create the input data. The script chooses different parse functions
        depending on this identifier.
    ptm_list (list): a list with identifiers used for filtering a specific data frame column
        for the strings within the list.
    require_im (bool): if True, requires ion mobility data; if False, makes ion mobility optional.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """
    try:
        # Special case for DIANN single-run which needs the file path directly
        if analysis_software == 'DIANN single-run':
            if library_name.split(".")[-1] == "parquet":
                analysis_software = 'DIANN library'
            else:
                return __parse_diann_single_run(library_name, ptm_list, require_im)

        
        # For all other software, load the dataframe first
        dataframe = __load_dataframe_from_file(library_name)
        parsers = {
            'AlphaPept': __parse_alpha_pept,
            'MaxQuant': __parse_max_quant,
            'MSFragger': __parse_ms_fragger,
            'Spectronaut single-run': __parse_spectronaut_single_shot,
            'Spectronaut library': __parse_spectronaut_library,
            'DIANN library': __parse_diann_lib,
            'AlphaPeptDeep library': __parse_alphadeep_lib,
            'OpenSWATH': __parse_openswath,
            'AlphaDIA': __parse_alphadia
        }
        
        if analysis_software not in parsers:
            raise Exception('Analysis software not supported.')
            
        return parsers[analysis_software](dataframe, ptm_list, require_im)
        
    except Exception as e:
        print(e)
        raise Exception("error while processing: Did you choose the correct analysis_software and is the modification present in the dataset?")


def __load_dataframe_from_file(
    library_name: str,
) -> pd.DataFrame:
    """Imports an output file of a proteomics analysis software as a data frame
        independent of the file format.

    Parameters:
    library_name (str): path, where the proteomics library is stored.

    Returns:
    pd.DataFrame: returns a data frame in the same way as the specific
        proteomics analysis software stored it.
    """

    if library_name.split(".")[-1] == "csv":
        return pd.read_csv(library_name, sep=',')
    if library_name.split(".")[-1] == "parquet":
        return pd.read_parquet(library_name, engine='fastparquet')
    else:
        return pd.read_csv(library_name, sep='\t')  # .xls, .tsv, .txt

class ColumnMapper:
    """Handles column name mapping across different software versions."""
    
    # Define all possible column variants as class attributes
    COLUMN_VARIANTS = {
        'decoy': ['decoy', 'Decoy', 'is_decoy'],
        'qvalue': ['QValue', 'Q.Value', 'q_value', 'Q_Value'],
        'mobility': [
            'PrecursorIonMobility',
            'IonMobility',
            'Ion Mobility',
            'ion_mobility',
            'IM',
            'Mobility',
            '1/K0'
        ],
        'mz': ['PrecursorMz', 'Precursor.Mz', 'Mz', 'PrecursorMZ', 'Calibrated Observed M/Z'],
        'charge': ['PrecursorCharge', 'Precursor.Charge', 'Charge'],
        'protein': ['ProteinId', 'ProteinName', 'Protein.Names', 'Protein', 'Protein ID'],
        'modified_peptide': [
            'ModifiedPeptideSequence',
            'ModifiedPeptide',
            'Modified.Sequence',
            'Modified Sequence',
            'Modified Peptide'
        ],
        'peptide': ['Peptide', 'PeptideSequence', 'Sequence']
    }

    def __init__(self, dataframe: pd.DataFrame):
        """Initialize with a dataframe and map its columns."""
        self.df = dataframe
        self.column_map = self._create_column_map()

    def _create_column_map(self) -> dict:
        """Create mapping of standard names to actual column names in dataframe."""
        column_map = {}
        for standard_name, variants in self.COLUMN_VARIANTS.items():
            found_col = next((col for col in variants if col in self.df.columns), None)
            column_map[standard_name] = found_col
        return column_map

    def get_column(self, standard_name: str) -> str:
        """Get the actual column name for a standard column identifier."""
        return self.column_map.get(standard_name)

    def validate_required_columns(self, required_columns: list) -> None:
        """Validate that all required columns exist."""
        missing = [col for col in required_columns if self.get_column(col) is None]
        if missing:
            raise ValueError(f"Required columns missing: {', '.join(missing)}")

    def has_column(self, standard_name: str) -> bool:
        """Check if a standard column exists in the dataframe."""
        return self.get_column(standard_name) is not None


def __parse_alpha_pept(
    dataframe: pd.DataFrame,
    ptm_list: list,
    require_im: bool = True
) -> pd.DataFrame:
    """Filters a data frame depending on the software specific requirements to
        only include valid precursors. Additionally, it parses the data frame
        to library_loader to filter for specific modified peptides and to unify
        the column names of the required columns.

    Parameters:
    dataframe (pd.DataFrame): imported output file from the analysis software
        "AlphaPept". File format: .csv, required columns: "q_value", "decoy",
        'mz', 'mobility', 'charge', 'protein', 'precursor'.
    ptm_list (list): a list with identifiers used for filtering a specific dataframe column.
    require_im (bool): if True, requires ion mobility data; if False, makes ion mobility optional.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """

    filtered_dataframe = dataframe[
        (dataframe["q_value"] <= 0.01) &
        (dataframe["decoy"] == False)
    ]
    
    # Check if mobility column exists
    im_col = 'mobility' if 'mobility' in dataframe.columns else None
    
    if require_im and im_col is None:
        raise Exception("Ion mobility data required but not found in AlphaPept output")
    
    return library_loader(
        filtered_dataframe,
        ptm_list,
        mz='mz',
        im=im_col,
        charge='charge',
        protein='protein',
        modified_peptide='precursor'
    )


def __parse_max_quant(
    dataframe: pd.DataFrame,
    ptm_list: list,
    require_im: bool = True
) -> pd.DataFrame:
    """Filters a data frame depending on the software specific requirements to
        only include valid precursors. Additionally, it parses the data frame
        to library_loader to filter for specific modified peptides and to unify
        the column names of the required columns.

    Parameters:
    dataframe (pd.DataFrame): imported output file from the analysis software
        "MaxQuant".
    File format: evidence.txt, required columns: "Reverse",  'm/z', '1/K0', 'Charge', 
        'Proteins', 'Modified sequence', '1/K0 length'.
    ptm_list (list): a list with identifiers used for filtering a specific dataframe column.
    require_im (bool): if True, requires ion mobility data; if False, makes ion mobility optional.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """

    filtered_dataframe = dataframe[
        dataframe["Reverse"] != "+"
    ]
    
    # Check if IM columns exist
    im_col = '1/K0' if '1/K0' in dataframe.columns else None
    im_length_col = '1/K0 length' if '1/K0 length' in dataframe.columns else None
    
    if require_im and im_col is None:
        raise Exception("Ion mobility data required but not found in MaxQuant output")
    
    return library_loader(
        filtered_dataframe,
        ptm_list,
        mz='m/z',
        im=im_col,
        charge='Charge',
        protein='Proteins',
        modified_peptide='Modified sequence',
        im_length=im_length_col
    )


def __parse_ms_fragger(
    dataframe: pd.DataFrame,
    ptm_list: list,
    require_im: bool = True
) -> pd.DataFrame:
    """It parses the data frame to library_loader to filter for specific
        modified peptides and to unify the column names of the required
        columns.

    Parameters:
    dataframe (pd.DataFrame): imported library or psm file from the analysis software
        "MSFragger". Required columns (supports multiple naming variants)
    File format: .tsv, required columns: 'PrecursorMz', 'PrecursorIonMobility',
        'PrecursorCharge', 'ProteinId', 'ModifiedPeptideSequence'.
    ptm_list (list): a list with identifiers used for filtering a specific dataframe column.
    require_im (bool): if True, requires ion mobility data; if False, makes ion mobility optional.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """

    mapper = ColumnMapper(dataframe)
    
    required_columns = ['mz', 'charge', 'protein', 'modified_peptide']
    if require_im:
        required_columns.append('mobility')
    
    mapper.validate_required_columns(required_columns)
    
    if mapper.has_column('peptide'):
        peptide_col = mapper.get_column('peptide')
        mod_peptide_col = mapper.get_column('modified_peptide')
        dataframe[mod_peptide_col] = dataframe[mod_peptide_col].replace('', pd.NA)
        dataframe[mod_peptide_col] = dataframe[mod_peptide_col].fillna(dataframe[peptide_col])
    
    return library_loader(
        dataframe,
        ptm_list,
        mz=mapper.get_column('mz'),
        im=mapper.get_column('mobility'),
        charge=mapper.get_column('charge'),
        protein=mapper.get_column('protein'),
        modified_peptide=mapper.get_column('modified_peptide')
    )


def __parse_spectronaut_single_shot(
    dataframe: pd.DataFrame,
    ptm_list: list,
    require_im: bool = True
) -> pd.DataFrame:
    """Filters a data frame depending on the software specific requirements to
        only include valid precursors. Additionally, it parses the data frame
        to library_loader to filter for specific modified peptides and to unify
        the column names of the required columns.

    Parameters:
    dataframe (pd.DataFrame): imported output file from the analysis software
        "Spectronaut - DIA analysis".
    File format: .xls, pre-filtered for FDR, required columns: 'PG.ProteinGroups', 'EG.PrecursorId',
        'FG.PrecMzCalibrated', 'FG.ApexIonMobility', 'FG.Charge', 'FG.IonMobilityPeakWidth'.
    ptm_list (list): a list with identifiers used for filtering a specific dataframe column.
    require_im (bool): if True, requires ion mobility data; if False, makes ion mobility optional.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """

    # Basic filtering for valid precursor m/z
    filter_selection = (dataframe['FG.PrecMzCalibrated'].isnull() == False)
    
    # Check if IM columns exist
    im_col = 'FG.ApexIonMobility' if 'FG.ApexIonMobility' in dataframe.columns else None
    im_length_col = 'FG.IonMobilityPeakWidth' if 'FG.IonMobilityPeakWidth' in dataframe.columns else None
    
    # Only apply IM filter if column exists and IM is required
    if im_col is not None:
        filter_selection &= (dataframe[im_col].isnull() == False)
    elif require_im:
        raise Exception("Ion mobility data required but not found in Spectronaut output")
    
    filtered_library = dataframe[filter_selection]
    
    return library_loader(
        filtered_library,
        ptm_list,
        mz='FG.PrecMzCalibrated',
        im=im_col,
        charge='FG.Charge',
        protein='PG.ProteinGroups',
        modified_peptide='EG.PrecursorId',
        im_length=im_length_col
    )


def __parse_spectronaut_library(
    dataframe: pd.DataFrame,
    ptm_list: list,
    require_im: bool = True
) -> pd.DataFrame:
    """Filters a data frame depending on the software specific requirements to
        only include valid precursors. Additionally, it parses the data frame
        to library_loader to filter for specific modified peptides and to unify
        the column names of the required columns.

    Parameters:
    dataframe (pd.DataFrame): imported output file from the analysis software
        "Spectronaut - library generation from Pulsar".
    File format: .xls, required columns: 'PrecursorMz', 'IonMobility',
        'PrecursorCharge', 'UniProtIds', 'ModifiedPeptide'.
    ptm_list (list): a list with identifiers used for filtering a specific dataframe column.
    require_im (bool): if True, requires ion mobility data; if False, makes ion mobility optional.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """
    filtered_dataframe = dataframe[
        dataframe["ExcludeFromAssay"] != True
    ]
    # Basic filtering for valid precursor m/z
    filter_selection = (filtered_dataframe['PrecursorMz'].isnull() == False)
    
    # Check if IM column exists and add to filter if present
    if 'IonMobility' in filtered_dataframe.columns:
        if require_im:
            filter_selection &= (filtered_dataframe['IonMobility'].isnull() == False)
        im_col = 'IonMobility'
    else:
        if require_im:
            raise Exception("Ion mobility data required but not found in Spectronaut library")
        im_col = None
    
    filtered_library = filtered_dataframe[filter_selection]
    
    return library_loader(
        filtered_library,
        ptm_list,
        mz='PrecursorMz',
        im=im_col,
        charge='PrecursorCharge',
        protein='UniProtIds',
        modified_peptide='ModifiedPeptide'
    )


def __parse_diann_lib(
    dataframe: pd.DataFrame,
    ptm_list: list,
    require_im: bool = True
) -> pd.DataFrame:
    """Filters a data frame depending on the software specific requirements to
    only include valid precursors. Additionally, it parses the data frame
    to library_loader to filter for specific modified peptides and to unify
    the column names of the required columns.

    Parameters:
    dataframe (pd.DataFrame): imported library file from the analysis software
        "DIANN". Required columns (supports multiple naming variants)
    ptm_list (list): a list with identifiers used for filtering a specific dataframe column.
    require_im (bool): if True, requires ion mobility data; if False, makes ion mobility optional.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """
        # Initialize column mapper
    mapper = ColumnMapper(dataframe)
    
    # Check required columns
    required_columns = ['decoy', 'qvalue', 'mz', 'charge', 'protein', 'modified_peptide']
    if require_im:
        required_columns.append('mobility')
    
    mapper.validate_required_columns(required_columns)
    
    # Filter dataframe
    filtered_dataframe = dataframe[
        (dataframe[mapper.get_column('decoy')] == 0) &
        (dataframe[mapper.get_column('qvalue')] <= 0.01)
    ]
    
    return library_loader(
        filtered_dataframe,
        ptm_list,
        mz=mapper.get_column('mz'),
        im=mapper.get_column('mobility'),
        charge=mapper.get_column('charge'),
        protein=mapper.get_column('protein'),
        modified_peptide=mapper.get_column('modified_peptide')
    )


from alphabase.psm_reader import psm_reader_provider

def __parse_diann_single_run(
    tsv_file: str,
    ptm_list: list,
    require_im: bool = True
) -> pd.DataFrame:
    """Filters a DIA-NN single-run output file and parses it to unify
    the column names of the required columns.

    Parameters:
    tsv_file (str): path to the DIA-NN single-run output file.
        the file should contain the following columns:
        'precursor_mz': precursor mass-to-charge ratio, calculated by alphabase
        'sequence': peptide sequence, Stripped.Sequence
        'charge': precursor charge state, Precursor.Charge
        'proteins': protein identifiers, Protein.Group
        'fdr': false discovery rate, Q.Value
        Optional columns:
        'mobility': ion mobility, IM
        'mod_sites': modification sites, Modified.Sequence
        'mods': modifications, Modified.Sequence
    ptm_list (list): a list with identifiers used for filtering a specific column
        for modifications.
    require_im (bool): if True, requires ion mobility data; if False, makes ion mobility optional.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
        The returned DataFrame contains filtered and processed data with standardized
        column names for downstream analysis.

    Raises:
    Exception: If ion mobility data is required but not found in the input file.
    """
    # Import the data using psm reader
    dataframe = psm_reader_provider.get_reader('diann').import_file(tsv_file)

    # Filter for high-quality identifications
    filtered_dataframe = dataframe[
        (dataframe["fdr"] <= 0.01)  # Filter for 1% FDR
    ]

    # Check if IM column exists
    im_col = 'mobility' if 'mobility' in dataframe.columns else None
    
    if require_im and im_col is None:
        raise Exception("Ion mobility data required but not found in DIA-NN output")

    # Create modified sequence identifier by combining sequence, mods, and mod_sites
    if 'mods' in dataframe.columns and 'mod_sites' in dataframe.columns:
        filtered_dataframe['modified_sequence'] = filtered_dataframe.apply(
            lambda row: f"{row['sequence']}_{row['mods']}_{row['mod_sites']}" 
            if pd.notna(row['mods']) else row['sequence'], 
            axis=1
        )
    else:
        filtered_dataframe['modified_sequence'] = filtered_dataframe['sequence']

    return library_loader(
        library=filtered_dataframe,
        ptm_list=ptm_list,
        mz='precursor_mz',
        charge='charge',
        protein='proteins',
        modified_peptide='modified_sequence',
        im=im_col
    )


def __parse_alphadeep_lib(
    dataframe: pd.DataFrame,
    ptm_list: list,
    require_im: bool = True
) -> pd.DataFrame:
    """Filters a data frame depending on the software specific requirements to
        only include valid precursors. Additionally, it parses the data frame
        to library_loader to filter for specific modified peptides and to unify
        the column names of the required columns.

    Parameters:
    dataframe (pd.DataFrame): imported output file from the analysis software
        "AlphaPeptDeep". File format: .csv, required columns: 
        'PrecursorMz',
        'IonMobility',
        'PrecursorCharge',
        'ProteinID',
        'ModifiedPeptide'.
    ptm_list (list): a list with identifiers used for filtering a specific dataframe column.
    require_im (bool): if True, requires ion mobility data; if False, makes ion mobility optional.
    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """
    im_col = 'IonMobility' if 'IonMobility' in dataframe.columns else None
    
    if require_im and im_col is None:
        raise Exception("Ion mobility data required but not found in AlphaPeptDeep library")
    
    return library_loader(
        dataframe,
        ptm_list,
        mz='PrecursorMz',
        im=im_col,
        charge='PrecursorCharge',
        protein='ProteinID',
        modified_peptide='ModifiedPeptide'
    )


def __parse_openswath(
    dataframe: pd.DataFrame,
    ptm_list: list,
    require_im: bool = True
) -> pd.DataFrame:
    """Filters a data frame from OpenSWATH output and parses it to unify
    the column names of the required columns.

    Parameters:
    dataframe (pd.DataFrame): imported output file from OpenSWATH. 
        File format: csv/tsv, required columns:
        'PrecursorMz',
        'PrecursorCharge',
        'ProteinName',
        'FullUniModPeptideName',
        Optional columns:
        'IonMobility'
    ptm_list (list): a list with identifiers used for filtering a specific dataframe column.
    require_im (bool): if True, requires ion mobility data; if False, makes ion mobility optional.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """
    # Check if IM column exists
    im_col = 'IonMobility' if 'IonMobility' in dataframe.columns else None
    
    if require_im and im_col is None:
        raise Exception("Ion mobility data required but not found in OpenSWATH data")
    
    # Drop duplicates to get unique precursors
    unique_precursors = dataframe.drop_duplicates(
        ['FullUniModPeptideName', 'PrecursorCharge']
    )
    
    return library_loader(
        library=unique_precursors,
        ptm_list=ptm_list,
        mz='PrecursorMz',
        charge='PrecursorCharge',
        protein='ProteinName',
        modified_peptide='FullUniModPeptideName',
        im=im_col
    )


def __parse_alphadia(
    dataframe: pd.DataFrame,
    ptm_list: list,
    require_im: bool = True
) -> pd.DataFrame:
    """Filters a data frame from AlphaDIA output and parses it to unify
    the column names of the required columns.

    Parameters:
    dataframe (pd.DataFrame): imported output file from AlphaDIA. 
        File format: csv/tsv, required columns:
        'mz_calibrated': calibrated precursor m/z
        'mobility_calibrated': calibrated ion mobility
        'charge': precursor charge state
        'proteins': protein identifiers
        'precursor_idx': precursor identifier
        'base_width_mobility': ion mobility peak width
        'decoy': decoy indicator (0 for targets, 1 for decoys)
        'qval': q-value for false discovery rate control
    ptm_list (list): a list with identifiers used for filtering a specific dataframe column.
    require_im (bool): if True, requires ion mobility data; if False, makes ion mobility optional.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """
    # Check if required columns exist
    required_cols = ['mz_calibrated', 'charge', 'proteins', 'precursor_idx', 'decoy', 'qval']
    missing_cols = [col for col in required_cols if col not in dataframe.columns]
    if missing_cols:
        raise Exception(f"Required columns missing from AlphaDIA output: {missing_cols}")
    
    # Filter for high-quality identifications
    filtered_dataframe = dataframe[
        (dataframe['decoy'] == 0) &  # Keep only target hits
        (dataframe['qval'] <= 0.01)  # Filter at 1% FDR
    ]

    # Check if IM columns exist
    im_col = 'mobility_calibrated' if 'mobility_calibrated' in filtered_dataframe.columns else None
    im_width_col = 'base_width_mobility' if 'base_width_mobility' in filtered_dataframe.columns else None
    
    if require_im and (im_col is None or im_width_col is None):
        raise Exception("Ion mobility data required but not found in AlphaDIA output")

    return library_loader(
        library=filtered_dataframe,
        ptm_list=ptm_list,
        mz='mz_calibrated',
        charge='charge',
        protein='proteins',
        modified_peptide='precursor_idx',
        im=im_col,
        im_length=im_width_col
    )


def library_loader(
    library: pd.DataFrame,
    ptm_list: list,
    mz: str,
    charge: str,
    protein: str,
    modified_peptide: str,
    im: str = None,
    im_length: str = None
) -> pd.DataFrame:
    """Filters a column of a data frame for a specific identifier (e.g.,
        "Phospho") and unifies the required column names.

    Parameters:
    library (pd.DataFrame): data frame, which should be modified. It needs to
        have columns, which indicate the precursor m/z values, precursor ion
        mobility, precursor charge, corresponding proteins, and modified
        peptide sequence (can contain also charge state information)
    ptm_list (list): a list with identifiers used for filtering the data frame 
        column "modified peptide sequence" for a specific modification.
    mz (str): column name of the column containing the precursor m/z value.
    charge (str): column name of the column containing the precursor charge
        state information.
    protein (str): column name of the column stating the corresponding proteins
        to a specific precursor.
    modified_peptide (str): column name containing the modified peptide
        sequence for each precursor (can additionally contain charge state
        information).
    im (str, optional): column name of the column containing the precursor specific ion
        mobility. If None, the IM column will contain NaN values.
    im_length (str, optional): column name containing the ion mobility length information
        per precursor. If None, this column will not be included in the output.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """
    # Filter for PTMs if specified
    if ptm_list != 'None':
        list_library_prefiltered = list()
        for ptm in ptm_list:
            library[ptm] = library[modified_peptide].apply(
                lambda x: find_PTM(x, ptm)
            )
            library_prefiltered = library[library[ptm] == True]
            list_library_prefiltered.append(library_prefiltered)
        library_filtered = pd.concat(list_library_prefiltered, ignore_index=True)
        library_subset = library_filtered.drop_duplicates(
            [modified_peptide, charge]
        )
    else:
        library_subset = library.drop_duplicates(
            [modified_peptide, charge]
        )  # use only unique precursors

    # Create base DataFrame with required columns
    library_small = pd.DataFrame()
    library_small['mz'] = library_subset[mz]
    library_small['Charge'] = library_subset[charge]
    library_small['Proteins'] = library_subset[protein]
    library_small['Peptide'] = library_subset[modified_peptide]
    
    # Add IM column if specified
    if im is not None:
        library_small['IM'] = library_subset[im]
    else:
        library_small['IM'] = pd.NA
    
    # Add IMlength column if specified
    if im_length is not None:
        library_small['IMlength'] = library_subset[im_length]
    
    return library_small



def find_PTM(column_value, ptm):
    """Identifies if a column value indicates a modification.

    Parameters:
    column_value (str): the value of the column.
    PTM (str): name of the item which is searched for, possibilities:
        'Phospho', 'STY' or 'Gly'

    Returns:
    boolean: Indicates if the specific string (ptm) is present in the column
        value, which allows filtering for modified peptides.
    """

    if ptm in column_value:
        return True
    else:
        next

        
def get_file_names_from_directory(
    directory: str,
    extensions_list: list
)-> list:
    """Search for files with the specified extension in the repository and
    return a list of all file names with that extention.

    Parameters
    ----------
    directory : str
        Path to the repository to search in.
    extensions_list : list
        A list of extensions, e.g. ['d', 'hdf'].

    Returns
    -------
    list
        The list of filtered file names based on their extensions.
    """
    file_names = [file for file in os.listdir(directory) if file.split('.')[-1] in extensions_list]
    return file_names


def create_opt_plot_df(
    filename: str
)-> pd.DataFrame:
    """Return the dataframe containing information about the optimized scan
    area coordinates and coverage based on the filename of the .png file.

    Parameters
    ----------
    filename : str
        The name of the .png file showing the kernel density estimation.

    Returns
    -------
    pd.DataFrame
        The data frame contains several columns:
            - parameters: ['A1', 'A2', 'B1', 'B2', 'coverage'];
            - values: showing the values for each parameter extracted from the
            filename.

    """

    values = re.findall(r'(\d+\.\d+)', filename)
    df = pd.DataFrame(
        {
            'parameters': ['A1', 'A2', 'B1', 'B2', 'coverage'],
            'values': values
        }
    )
    return df