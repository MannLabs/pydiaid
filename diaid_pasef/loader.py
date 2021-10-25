# for data manipulation:
import pandas as pd


def load_library(
    library_name: str,
    analysis_software: str,
    ptm: str,
) -> pd.DataFrame:
    """Loads an output file of a proteomics analysis software as a data frame,
        which should represent the diversity of possible peptides.

    Parameters:
    library_name (str): path where the proteomics library is stored.
    analysis_software (str): an identifier for the analysis software used to
        create the input data. The script chooses different parse functions
        depending on this identifier.
    ptm (str): an identifier used for filtering a specific data frame column
        for this string.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names
        independent of the input data frame analyzed with individual analysis
        software.
    """

    try:
        dataframe = __load_dataframe_from_file(library_name)
        if analysis_software == 'AlphaPept':
            return __parse_alpha_pept(dataframe, ptm)
        if analysis_software == 'MaxQuant':
            return __parse_max_quant(dataframe, ptm)
        if analysis_software == 'MS_Fragger':
            return __parse_ms_fragger(dataframe, ptm)
        if analysis_software == 'Spectronaut_single-shot':
            return __parse_spectronaut_single_shot(dataframe, ptm)
        if analysis_software == 'Spectronaut_library':
            return __parse_spectronaut_library(dataframe, ptm)
        raise Exception('Analysis software not supported.')
    except Exception as e:
        print(e)
        raise Exception("error while processing: Did you choose the correct analysis_software?")


def __load_dataframe_from_file(
    library_name: str,
) -> pd.DataFrame:
    """Imports an output file of a proteomics analysis software as a data frame
        independent of the file format.

    Parameters:
    library_name (str): path where the proteomics library is stored.

    Returns:
    pd.DataFrame: returns a data frame in the same way as the specific
        proteomics analysis software stored it.
    """

    if library_name.split(".")[-1] == "csv":
        return pd.read_csv(library_name, sep=',')
    else:
        return pd.read_csv(library_name, sep='\t')  # .xls, .tsv, .txt


def __parse_alpha_pept(
    dataframe: pd.DataFrame,
    ptm: str,
) -> pd.DataFrame:
    """Filters a data frame depending on the software specific requirements to
        only include valid precursors. Additionally, it parses the data frame
        to library_loader to filter for specific modified peptides and to unify
        the column names of the required columns.

    Parameters:
    dataframe (pd.DataFrame): imported output file from the analysis software
        "AlphaPept". File format: .csv, required columns: "q_value", "decoy",
        'mz', 'mobility', 'charge', 'protein', 'precursor'.
    ptm (str): an identifier used for filtering a specific data frame column.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """

    filtered_dataframe = dataframe[
        (dataframe["q_value"] <= 0.01) &
        (dataframe["decoy"] is False)
    ]
    library_subset = library_loader(
        filtered_dataframe,
        ptm,
        'mz',
        'mobility',
        'charge',
        'protein',
        'precursor'
    )
    return library_subset


def __parse_max_quant(
    dataframe: pd.DataFrame,
    ptm: str,
) -> pd.DataFrame:
    """Filters a data frame depending on the software specific requirements to
        only include valid precursors. Additionally, it parses the data frame
        to library_loader to filter for specific modified peptides and to unify
        the column names of the required columns.

    Parameters:
    dataframe (pd.DataFrame): imported output file from the analysis software
        "MaxQuant".
    File format: evidence.txt, required columns: "Reverse", "Potential
        contaminant", 'm/z', '1/K0', 'Charge', 'Proteins', 'Modified sequence'.
    ptm (str): an identifier used for filtering a specific data frame column.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """

    filtered_dataframe = dataframe[
        (dataframe["Reverse"] != "+") &
        (dataframe["Potential contaminant"] != "+")
    ]
    library_subset = library_loader(
        filtered_dataframe,
        ptm,
        'm/z',
        '1/K0',
        'Charge',
        'Proteins',
        'Modified sequence'
    )
    return library_subset

# todo: Maria told me about new MSFragger versions, which might have a
# different file layout.


def __parse_ms_fragger(
    dataframe: pd.DataFrame,
    ptm: str,
) -> pd.DataFrame:
    """It parses the data frame to library_loader to filter for specific
        modified peptides and to unify the column names of the required
        columns.

    Parameters:
    dataframe (pd.DataFrame): imported output file from the analysis software
        "MSFragger".
    File format: .tsv, required columns: 'PrecursorMz', 'PrecursorIonMobility',
        'PrecursorCharge', 'ProteinId', 'ModifiedPeptideSequence'.
    ptm (str): an identifier used for filtering a specific data frame column.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """
    #dataframe["ModifiedPeptideSequence"] = dataframe.apply(
    #    lambda x: combine_columns(x, dataframe),
    #    axis=1
    #)

    library_subset = library_loader(
        dataframe,
        ptm,
        #'Calibrated Observed M/Z',
        #'Ion Mobility',
        #'Charge',
        #'Protein ID',
        #"ModifiedPeptideSequence"
         'PrecursorMz',
         'PrecursorIonMobility',
         'PrecursorCharge',
         'ProteinId',
         'ModifiedPeptideSequence'
    )
    return library_subset


def combine_columns(x, dataframe):
    if pd.isna(x['Modified Peptide']) is True:
        return x['Peptide']
    else:
        return x['Modified Peptide']


def __parse_spectronaut_single_shot(
    dataframe: pd.DataFrame,
    ptm: str,
) -> pd.DataFrame:
    """Filters a data frame depending on the software specific requirements to
        only include valid precursors. Additionally, it parses the data frame
        to library_loader to filter for specific modified peptides and to unify
        the column names of the required columns.

    Parameters:
    dataframe (pd.DataFrame): imported output file from the analysis software
        "Spectronaut - DIA analysis".
    File format: .xls, required columns: 'PG.ProteinGroups', 'EG.PrecursorId',
        'FG.PrecMzCalibrated', 'FG.ApexIonMobility', 'FG.Charge'.
    ptm (str): an identifier used for filtering a specific data frame column.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """

    filter_selection = (dataframe['FG.PrecMzCalibrated'].isnull() == False)
    filter_selection &= (dataframe['FG.ApexIonMobility'].isnull() == False)
    filtered_library = dataframe[filter_selection]
    library_subset = library_loader(
        filtered_library,
        ptm,
        'FG.PrecMzCalibrated',
        'FG.ApexIonMobility',
        'FG.Charge',
        'PG.ProteinGroups',
        'EG.PrecursorId'
    )
    return library_subset


def __parse_spectronaut_library(
    dataframe: pd.DataFrame,
    ptm: str,
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
    ptm (str): an identifier used for filtering a specific data frame column.

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """
    filter_selection = (dataframe['PrecursorMz'].isnull() == False)
    filter_selection &= (dataframe['IonMobility'].isnull() == False)
    filtered_library = dataframe[filter_selection]

    library_subset = library_loader(
        filtered_library,
        ptm,
        'PrecursorMz',
        'IonMobility',
        'PrecursorCharge',
        'UniProtIds',
        'ModifiedPeptide'
    )

    return library_subset


def library_loader(
    library: pd.DataFrame,
    ptm: str,
    mz: str,
    im: str,
    charge: str,
    protein: str,
    modified_peptide: str,
) -> pd.DataFrame:
    """Filters a column of a data frame for a specific identifier (e.g.,
        "Phospho") and unifies the required column names.

    Parameters:
    library (pd.DataFrame): data frame, which should be modified. It needs to
        have columns, which indicate the precursor m/z values, precursor ion
        mobility, precursor charge, corresponding proteins, and modified
        peptide sequence (can contain also charge state information)
    ptm (str): an identifier used for filtering the data frame column "modified
        peptide sequence" for a specific modification.
    mz (str): column name of the column containing the precursor m/z value.
    im (str): column name of the column containing the precursor specific ion
        mobility.
    charge (str): column name of the column containing the precursor charge
        state information.
    protein (str): column name of the column stating the corresponding proteins
        to a specific precursor.
    modified_peptide (str):  column name containing the modified peptide
        sequence for each precursor (can additionally contain charge state
        information).

    Returns:
    pd.DataFrame: returns a pre-filtered data frame with unified column names.
    """

    if ptm != 'None':
        library[ptm] = library[modified_peptide].apply(
            lambda x: find_PTM(x, ptm)
        )
        library_filtered = library[library[ptm] is True]
        library_subset = library_filtered.drop_duplicates(
            [modified_peptide, charge]
        )
    if ptm == 'None':
        library_subset = library.drop_duplicates(
            [modified_peptide, charge]
        )  # use only unique precursors

    library_small = pd.DataFrame()
    library_small['mz'] = library_subset[mz]
    library_small['IM'] = library_subset[im]
    library_small['Charge'] = library_subset[charge]
    library_small['Proteins'] = library_subset[protein]
    library_small['Peptide'] = library_subset[modified_peptide]

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
