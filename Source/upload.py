import os
import time
import shutil
import common
import streamlit as st
import pandas as pd


type_list_csv = ['csv', 'tsv']
type_list_zip = []
for i in shutil.get_unpack_formats():
    type_list_zip += i[1]

path_uploaded = os.path.join('..', 'Data', 'uploaded')
path_uploaded_genomics = os.path.join(path_uploaded, 'genomics')
path_uploaded_proteomics = os.path.join(path_uploaded, 'proteomics')
path_uploaded_transcriptomics = os.path.join(path_uploaded, 'transcriptomics')
path_uploaded_metabolomics = os.path.join(path_uploaded, 'metabolomics')
path_uploaded_phy_che = os.path.join(path_uploaded, 'phy_che')

path_uploaded_dict = {
    'Genomics': path_uploaded_genomics,
    'Proteomics': path_uploaded_proteomics,
    'Transcriptomics': path_uploaded_transcriptomics,
    'Metabolomics': path_uploaded_metabolomics,
    'Physico-chemical': path_uploaded_phy_che
}


def remove_cached_data():

    SAVE_TIME = 5400  # Time is given in seconds, 5400s = 1.5h
    time_limit = time.time() - SAVE_TIME

    # Traversing every cache directory and deleting everything from it
    for omic in path_uploaded_dict:
        path_omic = path_uploaded_dict[omic]

        for file_name in os.listdir(path_omic):
            file_path = os.path.join(path_omic, file_name)

            # Check if file is older than 1h, and if yes, delete it
            if os.stat(file_path).st_mtime < time_limit:
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except OSError as e:
                    print('Failed to delete ' + file_path + '. Reason: ' + e)

                print('Successfully removed everything from ' + omic +
                      ' data set')

            else:
                pass

    return None


def upload_csv(key_suffix):
    upload_text_csv = {
        'Metabolomics': '''Upload your data set here. Maximum size is 200MB''',
        'Proteomics': '''Upload your data set here. Maximum size is 200MB''',
        'Physico-chemical': '''Upload your data set here. Maximum size is
                               200MB''',
        'Transcriptomics': '''Upload one or more data sets here. Maximum size
                              is 200MB. Data sets must have exactly the same
                              features(columns)'''}

    imported_files = []

    if key_suffix == 'Transcriptomics':
        # This one returns a list because of the accept_multiple_files=True
        imported_files = st.file_uploader(
            upload_text_csv[key_suffix], type=type_list_csv,
            accept_multiple_files=True, key='Upload_file_multi' + key_suffix)

    else:
        # This one returns a file, so it has to be appended
        imported_files.append(st.file_uploader(
            upload_text_csv[key_suffix], type=type_list_csv,
            accept_multiple_files=False,
            key='Upload_file_single' + key_suffix))

    delimiter_dict = {
        'Comma (,)': ',', 'Semicolon (;)': ';', 'Tab (\\t)': '\t'}

    # TODO: Check what happens with semicolon
    default_delimiter_dict = {'csv': 0, 'tsv': 2}

    number_of_files = len(imported_files)

    if imported_files == [] or imported_files[0] is None:
        if key_suffix == 'Transcriptomics':
            return None, None
        else:
            return None

    else:

        # Checking file extension for the first file only
        imported_file_extension = os.path.splitext(
            imported_files[0].name)[1][1:].strip().lower()

        delimiter = st.selectbox(
            'Select the delimiter for your data set',
            list(delimiter_dict.keys()),
            index=default_delimiter_dict[imported_file_extension],
            key='Upload_delim_1_' + key_suffix)

        if number_of_files <= 1:
            warning_text = 'Please choose the right delimiter'
            success_text = 'Data set succesfully uploaded'
        else:
            warning_text = 'Data sets do not have the same delimiter'
            success_text = 'Data sets succesfully uploaded'

        df_list = []
        if key_suffix == 'Transcriptomics':
            selected_df_names = []

        for i in range(number_of_files):
            try:
                # TODO: Implement data imputation, maybe
                df = pd.read_csv(
                    imported_files[i],
                    delimiter=delimiter_dict[delimiter],
                    low_memory=False)
                # This should be used if one wants to generate example 1
                # phy_che data set. First
                # df = common.example_1_fix_double_header(df)
                df.dropna(inplace=True)
                df.reset_index(inplace=True, drop=True)
                df = common.fix_dataframe_columns(df)
            except ValueError:
                st.warning(warning_text)

            df = df.convert_dtypes()
            df_list.append(df)

            if key_suffix == 'Transcriptomics':
                selected_df_names.append(os.path.splitext(
                    imported_files[i].name)[0].strip())

        st.success(success_text)

        if key_suffix == 'Transcriptomics':
            return df_list, selected_df_names

        else:
            return df_list[0]


def upload_multiple(key_suffix):
    # TODO: Check text messages and change where neccessary
    upload_text_zip_fasta = {
        'Genomics': '''Upload your archive here. Archive should
                       contain only FASTA (.fa) files. Possible file names are
                       given as help, on the right.''',
        'Proteomics': '''Upload your archive here. Archive should
                         contain only FASTA (.faa) files. Possible file names
                         are given as help, on the right.'''}

    upload_text_zip_kegg = {
        'Genomics': '''Upload your archive here. Archive should
                       contain only KO besthits (.besthits) files. Possible
                       file names are given as help, on the right.'''}

    upload_help_zip_fasta = {
        'Genomics': '''File names can be given in two formats:
                       1. D03.fa for FASTA file collected on the third day, or
                       W03.fa for FASTA file collected on the third week.
                       You will be given an option to select the start date.
                       2. 2019-03-15.fa for FASTA file collected on 15.03.2019.
                       You should use either the first or the second option,
                       mixing name options is not allowed.''',
        'Proteomics': '''File names can be given in two formats:
                         1. D03.fa[a] for FASTA file collected on the third
                         day, or W03.fa[a] for FASTA file collected the on
                         third week. You will be given an option to select the
                         start date. 2. 2019-03-15.fa[a] for FASTA file
                         collected on 15.03.2019. You should use either the
                         first or the second option, mixing name options is not
                         allowed.'''}

    upload_help_zip_kegg = {
        'Genomics': '''File names can be given in two formats:
                       1. D03.KOs.besthits for FASTA file collected on the
                       third day, or W03.KOs.besthits for FASTA file collected
                       on the third week. You will be given an option to select
                       the start date.
                       2. 2019-03-15.KOs.besthits for FASTA file collected on
                       15.03.2019. You should use either the first or the
                       second option, mixing name options is not allowed.
                       Delimiter in this file should be tab ("\\t").'''}

    upload_text_zip_bins = {
        'Genomics': '''Upload your archive here. Archive should
                       contain only annotation (.gff) files. Possible file
                       names are given as help, on the right.'''}

    upload_help_zip_bins = {
        'Genomics': '''File names can be given in two formats:
                       1. D03.gff for samples collected on the
                       third day, or W03.gff for samples collected
                       on the third week. You will be given an option to select
                       the start date.
                       2. 2019-03-15.gff for samples collected on 15.03.2019.
                       You should use either the first or the second option,
                       mixing name options is not allowed.'''}

    # TODO: Change for transcriptomics and add more type options
    available_data_set_types = {
        'Genomics': {
            'Raw FASTA files': 'FASTA',
            'KEGG annotation files': 'KEGG',
            'BINS annotation files': 'BINS'},
        'Proteomics': {
            'Raw FASTA files': 'FASTA',
            'Calculated data set': 'Calculated'}
    }

    selected_data_set_type = st.selectbox(
        'What kind of data set do you want to upload?',
        list(available_data_set_types[key_suffix].keys())
    )

    if key_suffix == 'Genomics':
        if selected_data_set_type == 'FASTA':

            label_text = upload_text_zip_fasta[key_suffix]
            help_text = upload_help_zip_fasta[key_suffix]

        elif selected_data_set_type == 'KEGG':

            label_text = upload_text_zip_kegg[key_suffix]
            help_text = upload_help_zip_kegg[key_suffix]

        else:
            label_text = upload_text_zip_bins[key_suffix]
            help_text = upload_help_zip_bins[key_suffix]

        imported_file = st.file_uploader(
            label_text, type=type_list_zip, accept_multiple_files=False,
            help=help_text, key='Upload_file_' + key_suffix)

    elif key_suffix == 'Proteomics':
        if selected_data_set_type == 'FASTA':
            imported_file = st.file_uploader(
                upload_text_zip_fasta[key_suffix], type=type_list_zip,
                accept_multiple_files=False,
                help=upload_help_zip_fasta[key_suffix],
                key='Upload_file_' + key_suffix)

        else:
            return (upload_csv(key_suffix),
                    available_data_set_types[key_suffix]
                    [selected_data_set_type])

    else:
        imported_file = None

    if imported_file is not None:
        return (common.import_archive(imported_file,
                                      path_uploaded_dict[key_suffix]),
                available_data_set_types[key_suffix][selected_data_set_type])

    else:
        return None, None


def upload_intro(folder_path, key_suffix):
    st.header(key_suffix + ' data')
    st.markdown('')

    df = None

    if key_suffix in ['Metabolomics', 'Physico-chemical', 'Transcriptomics']:

        # TODO: This should be revised, because there is not such
        # functionality for genomics and proteomics. Maybe it is better
        # that the user uploads the data set each time
        # CALCULATED_DATA_SET_NAME = 'calculated.pkl'
        # CALCULATED_DATA_SET_PATH = os.path.join(
        #     folder_path, CALCULATED_DATA_SET_NAME)

        # if os.path.exists(CALCULATED_DATA_SET_PATH):
        #     df = common.get_cached_dataframe(CALCULATED_DATA_SET_PATH)
        # else:
        #     df = upload_csv(key_suffix)
        #     common.cache_dataframe(df, CALCULATED_DATA_SET_PATH)

        df_list = upload_csv(key_suffix)

        if df_list is None:
            st.warning('Upload your data set')

        return df_list

    else:
        df, data_set_type = upload_multiple(key_suffix)

        if df is None:
            st.warning('Upload your data set')

        return df, data_set_type


def upload_genomics():
    key_suffix = 'Genomics'
    cache_folder_path = path_uploaded_genomics

    folder_path_or_df, data_set_type = upload_intro(
        cache_folder_path, key_suffix)

    return common.work_with_zip(
        folder_path_or_df, data_set_type, cache_folder_path, key_suffix)


def upload_proteomics():

    key_suffix = 'Proteomics'
    cache_folder_path = path_uploaded_proteomics

    folder_path_or_df, data_set_type = upload_intro(
        cache_folder_path, key_suffix)

    return common.work_with_zip(
        folder_path_or_df, data_set_type, cache_folder_path, key_suffix)


def upload_transcriptomics():

    key_suffix = 'Transcriptomics'
    cache_folder_path = path_uploaded_transcriptomics

    df_list, selected_df_names = upload_intro(cache_folder_path, key_suffix)

    df = common.work_with_multi_transcriptomics(
        df_list, selected_df_names)

    return common.work_with_csv(df, cache_folder_path, key_suffix)


def upload_metabolomics():

    key_suffix = 'Metabolomics'
    cache_folder_path = path_uploaded_metabolomics

    df = upload_intro(cache_folder_path, key_suffix)

    return common.work_with_csv(df, cache_folder_path, key_suffix)


def upload_phy_che():

    key_suffix = 'Physico-chemical'
    cache_folder_path = path_uploaded_phy_che

    df = upload_intro(cache_folder_path, key_suffix)

    return common.work_with_csv(df, cache_folder_path, key_suffix)


def create_main_upload():

    st.header('Dataset')

    omics_list = ['Genomics', 'Metabolomics', 'Proteomics',
                  'Physico-chemical', 'Transcriptomics']

    choose_omics = st.multiselect('Which omic data do you want to upload:',
                                  omics_list)

    num_of_columns = len(choose_omics)
    charts = []  # An empty list to hold all pairs (visualizations, key)

    with st.beta_expander('Show/hide data sets and related info',
                          expanded=True):
        if num_of_columns >= 2:
            column_list = st.beta_columns(num_of_columns)
            curr_pos = 0

            for i in choose_omics:
                if i == 'Genomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        charts += upload_genomics()
                elif i == 'Metabolomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        charts += upload_metabolomics()
                elif i == 'Proteomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        charts += upload_proteomics()
                elif i == 'Transcriptomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        charts += upload_transcriptomics()
                else:
                    with column_list[curr_pos]:
                        curr_pos += 1
                        charts += upload_phy_che()

        else:
            for i in choose_omics:
                if i == 'Genomics':
                    charts += upload_genomics()
                elif i == 'Metabolomics':
                    charts += upload_metabolomics()
                elif i == 'Proteomics':
                    charts += upload_proteomics()
                elif i == 'Transcriptomics':
                    charts += upload_transcriptomics()
                else:
                    charts += upload_phy_che()

    with st.beta_expander('Show/hide visualizations', expanded=True):
        for i in charts:
            type_of_chart = type(i[0])

            with st.spinner('Visualizing...'):
                if 'altair' in str(type_of_chart):
                    st.altair_chart(i[0], use_container_width=True)
                else:
                    pass

    return None
