import os
import shutil
import common
from tempfile import NamedTemporaryFile
import streamlit as st
import pandas as pd
import datetime as dt
from gensim.models import Word2Vec

path_uploaded = 'uploaded'
path_uploaded_genomics = os.path.join(path_uploaded, 'genomics')
path_uploaded_proteomics = os.path.join(path_uploaded, 'proteomics')
path_uploaded_transcriptomics = os.path.join(path_uploaded, 'transcriptomics')
path_uploaded_metabolomics = os.path.join(path_uploaded, 'metabolomics')
path_uploaded_phy_che = os.path.join(path_uploaded, 'phy_che')
path_ramdisk = os.path.join('dev', 'shm')

# TODO: Implement removing previously extracted archives and files
# it should be done at the end of tool in main
# Take care of everything

type_list_csv = ['csv', 'tsv']
type_list_zip = []
for i in shutil.get_unpack_formats():
    type_list_zip += i[1]


@st.cache(suppress_st_warning=True)
def import_archive(imported_file, key_suffix):

    # Creating the file from BytesIO stream
    tmp_file = NamedTemporaryFile(delete=False, suffix=imported_file.name)
    tmp_file_path = tmp_file.name
    tmp_file.write(imported_file.getvalue())
    tmp_file.flush()
    tmp_file.close()

    index_of_dot = imported_file.name.index('.')
    return_file_name_no_ext = imported_file.name[:index_of_dot]

    try:
        if key_suffix == 'Genomics':
            extract_path = path_uploaded_genomics
        elif key_suffix == 'Transcriptomics':
            extract_path = path_uploaded_transcriptomics
        elif key_suffix == 'Proteomics':
            extract_path = path_uploaded_proteomics
        else:
            st.error('Bad key suffix for archive unpacking')
            return None

        shutil.unpack_archive(tmp_file_path, extract_dir=extract_path)
        return os.path.join(extract_path, return_file_name_no_ext)

    except ValueError:
        st.error('Error while unpacking the archive')
        return None

    finally:
        st.success('Data set succesfully uploaded')
        os.remove(tmp_file_path)


def import_csv(key_suffix):
    upload_text_csv = '''Upload your data set here. Maximum size is 200MB'''

    imported_file = st.file_uploader(
        upload_text_csv, type=type_list_csv, accept_multiple_files=False,
        key='Upload_file_' + key_suffix)

    delimiter_dict = {
        'Comma (,)': ',', 'Semicolon (;)': ';', 'Tab (\\t)': '\t'}

    # TODO: Check what happens with semicolon
    default_delimiter_dict = {'csv': 0, 'tsv': 2}

    if imported_file is not None:
        df = None
        imported_file_extension = os.path.splitext(
            imported_file.name)[1][1:].strip().lower()

        delimiter = st.selectbox(
            'Select the delimiter for your data set',
            list(delimiter_dict.keys()),
            index=default_delimiter_dict[imported_file_extension],
            key='Upload_delim_' + key_suffix)
        try:
            # TODO: Implement data imputation, maybe
            df = pd.read_csv(
                imported_file,
                delimiter=delimiter_dict[delimiter],
                low_memory=False).dropna()
            df.reset_index(inplace=True, drop=True)
            df = common.fix_dataframe_columns(df)
        except ValueError:
            st.warning('Please choose the right delimiter')

        df = df.convert_dtypes()
        st.success('Data set succesfully uploaded')

        return df

    else:
        return None


def import_multiple(key_suffix):
    # TODO: Check text messages and change where neccessary
    upload_text_zip_fasta = {
        'Genomics': '''Upload your archive here. Archive should
                       contain only FASTA (.fa) files. Possible file names are
                       given as help, on the right.''',
        'Proteomics': '''Upload your archive here. Archive should
                         contain only FASTA (.faa) files. Possible file names
                         are given as help, on the right.''',
        'Transcriptomics': '''Upload your archive here. Archive should
                              contain only FASTA (.fa) files. Possible file
                              names are given as help, on the right.'''}

    upload_text_zip_kegg = {
        'Genomics': '''Upload your archive here. Archive should
                       contain only KO besthits (.besthits) files. Possible
                       file names are given as help, on the right.''',
        'Proteomics': '''Upload your archive here. Archive should
                         contain only KO besthits (.besthits) files. Possible
                         file names are given as help, on the right.''',
        'Transcriptomics': '''Upload your archive here. Archive should
                              contain only KO besthits (.besthits) files.
                              Possible file names are given as help, on the
                              right.'''}

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
                         allowed.''',
        'Transcriptomics': '''File names can be given in two formats:
                              1. D03.fa for FASTA file collected on the third
                              day, or W03.fa for FASTA file collected on the
                              third week. You will be given an option to select
                              the start date.
                              2. 2019-03-15.fa for FASTA file collected on
                              15.03.2019. You should use either the first or
                              the second option, mixing name options is not
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
                       Delimiter in this file should be tab ("\\t").''',
        'Proteomics': '''File names can be given in two formats:
                         1. D03.KOs.besthits for FASTA file collected on the
                         third day, or W03.KOs.besthits for FASTA file
                         collected on the third week. You will be given an
                         option to select the start date.
                         2. 2019-03-15.fa for FASTA file collected on
                         15.03.2019. You should use either the first or the
                         second option, mixing name options is not allowed.
                         Delimiter in this file should be tab ("\t").''',
        'Transcriptomics': '''File names can be given in two formats:
                              1. D03.KOs.besthits for FASTA file collected on
                              the third day, or W03.KOs.besthits for FASTA file
                              collected on the third week. You will be given an
                              option to select the start date.
                              2. 2019-03-15.fa for FASTA file collected on
                              15.03.2019. You should use either the first or
                              the second option, mixing name options is not
                              allowed. Delimiter in this file should be tab
                              ("\t").'''}

    # TODO: Change for transcriptomics and add more type options
    available_data_set_types = {
        'Genomics': {
            'Raw FASTA files': 'FASTA',
            'KEGG annotation files': 'KEGG'},
        'Proteomics': {
            'Raw FASTA files': 'FASTA',
            'Calculated data set': 'Calculated'},
        'Transcriptomics': {
            'Raw FASTA files': 'FASTA',
            'KEGG annotation files': 'KEGG'}
    }

    selected_data_set_type = st.selectbox(
        'What kind of data set do you want to upload?',
        list(available_data_set_types[key_suffix].keys())
    )

    if key_suffix == 'Genomics':
        if selected_data_set_type ==\
         list(available_data_set_types[key_suffix].keys())[0]:
            imported_file = st.file_uploader(
                upload_text_zip_fasta[key_suffix], type=type_list_zip,
                accept_multiple_files=False,
                help=upload_help_zip_fasta[key_suffix],
                key='Upload_file_' + key_suffix)

        else:
            imported_file = st.file_uploader(
                upload_text_zip_kegg[key_suffix], type=type_list_zip,
                accept_multiple_files=False,
                help=upload_help_zip_kegg[key_suffix],
                key='Upload_file_' + key_suffix)

    elif key_suffix == 'Proteomics':
        if selected_data_set_type ==\
         list(available_data_set_types[key_suffix].keys())[0]:
            imported_file = st.file_uploader(
                upload_text_zip_fasta[key_suffix], type=type_list_zip,
                accept_multiple_files=False,
                help=upload_help_zip_fasta[key_suffix],
                key='Upload_file_' + key_suffix)

        else:
            return (import_csv(key_suffix),
                    available_data_set_types[key_suffix]
                    [selected_data_set_type])

    # TODO: Deal with this. This is if key_suffix == 'Transcriptomics
    else:
        imported_file = None

    if imported_file is not None:
        return (import_archive(imported_file, key_suffix),
                available_data_set_types[key_suffix][selected_data_set_type])

    else:
        return None, None


def modify_data_set(df, temporal_column, feature_list, key_suffix):

    columns_to_remove = st.multiselect(
        'Select columns to remove', feature_list,
        key='Col_remove_' + key_suffix)

    if len(columns_to_remove) != 0:
        df.drop(columns_to_remove, axis=1, inplace=True)
        feature_list = [i for i in feature_list if i not in columns_to_remove]

    rows_to_remove_text = st.text_input(
        'Insert row numbers to remove, seperated by comma. See help (right)\
         for example.', value='', key='Row_remove_' + key_suffix,
        help='Example: 42 or 2, 3, 15, 55')

    if rows_to_remove_text != '':
        rows_to_remove = [i.strip() for i in rows_to_remove_text.split(',')]
        # First we check if input is good or not
        if any(not row.isnumeric() for row in rows_to_remove):
            st.error('Wrong number input')
            st.stop()

        rows_to_remove = [int(i) for i in rows_to_remove_text.split(',')]
        df.drop(rows_to_remove, axis=0, inplace=True)

    time_to_remove_text = st.text_input(
        'Insert the begining and the end of a time period to keep, seperated\
         by comma. See help (right) for example.',
        value='2011-03-21, 2012-05-03', key='Row_remove_' + key_suffix,
        help='Example: 2011-03-21, 2012-05-03')

    if time_to_remove_text != '':
        try:
            time_to_remove = [dt.datetime.strptime(
                i.strip(), "%Y-%m-%d") for i in time_to_remove_text.split(',')]
            df.drop(df[(df[temporal_column] < time_to_remove[0]) |
                       (df[temporal_column] > time_to_remove[1])].index,
                    inplace=True)
        except ValueError:
            st.error('Wrong date input')
            st.stop()

    df.dropna(inplace=True)
    df.reset_index(inplace=True, drop=True)
    df[feature_list] = df[feature_list].apply(
        pd.to_numeric)
    df = df.convert_dtypes()

    return df, feature_list


# This function changes all file names of 'D' or 'W' type into timestamp type
# This is done in place for the unpacked uploaded directory
def create_zip_temporality(folder_path, file_name_type, key_suffix):

    if file_name_type in ['D', 'W']:
        start_date = st.date_input(
            'Insert start date for your data set:',
            dt.datetime.strptime("2011-03-21", "%Y-%m-%d"),
            key='Date_input_' + key_suffix)

        # This is only needed for example 1 data set
        # TODO: Change this for production data sets
        common.example_1_fix_archive_file_names(start_date, folder_path)
        # common.fix_archive_file_names(start_date, folder_path)
        # After this step, every file has a timestamp as a name

    else:
        # Check if file name starts with a valid timestamp
        tmp_file_name = os.listdir(folder_path)[0].split('.')[0]
        try:
            dt.datetime.strptime(tmp_file_name, '%Y-%m-%d')

        except ValueError:
            st.error(
                '''File names are not valid. File names should start with "D"
                   or "W", or be of a timestamp type (%Y-%m-%d.fa[a])''')
            st.stop()

    return None


@st.cache
def work_with_fasta(data_set_type, folder_path, key_suffix):

    # TODO: Allow user to change number of epochs for training
    MODEL_NAME = 'w2v_model.saved'
    MODEL_PATH = os.path.join(os.path.split(folder_path)[0], MODEL_NAME)
    fasta_files = os.listdir(folder_path)
    num_of_fasta_files = len(fasta_files)

    if os.path.exists(MODEL_PATH):
        w2v_model = Word2Vec.load(MODEL_PATH)

    # If we already created a model, we won't do it again
    else:
        w2v_model, fasta_files = common.import_mags_and_build_model(
            num_of_fasta_files, folder_path)
        w2v_model = common.train_model(
            w2v_model, path_fasta=folder_path, end=num_of_fasta_files)
        w2v_model.save(MODEL_PATH)

    list_of_vectors = common.vectorize_mags(
        w2v_model, path_fasta=folder_path, end=num_of_fasta_files)
    df = pd.DataFrame(list_of_vectors)
    list_of_dates = common.create_temporal_column(
        fasta_files, None, None, 'TIMESTAMP')
    df.insert(0, 'DateTime', list_of_dates)

    return df


@st.cache
def work_calculate_proteomics(data_set_type, folder_path, key_suffix):

    fasta_files = os.listdir(folder_path)
    num_of_fasta_files = len(fasta_files)

    tmp_df = common.import_proteomics(
        path_proteomics=folder_path, end=num_of_fasta_files)

    list_of_dates = common.create_temporal_column(
        fasta_files, None, None, 'TIMESTAMP')
    tmp_df.insert(0, 'DateTime', list_of_dates)

    return tmp_df


def work_with_data_set(df, data_set_type, folder_path, key_suffix):

    chosen_charts = []

    if data_set_type == 'FASTA':
        VECTORIZED_DATA_SET_NAME = 'vectorized.pkl'
        VECTORIZED_DATA_SET_PATH = os.path.join(
            os.path.split(folder_path)[0], VECTORIZED_DATA_SET_NAME)

        if os.path.exists(VECTORIZED_DATA_SET_PATH):
            df = common.get_cached_dataframe(None, VECTORIZED_DATA_SET_PATH)

        else:
            with st.spinner('Vectorizing FASTA files using W2V...'):
                df = work_with_fasta(
                    data_set_type, folder_path, key_suffix)
                common.cache_dataframe(df, None, VECTORIZED_DATA_SET_PATH)

        common.show_calculated_data_set(df, 'Embedded FASTA files')
        labels_list = common.show_clustering_info(df, key_suffix)

        # Traversing pairs in list
        for i in labels_list:
            temporal_feature, feature_list = common.find_temporal_feature(df)
            feature_list = i[0]
            df[i[0]] = i[1]
            chosen_charts += common.visualize_data_set(
                    df, temporal_feature, feature_list,
                    'Cluster_FASTA_' + key_suffix + '_' + i[0])

    elif data_set_type == 'KEGG':
        KEGG_DATA_SET_NAME = 'kegg.pkl'
        KEGG_DATA_SET_PATH = os.path.join(
            os.path.split(folder_path)[0], KEGG_DATA_SET_NAME)

        if os.path.exists(KEGG_DATA_SET_PATH):
            df = common.get_cached_dataframe(None, KEGG_DATA_SET_PATH)

        else:
            with st.spinner('Creating KO matrix...'):
                besthits_files = os.listdir(folder_path)
                num_of_besthits_files = len(besthits_files)
                df = common.import_kegg_and_create_df(
                    end=num_of_besthits_files, path_all_keggs=folder_path)
                common.cache_dataframe(df, None, KEGG_DATA_SET_PATH)

        common.show_calculated_data_set(df, 'Calculated KO matrix')
        labels_list = common.show_clustering_info(df, key_suffix)

        # Traversing pairs in list
        for i in labels_list:
            temporal_feature, feature_list = common.find_temporal_feature(df)
            feature_list = i[0]
            df[i[0]] = i[1]
            chosen_charts += common.visualize_data_set(
                    df, temporal_feature, feature_list,
                    'Cluster_KEGG_' + key_suffix + '_' + i[0])

    # TODO: This elif was never used and has to be checked ASAP
    elif data_set_type == 'Calculated':
        CALCULATED_DATA_SET_NAME = 'calculated.pkl'
        CALCULATED_DATA_SET_PATH = os.path.join(
            os.path.split(folder_path)[0], CALCULATED_DATA_SET_NAME)

        df = common.fix_data_set(df)
        temporal_feature, feature_list = common.find_temporal_feature(df)
        df, feature_list = modify_data_set(
            df, temporal_feature, feature_list, key_suffix)

        if not os.path.exists(CALCULATED_DATA_SET_PATH):
            common.cache_dataframe(df, None, CALCULATED_DATA_SET_PATH)

        chosen_charts = common.visualize_data_set(
            df, temporal_feature, feature_list, key_suffix)

    elif data_set_type == 'Calculate_now':
        CALCULATED_NOW_DATA_SET_NAME = 'calculated_now.pkl'
        CALCULATED_NOW_DATA_SET_PATH = os.path.join(
            os.path.split(folder_path)[0], CALCULATED_NOW_DATA_SET_NAME)

        if os.path.exists(CALCULATED_NOW_DATA_SET_PATH):
            df = common.get_cached_dataframe(
                None, CALCULATED_NOW_DATA_SET_PATH)

        else:
            with st.spinner('Calculating additional properties...'):
                df = work_calculate_proteomics(
                    data_set_type, folder_path, key_suffix)
                common.cache_dataframe(df, None, CALCULATED_NOW_DATA_SET_PATH)

                common.show_calculated_data_set(df, 'Additional properties')

        df = common.fix_data_set(df)
        temporal_feature, feature_list = common.find_temporal_feature(df)
        df, feature_list = modify_data_set(
            df, temporal_feature, feature_list, key_suffix)
        chosen_charts = common.visualize_data_set(
            df, temporal_feature, feature_list, key_suffix)

    else:
        pass

    return chosen_charts


def upload_intro(folder_path, key_suffix):
    st.header(key_suffix)
    st.markdown('')

    df = None

    if key_suffix in ['Metabolomics', 'Physico-chemical']:

        CALCULATED_DATA_SET_NAME = 'calculated.pkl'
        CALCULATED_DATA_SET_PATH = os.path.join(
            folder_path, CALCULATED_DATA_SET_NAME)

        if os.path.exists(CALCULATED_DATA_SET_PATH):
            df = common.cache_dataframe(None, CALCULATED_DATA_SET_PATH)
        else:
            df = import_csv(key_suffix)

        if df is None:
            st.warning('Upload your data set')
            st.stop()

        return df

    else:
        df, data_set_type = import_multiple(key_suffix)

        if df is None:
            st.warning('Upload your data set')
            st.stop()

        return df, data_set_type


def upload_genomics():
    # TODO: Find out if there is a possibility to calculate physico-chemical
    # properties for genes as it is done with proteins
    key_suffix = 'Genomics'

    folder_path_or_df, data_set_type = upload_intro(
        path_uploaded_genomics, key_suffix)

    if data_set_type == 'FASTA':
        file_name_type = common.show_folder_structure(folder_path_or_df)
        create_zip_temporality(folder_path_or_df, file_name_type, key_suffix)

        chosen_charts = work_with_data_set(
            None, data_set_type, folder_path_or_df, key_suffix)

    elif data_set_type == 'KEGG':
        file_name_type = common.show_folder_structure(folder_path_or_df)
        create_zip_temporality(folder_path_or_df, file_name_type, key_suffix)

        chosen_charts = work_with_data_set(
            None, data_set_type, folder_path_or_df, key_suffix)

    else:
        common.show_data_set(folder_path_or_df)
        chosen_charts = work_with_data_set(
            folder_path_or_df, 'Calculated', path_uploaded_genomics,
            key_suffix)

    return chosen_charts


def upload_proteomics():

    key_suffix = 'Proteomics'

    folder_path_or_df, data_set_type = upload_intro(
        path_uploaded_proteomics, key_suffix)

    # TODO: If we have proteomics FASTA files, we can calculate different
    # physico-chemical properties and show them as well
    if data_set_type == 'FASTA':
        file_name_type = common.show_folder_structure(folder_path_or_df)
        create_zip_temporality(folder_path_or_df, file_name_type, key_suffix)

        # Calculating additional physico-chemical properties
        chosen_charts = []
        additional_check = st.checkbox(
            'Calculate additional physico-chemical properties?',
            value=False, key='Additional_check_' + key_suffix)

        if additional_check:
            chosen_charts = work_with_data_set(
                None, 'Calculate_now', folder_path_or_df, key_suffix)

        chosen_charts += work_with_data_set(
            None, data_set_type, folder_path_or_df, key_suffix)

    else:
        common.show_data_set(folder_path_or_df)
        chosen_charts = work_with_data_set(
            folder_path_or_df, 'Calculated', path_uploaded_proteomics,
            key_suffix)

    return chosen_charts


def upload_transcriptomics():

    key_suffix = 'Transcriptomics'

    folder_path_or_df, data_set_type = upload_intro(
        path_uploaded_transcriptomics, key_suffix)

    if data_set_type == 'FASTA':
        file_name_type = common.show_folder_structure(folder_path_or_df)
        create_zip_temporality(folder_path_or_df, file_name_type, key_suffix)

        chosen_charts = work_with_data_set(
            None, data_set_type, folder_path_or_df, key_suffix)

    else:
        common.show_data_set(folder_path_or_df)
        chosen_charts = work_with_data_set(
            folder_path_or_df, 'Calculated', path_uploaded_transcriptomics,
            key_suffix)

    return chosen_charts


def upload_metabolomics():

    key_suffix = 'Metabolomics'

    df = upload_intro(path_uploaded_metabolomics, key_suffix)
    common.show_data_set(df)
    chosen_charts = work_with_data_set(
        df, 'Calculated', path_uploaded_metabolomics, key_suffix)

    return chosen_charts


def upload_phy_che():

    key_suffix = 'Physico-chemical'

    df = upload_intro(path_uploaded_phy_che, key_suffix)
    common.show_data_set(df)
    chosen_charts = work_with_data_set(
        df, 'Calculated', path_uploaded_phy_che, key_suffix)

    return chosen_charts


def create_main_upload():

    st.header('Dataset')

    upload_omics_list = ['Genomics', 'Metabolomics', 'Proteomics',
                         'Physico-chemical', 'Transcriptomics']

    choose_omics = st.multiselect('Which omic data do you want to upload:',
                                  upload_omics_list)

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
