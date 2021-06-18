import os
import shutil
import common
from tempfile import NamedTemporaryFile
import streamlit as st
import pandas as pd


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

type_list_csv = ['csv']
type_list_zip = []
for i in shutil.get_unpack_formats():
    type_list_zip += i[1]


def show_data_set(df):
    with st.beta_expander('Show the data set and related info', expanded=True):
        st.markdown('First 100 entries')
        st.dataframe(df.head(100))
        st.dataframe(df.describe())

    return None


def import_archive(imported_file, key_suffix):

    # Creating the file from BytesIO stream
    tmp_file = NamedTemporaryFile(delete=False, suffix=imported_file.name)
    tmp_file_path = tmp_file.name
    tmp_file.write(imported_file.getvalue())
    tmp_file.flush()
    tmp_file.close()

    try:
        if key_suffix == 'genomics':
            shutil.unpack_archive(
                tmp_file_path, extract_dir=path_uploaded_genomics)
        elif key_suffix == 'transcriptomics':
            shutil.unpack_archive(
                tmp_file_path, extract_dir=path_uploaded_transcriptomics)
        elif key_suffix == 'proteomics':
            shutil.unpack_archive(
                tmp_file_path, extract_dir=path_uploaded_proteomics)
        else:
            pass

    except ValueError:
        st.error('Error while unpacking the archive')
        return None

    finally:
        os.remove(tmp_file_path)


def upload_data_set(file_types, key_suffix):

    upload_text_csv = '''Upload your data set here. Maximum size is 200MB'''

    # TODO: Change the text for transcriptomics
    upload_text_zip_genomics = '''Upload your archive here. Archive should
                                  contain only FASTA (.fa) files named:
                                  D[day_number].fa\nFor example: D01.fa'''
    upload_text_zip_proteomics = '''Upload your archive here. Archive should
                                    contain only FASTA (.faa) files named:
                                    D[day_number].fa\nFor example: D01.faa'''
    upload_text_zip_transcriptomics = '''Upload your archive here. Archive
                                         should contain only FASTA (.fa) files
                                         named: D[day_number].fa\nFor example:
                                         D01.fa'''

    if file_types == type_list_csv:
        imported_file = st.file_uploader(upload_text_csv, type=file_types,
                                         accept_multiple_files=False,
                                         key='Upload_file_' + key_suffix)

        delimiter_dict = {
            'Comma (,)': ',', 'Semicolon (;)': ';', 'Tab (\\t)': '\t'}

        delimiter = st.selectbox('Select the delimiter in your data set',
                                 list(delimiter_dict.keys()),
                                 key='Upload_delim_' + key_suffix)

        if imported_file is not None:
            df = None
            try:
                df = pd.read_csv(
                    imported_file, delimiter=delimiter_dict[delimiter])
            except ValueError:
                st.warning('Please choose the right delimiter')

            return df

        else:
            return None

    else:
        if key_suffix == 'genomics':
            imported_file = st.file_uploader(
                upload_text_zip_genomics, type=file_types,
                accept_multiple_files=False, key='Upload_file_' + key_suffix)
        elif key_suffix == 'proteomics':
            imported_file = st.file_uploader(
                upload_text_zip_proteomics, type=file_types,
                accept_multiple_files=False, key='Upload_file_' + key_suffix)
        else:
            imported_file = st.file_uploader(
                upload_text_zip_transcriptomics, type=file_types,
                accept_multiple_files=False, key='Upload_file_' + key_suffix)

        if imported_file is not None:
            return import_archive(imported_file, key_suffix)

        else:
            return None


def upload_genomics():
    st.header('Genomics')
    st.markdown('')

    # TODO: Deal with other types of genomics data, like KEGG etc.
    df = upload_data_set(type_list_zip, 'genomics')

    if df is None:
        st.warning('Please upload your data set')
        return []

    st.success('Data set succesfully uploaded')
    show_data_set(df)

    return []


def upload_proteomics():
    st.header('Proteomics')
    st.markdown('')

    df = upload_data_set(type_list_zip, 'proteomics')

    if df is None:
        st.warning('Upload your data set')
        return []

    st.success('Data set succesfully uploaded')
    show_data_set(df)

    return []


def upload_metabolomics():
    st.header('Metabolomics')
    st.markdown('')

    df = upload_data_set(type_list_csv, 'metabolomics')

    if df is None:
        st.warning('Upload your data set')
        return []

    st.success('Data set succesfully uploaded')
    show_data_set(df)
    temporal_feature, feature_list = common.find_temporal_feature(df)

    chosen_charts = common.visualize_data_set(
        df, temporal_feature, feature_list, 'Metabolomics')

    return chosen_charts


def upload_transcriptomics():
    st.header('Transcriptomics')
    st.markdown('')

    df = upload_data_set(type_list_zip, 'transcriptomics')

    if df is None:
        st.warning('Upload your data set')
        return []

    st.success('Data set succesfully uploaded')
    show_data_set(df)

    return []


def upload_phy_che():
    st.header('Physico-chemical')
    st.markdown('')

    df = upload_data_set(type_list_csv, 'phy_che')

    if df is None:
        st.warning('Upload your data set')
        return []

    st.success('Data set succesfully uploaded')
    show_data_set(df)

    return []


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

    # TODO: Uncomment when you deal with all upload_* functions
    # with st.beta_expander('Show/hide visualizations', expanded=True):
    #     for i in charts:
    #         type_of_chart = type(i[0])

    #         with st.spinner('Visualizing...'):
    #             if 'altair' in str(type_of_chart):
    #                 st.altair_chart(i[0], use_container_width=True)
    #             else:
    #                 pass

    return None
