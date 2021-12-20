import os
import streamlit as st
import pandas as pd
import common


__author__ = 'Aleksandar Anžel'
__copyright__ = ''
__credits__ = ['Aleksandar Anžel', 'Georges Hattab']
__license__ = 'GNU General Public License v3.0'
__version__ = '1.0'
__maintainer__ = 'Aleksandar Anžel'
__email__ = 'aleksandar.anzel@uni-marburg.de'
__status__ = 'Dev'


path_example_1_root_data = os.path.join('..', 'Data', 'cached', 'example_1')
path_example_1_genomics = os.path.join(path_example_1_root_data, 'genomics')
path_example_1_proteomics = os.path.join(
    path_example_1_root_data, 'proteomics')
path_example_1_transcriptomics = os.path.join(
    path_example_1_root_data, 'transcriptomics')
path_example_1_metabolomics = os.path.join(
    path_example_1_root_data, 'metabolomics')
path_example_1_phy_che = os.path.join(
    path_example_1_root_data, 'phy_che')
path_example_1_viz = os.path.join(
    path_example_1_root_data, 'visualizations')

path_example_1_genomics_fasta = os.path.join(
    path_example_1_genomics, 'rmags_filtered')
path_example_1_genomics_kegg = os.path.join(path_example_1_genomics, 'KEGG')
path_example_1_genomics_bins = os.path.join(path_example_1_genomics, 'Bins')
path_example_1_genomics_depths = os.path.join(
    path_example_1_genomics, 'MG_Depths')
path_example_1_transcriptomics_depths = os.path.join(
    path_example_1_transcriptomics, 'MT_Depths')
path_example_1_proteomics_fasta = os.path.join(
    path_example_1_proteomics, 'set_of_78')


def upload_multiple(key_suffix):

    # TODO: Add KEGG data set, since it is not yet ready for production
    # 'KEGG annotation files': 'KEGG',
    available_data_set_types = {
        'Metagenomics': {
            'Raw FASTA files': 'FASTA',
            'BINS annotation files': 'BINS',
            'Depth-of-coverage': 'DEPTH'},
        'Metaproteomics': {
            'Raw FASTA files': 'FASTA'},
        'Metatranscriptomics': {
            'Depth-of-coverage': 'DEPTH'}
    }

    selected_data_set_type = st.selectbox(
        'What kind of data set do you want to see?',
        list(available_data_set_types[key_suffix].keys())
    )

    if key_suffix == 'Metagenomics':
        if selected_data_set_type == 'Raw FASTA files':

            return_path = path_example_1_genomics_fasta

        # elif selected_data_set_type == 'KEGG annotation files':
        #     return_path = path_example_1_genomics_kegg
        elif selected_data_set_type == 'Depth-of-coverage':
            return_path = path_example_1_genomics_depths

        else:
            return_path = path_example_1_genomics_bins

    elif key_suffix == 'Metaproteomics':
        return_path = path_example_1_proteomics_fasta

    elif key_suffix == 'Metatranscriptomics':
        return_path = path_example_1_transcriptomics_depths

    else:
        pass

    return (return_path,
            available_data_set_types[key_suffix][selected_data_set_type])


def upload_intro(folder_path, key_suffix):
    st.header(key_suffix + ' data')
    st.markdown('')

    df = None

    if key_suffix in ['Metabolomics', 'Physico-chemical']:

        CALCULATED_DATA_SET_NAME = 'calculated.pkl'
        CALCULATED_DATA_SET_PATH = os.path.join(
            folder_path, CALCULATED_DATA_SET_NAME)

        if os.path.exists(CALCULATED_DATA_SET_PATH):
            df = common.get_cached_dataframe(CALCULATED_DATA_SET_PATH)
        else:
            st.error('Wrong cache path')
            st.stop()

        return df

    else:
        df, data_set_type = upload_multiple(key_suffix)

        if df is None:
            st.warning('Upload your data set')

        return df, data_set_type


def example_1_genomics():
    key_suffix = 'Metagenomics'
    cache_folder_path = path_example_1_genomics

    folder_path_or_df, data_set_type = upload_intro(
        cache_folder_path, key_suffix)

    return common.work_with_zip(
        folder_path_or_df, data_set_type, cache_folder_path, key_suffix)


def example_1_proteomics():
    key_suffix = 'Metaproteomics'
    cache_folder_path = path_example_1_proteomics

    folder_path_or_df, data_set_type = upload_intro(
        cache_folder_path, key_suffix)

    return common.work_with_zip(
        folder_path_or_df, data_set_type, cache_folder_path, key_suffix)


def example_1_metabolomics():
    key_suffix = 'Metabolomics'
    cache_folder_path = path_example_1_metabolomics

    df = upload_intro(cache_folder_path, key_suffix)

    return common.work_with_csv(df, cache_folder_path, key_suffix)


def example_1_transcriptomics():
    key_suffix = 'Metatranscriptomics'
    cache_folder_path = path_example_1_transcriptomics

    folder_path_or_df, data_set_type = upload_intro(
        cache_folder_path, key_suffix)

    return common.work_with_zip(
        folder_path_or_df, data_set_type, cache_folder_path, key_suffix)


def example_1_phy_che():
    key_suffix = 'Physico-chemical'
    cache_folder_path = path_example_1_phy_che

    df = upload_intro(cache_folder_path, key_suffix)

    return common.work_with_csv(df, cache_folder_path, key_suffix)


def create_main_example_1():

    col_1, col_2 = st.columns([1, 2])

    col_1.info('''
            This data set comes from the following paper:

            **Herold, M., Martínez Arbas, S., Narayanasamy, S. et al.\
            Integration of time-series meta-omics data reveals how microbial\
            ecosystems respond to disturbance. Nat Commun 11, 5281(2020).\
            https://doi.org/10.1038/s41467-020-19006-2**. Analyzed samples\
            were collected from a biological wastewater treatment plant in\
            Schifflange, Luxembourg (49.513414, 6.017925). A precise location
            is shown on the map located on the right.

            It contains **genomics**, **metabolomics**, **proteomics**, and\
            **physico-chemical** data. The code used to parse the data can be\
            found here: [GitLab]
            (https://git-r3lab.uni.lu/malte.herold/laots_niche_ecology_analysis)
            ''')

    col_2.map(pd.DataFrame({'lat': [49.513414], 'lon': [6.017925]}),
              zoom=8, use_container_width=True)

    example_1_omics_list = ['Metagenomics', 'Metabolomics', 'Metaproteomics',
                            'Metatranscriptomics', 'Physico-chemical']
    choose_omics = st.multiselect(
        'What kind of data set do you want to see?', example_1_omics_list)

    num_of_columns = len(choose_omics)

    charts = []  # An empty list to hold all pairs (visualizations, key)

    if num_of_columns >= 2:
        column_list = st.columns(num_of_columns)
        curr_pos = 0

        for i in choose_omics:
            if i == 'Metagenomics':
                with column_list[curr_pos]:
                    curr_pos += 1
                    charts += example_1_genomics()
            elif i == 'Metabolomics':
                with column_list[curr_pos]:
                    curr_pos += 1
                    charts += example_1_metabolomics()
            elif i == 'Metaproteomics':
                with column_list[curr_pos]:
                    curr_pos += 1
                    charts += example_1_proteomics()
            elif i == 'Metatranscriptomics':
                with column_list[curr_pos]:
                    curr_pos += 1
                    charts += example_1_transcriptomics()
            elif i == 'Physico-chemical':
                with column_list[curr_pos]:
                    curr_pos += 1
                    charts += example_1_phy_che()
            else:
                pass

    else:
        for i in choose_omics:
            if i == 'Metagenomics':
                charts += example_1_genomics()
            elif i == 'Metabolomics':
                charts += example_1_metabolomics()
            elif i == 'Metaproteomics':
                charts += example_1_proteomics()
            elif i == 'Metatranscriptomics':
                charts += example_1_transcriptomics()
            elif i == 'Physico-chemical':
                charts += example_1_phy_che()
            else:
                pass

    for i in charts:
        type_of_chart = type(i[0])

        with st.spinner('Visualizing...'):
            if 'altair' in str(type_of_chart):
                st.altair_chart(i[0], use_container_width=True)
                common.save_chart(i[0], path_example_1_viz, i[1])
            else:
                pass

    return None
