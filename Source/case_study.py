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


path_case_study_root_data = os.path.join('..', 'Data', 'cached', 'example_1')
path_case_study_genomics = os.path.join(path_case_study_root_data, 'genomics')
path_case_study_proteomics = os.path.join(
    path_case_study_root_data, 'proteomics')
path_case_study_transcriptomics = os.path.join(
    path_case_study_root_data, 'transcriptomics')
path_case_study_metabolomics = os.path.join(
    path_case_study_root_data, 'metabolomics')
path_case_study_phy_che = os.path.join(
    path_case_study_root_data, 'phy_che')
path_case_study_viz = os.path.join(
    path_case_study_root_data, 'visualizations')

CALCULATED_DATA_SET_NAME = 'calculated.pkl'
CALCULATED_NOW_DATA_SET_NAME = 'calculated_now.pkl'
path_case_study_genomics_fasta = os.path.join(
    path_case_study_genomics, 'rmags_filtered')
path_case_study_genomics_kegg = os.path.join(path_case_study_genomics, 'KEGG')
path_case_study_genomics_bins = os.path.join(path_case_study_genomics, 'Bins')
path_case_study_genomics_depths = os.path.join(
    path_case_study_genomics, 'MG_Depths')
path_case_study_transcriptomics_depths = os.path.join(
    path_case_study_transcriptomics, 'MT_Depths')
path_case_study_proteomics_fasta = os.path.join(
    path_case_study_proteomics, 'set_of_78')
path_case_study_metabolomics_prec_1 = os.path.join(
    path_case_study_metabolomics, CALCULATED_DATA_SET_NAME)
path_case_study_metabolomics_prec_2 = os.path.join(
    path_case_study_metabolomics, CALCULATED_NOW_DATA_SET_NAME)
path_case_study_phy_che_prec_1 = os.path.join(
    path_case_study_phy_che, CALCULATED_DATA_SET_NAME)


def upload_multiple(key_suffix):

    available_data_set_types = {
        'Metagenomics': {
            'Raw FASTA files': 'FASTA',
            'BINS annotation files': 'BINS',
            'Depth-of-coverage': 'DEPTH'},
        'Metaproteomics': {
            'Raw FASTA files': 'FASTA'},
        'Metatranscriptomics': {
            'Depth-of-coverage': 'DEPTH'},
        'Metabolomics': {
            'Processed data set 1': 'CALC',
            'Processed data set 2': 'CALC'},
        'Physico-chemical': {
            'Processed data set 1': 'CALC'}
    }

    default_case_study_data_set_types = {
        'Metagenomics': 2,
        'Metaproteomics': 0,
        'Metatranscriptomics': 0,
        'Metabolomics': 1,
        'Physico-chemical': 0
    }

    selected_data_set_type = st.selectbox(
        'What kind of data set do you want to see?',
        options=list(available_data_set_types[key_suffix].keys()),
        index=default_case_study_data_set_types[key_suffix],
        key='Case_study_' + key_suffix)

    if key_suffix == 'Metagenomics':
        if selected_data_set_type == 'Raw FASTA files':
            return_path = path_case_study_genomics_fasta
        # elif selected_data_set_type == 'KEGG annotation files':
        #     return_path = path_case_study_genomics_kegg
        elif selected_data_set_type == 'Depth-of-coverage':
            return_path = path_case_study_genomics_depths
        else:
            return_path = path_case_study_genomics_bins

    elif key_suffix == 'Metaproteomics':
        return_path = path_case_study_proteomics_fasta

    elif key_suffix == 'Metatranscriptomics':
        return_path = path_case_study_transcriptomics_depths

    elif key_suffix == 'Metabolomics':
        if selected_data_set_type == 'Processed data set 1':
            return_path = path_case_study_metabolomics_prec_1
        elif selected_data_set_type == 'Processed data set 2':
            return_path = path_case_study_metabolomics_prec_2
        else:
            pass

    elif key_suffix == 'Physico-chemical':
        if selected_data_set_type == 'Processed data set 1':
            return_path = path_case_study_phy_che_prec_1
        else:
            pass

    else:
        pass

    return (return_path,
            available_data_set_types[key_suffix][selected_data_set_type])


def upload_intro(folder_path, key_suffix):
    st.header(key_suffix + ' data')
    st.markdown('')

    return_path = None
    return_path, data_set_type = upload_multiple(key_suffix)

    if return_path is None:
        st.warning('Upload your data set')

    # We return DataFrame if we work with tabular data format or precalculated
    # We return folder_path if we work with archived data
    # Data_set_type is always returned
    if data_set_type == 'CALC':
        return_path_or_df = common.get_cached_dataframe(return_path)
    else:
        return_path_or_df = return_path

    return return_path_or_df, data_set_type


def case_study_genomics():
    key_suffix = 'Metagenomics'
    cache_folder_path = path_case_study_genomics

    folder_path_or_df, data_set_type = upload_intro(
        cache_folder_path, key_suffix)
    key_suffix += '_CASE_STUDY'

    return common.work_with_zip(
        folder_path_or_df, data_set_type, cache_folder_path, key_suffix)


def case_study_proteomics():
    key_suffix = 'Metaproteomics'
    cache_folder_path = path_case_study_proteomics

    folder_path_or_df, data_set_type = upload_intro(
        cache_folder_path, key_suffix)
    key_suffix += '_CASE_STUDY'

    return common.work_with_zip(
        folder_path_or_df, data_set_type, cache_folder_path, key_suffix)


def case_study_metabolomics():
    key_suffix = 'Metabolomics'
    cache_folder_path = path_case_study_metabolomics

    folder_path_or_df, data_set_type = upload_intro(
        cache_folder_path, key_suffix)
    key_suffix += '_CASE_STUDY'

    return common.work_with_csv(
        folder_path_or_df, cache_folder_path, key_suffix)


def case_study_transcriptomics():
    key_suffix = 'Metatranscriptomics'
    cache_folder_path = path_case_study_transcriptomics

    folder_path_or_df, data_set_type = upload_intro(
        cache_folder_path, key_suffix)
    key_suffix += '_CASE_STUDY'

    return common.work_with_zip(
        folder_path_or_df, data_set_type, cache_folder_path, key_suffix)


def case_study_phy_che():
    key_suffix = 'Physico-chemical'
    cache_folder_path = path_case_study_phy_che

    folder_path_or_df, data_set_type = upload_intro(
        cache_folder_path, key_suffix)
    key_suffix += '_CASE_STUDY'

    return common.work_with_csv(
        folder_path_or_df, cache_folder_path, key_suffix)


def create_main_case_study():

    col_1, col_2 = st.columns([1, 2])

    col_1.info('''
            This data set comes from the following paper:

            **Herold, M., Martínez Arbas, S., Narayanasamy, S. et al.
            Integration of time-series meta-omics data reveals how microbial
            ecosystems respond to disturbance. Nat Commun 11, 5281(2020).
            https://doi.org/10.1038/s41467-020-19006-2**. Analyzed samples
            were collected from a biological wastewater treatment plant in
            Schifflange, Luxembourg (49.513414, 6.017925). A precise location
            is shown on the map located on the right.

            It contains **metagenomics**, **metabolomics**, **metaproteomics**,
            and **physico-chemical** data. The code used to parse the data can
            be found here: [GitLab]
            (https://git-r3lab.uni.lu/malte.herold/laots_niche_ecology_analysis)
            ''')

    col_2.map(pd.DataFrame({'lat': [49.513414], 'lon': [6.017925]}),
              zoom=8, use_container_width=True)

    case_study_omics_list = ['Metagenomics', 'Metabolomics', 'Metaproteomics',
                             'Metatranscriptomics', 'Physico-chemical']
    # We use this specific omics for case study, as default values for
    # multiselect widget
    default_case_study_omics = ['Physico-chemical', 'Metaproteomics',
                                'Metabolomics', 'Metagenomics']
    choose_omics = st.multiselect(
        'What kind of omic data do you want to explore?',
        case_study_omics_list, default=default_case_study_omics)

    num_of_columns = len(choose_omics)

    charts = []  # An empty list to hold all pairs (visualizations, key)

    if num_of_columns >= 2:
        column_list = st.columns(num_of_columns)
        curr_pos = 0

        for i in choose_omics:
            if i == 'Metagenomics':
                with column_list[curr_pos]:
                    curr_pos += 1
                    charts += case_study_genomics()
            elif i == 'Metabolomics':
                with column_list[curr_pos]:
                    curr_pos += 1
                    charts += case_study_metabolomics()
            elif i == 'Metaproteomics':
                with column_list[curr_pos]:
                    curr_pos += 1
                    charts += case_study_proteomics()
            elif i == 'Metatranscriptomics':
                with column_list[curr_pos]:
                    curr_pos += 1
                    charts += case_study_transcriptomics()
            elif i == 'Physico-chemical':
                with column_list[curr_pos]:
                    curr_pos += 1
                    charts += case_study_phy_che()
            else:
                pass

    else:
        for i in choose_omics:
            if i == 'Metagenomics':
                charts += case_study_genomics()
            elif i == 'Metabolomics':
                charts += case_study_metabolomics()
            elif i == 'Metaproteomics':
                charts += case_study_proteomics()
            elif i == 'Metatranscriptomics':
                charts += case_study_transcriptomics()
            elif i == 'Physico-chemical':
                charts += case_study_phy_che()
            else:
                pass

    st.markdown('---')

    for i in charts:
        type_of_chart = type(i[0])

        with st.spinner('Visualizing...'):
            if 'altair' in str(type_of_chart):
                st.altair_chart(i[0], use_container_width=True)
                common.save_chart(i[0], path_case_study_viz, i[1])
            else:
                pass

    return None
