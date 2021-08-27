import os
import streamlit as st
import pandas as pd
import common


path_example_2_root_data = os.path.join('..', 'Data', 'cached', 'example_2')
path_example_2_transcriptomics = os.path.join(path_example_2_root_data,
                                              'transcriptomics')
path_example_2_transcriptomics_control_1 = os.path.join(
    path_example_2_transcriptomics, 'CONTROL_1.csv')
path_example_2_transcriptomics_control_2 = os.path.join(
    path_example_2_transcriptomics, 'CONTROL_2.csv')
path_example_2_transcriptomics_control_3 = os.path.join(
    path_example_2_transcriptomics, 'CONTROL_3.csv')

path_example_2_transcriptomics_benz_1 = os.path.join(
    path_example_2_transcriptomics, 'BENZ_1.csv')
path_example_2_transcriptomics_benz_2 = os.path.join(
    path_example_2_transcriptomics, 'BENZ_2.csv')
path_example_2_transcriptomics_benz_3 = os.path.join(
    path_example_2_transcriptomics, 'BENZ_3.csv')

path_example_2_transcriptomics_phe_1 = os.path.join(
    path_example_2_transcriptomics, 'PHE_1.csv')
path_example_2_transcriptomics_phe_2 = os.path.join(
    path_example_2_transcriptomics, 'PHE_1.csv')
path_example_2_transcriptomics_phe_3 = os.path.join(
    path_example_2_transcriptomics, 'PHE_3.csv')

path_example_2_transcriptomics_pov_1 = os.path.join(
    path_example_2_transcriptomics, 'POV_1.csv')
path_example_2_transcriptomics_pov_2 = os.path.join(
    path_example_2_transcriptomics, 'POV_2.csv')
path_example_2_transcriptomics_pov_3 = os.path.join(
    path_example_2_transcriptomics, 'POV_3.csv')


def upload_multiple(key_suffix):

    available_data_sets = {
        'Control replicate 1': path_example_2_transcriptomics_control_1,
        'Control replicate 2': path_example_2_transcriptomics_control_2,
        'Control replicate 3': path_example_2_transcriptomics_control_3,
        'Benz. chl. replicate 1': path_example_2_transcriptomics_benz_1,
        'Benz. chl. replicate 2': path_example_2_transcriptomics_benz_2,
        'Benz. chl. replicate 3': path_example_2_transcriptomics_benz_3,
        'Pov. iod. replicate 1': path_example_2_transcriptomics_pov_1,
        'Pov. iod. replicate 2': path_example_2_transcriptomics_pov_2,
        'Pov. iod. replicate 3': path_example_2_transcriptomics_pov_3,
        'Chloroph. replicate 1': path_example_2_transcriptomics_phe_1,
        'Chloroph. replicate 2': path_example_2_transcriptomics_phe_2,
        'Chloroph. replicate 3': path_example_2_transcriptomics_phe_3
    }

    selected_data_sets = st.multiselect(
        'What data sets do you want to see?', list(available_data_sets.keys()),
        default=['Control replicate 1', 'Control replicate 2',
                 'Control replicate 3'])

    df_list = []
    selected_df_names = []

    for data_set in selected_data_sets:
        tmp_df = pd.read_csv(available_data_sets[data_set])
        tmp_df = common.fix_dataframe_columns(tmp_df)
        tmp_df = tmp_df.convert_dtypes()

        df_list.append(tmp_df)
        selected_df_names.append(data_set)

    return df_list, selected_df_names


def upload_intro(folder_path, key_suffix):
    st.header(key_suffix + ' data')
    st.markdown('')
    st.info('This data set has many features and takes time to load')
    st.markdown('')

    df_list = None
    df_list, selected_df_names = upload_multiple(key_suffix)

    if df_list is None:
        st.warning('Choose your data set(s)')
        st.stop()

    return df_list, selected_df_names


def example_2_transcriptomics():
    key_suffix = 'Transcriptomics'
    cache_folder_path = path_example_2_transcriptomics

    df_list, selected_df_names = upload_intro(cache_folder_path, key_suffix)

    if len(df_list) > 1:
        df = common.work_with_multi_transcriptomics(
            df_list, selected_df_names)

    return common.work_with_csv(df, cache_folder_path, key_suffix)


def create_main_example_2():

    col_1, col_2 = st.beta_columns([1, 2])

    col_1.info('''
            This data set comes from the following paper:
            **Merchel Piovesan Pereira, B., Wang, X., & Tagkopoulos, I. (2020).
            Short- and Long-Term Transcriptomic Responses of Escherichia coli
            to Biocides: a Systems Analysis. Applied and environmental
            microbiology, 86(14), e00708-20.
            https://doi.org/10.1128/AEM.00708-20**. Analyzed samples
            were processed in DNA Technologies Core, Genome and Biomedical
            Sciences Facility (GBSF) (38.535244, -121.764920).
            A precise location is shown on the map located on the right.

            It contains **transcriptomics** data.
            ''')

    col_2.map(pd.DataFrame({'lat': [38.53524478359195],
                            'lon': [-121.76492092285179]}),
              zoom=8, use_container_width=True)

    example_2_omics_list = ['Transcriptomics']
    choose_omics = st.multiselect(
        'What kind of data set do you want to see?', example_2_omics_list)

    charts = []  # An empty list to hold all pairs (visualizations, key)

    with st.beta_expander('Show/hide data sets and related info',
                          expanded=True):
        for i in choose_omics:
            if i == 'Transcriptomics':
                charts += example_2_transcriptomics()
            else:
                pass

    with st.beta_expander('Show/hide visualizations', expanded=True):
        for i in charts:
            type_of_chart = type(i[0])

            with st.spinner('Visualizing...'):
                if 'altair' in str(type_of_chart):
                    st.altair_chart(i[0], use_container_width=True)
                else:
                    pass

    return None
