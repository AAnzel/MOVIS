import os
import streamlit as st
import pandas as pd
import common


path_example_2_root_data = os.path.join('..', 'Data', 'cached', 'example_2')
path_example_2_transcriptomics = os.path.join(path_example_2_root_data,
                                              'transcriptomics')


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


def example_2_transcriptomics():
    key_suffix = 'Transcriptomics'
    cache_folder_path = path_example_2_transcriptomics

    df = upload_intro(cache_folder_path, key_suffix)

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
