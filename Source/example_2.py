import os
import streamlit as st
import pandas as pd
import common


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

    num_of_columns = len(choose_omics)

    charts = []  # An empty list to hold all pairs (visualizations, key)

    # with st.beta_expander('Show/hide data sets and related info',
    #                       expanded=True):
    #     if num_of_columns >= 2:
    #         column_list = st.beta_columns(num_of_columns)
    #         curr_pos = 0

    #         for i in choose_omics:
    #             if i == 'Proteomics':
    #                 with column_list[curr_pos]:
    #                     curr_pos += 1
    #                     charts += example_2_proteomics()
    #             else:
    #                 with column_list[curr_pos]:
    #                     curr_pos += 1
    #                     charts += example_2_metabolomics()

    #     else:
    #         for i in choose_omics:
    #             if i == 'Proteomics':
    #                 charts += example_2_proteomics()
    #             else:
    #                 charts += example_2_metabolomics()

    with st.beta_expander('Show/hide visualizations', expanded=True):
        for i in charts:
            type_of_chart = type(i[0])

            with st.spinner('Visualizing...'):
                if 'altair' in str(type_of_chart):
                    st.altair_chart(i[0], use_container_width=True)
                else:
                    pass

    return None
