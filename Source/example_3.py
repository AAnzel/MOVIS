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


def create_main_example_3():

    col_1, col_2 = st.beta_columns([1, 2])

    col_1.info('''
            This data set comes from the following paper:
            **Zak Costello and Hector Garcia Martin
            A machine learning approach to predict metabolic pathway dynamics
            from time-series multiomics data. npj Syst Biol Appl 4, 19 (2018).
            https://doi.org/10.1038/s41540-018-0054-3**. Analyzed samples
            were processed in Powell-Focht Bioengineering Hall
            (32.881803, -117.233881).
            A precise location is shown on the map located on the right.

            It contains **metabolomics** and **proteomics** data. The code used
            to parse the data can be found here: [GitHub]
            (https://github.com/JBEI/KineticLearning)
            ''')

    col_2.map(pd.DataFrame({'lat': [32.88180365509607],
                            'lon': [-117.23388140034122]}),
              zoom=8, use_container_width=True)

    example_3_omics_list = ['Metabolomics', 'Proteomics']
    choose_omics = st.multiselect(
        'What kind of data set do you want to see?', example_3_omics_list)

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
    #                     charts += example_3_proteomics()
    #             else:
    #                 with column_list[curr_pos]:
    #                     curr_pos += 1
    #                     charts += example_3_metabolomics()

    #     else:
    #         for i in choose_omics:
    #             if i == 'Proteomics':
    #                 charts += example_3_proteomics()
    #             else:
    #                 charts += example_3_metabolomics()

    with st.beta_expander('Show/hide visualizations', expanded=True):
        for i in charts:
            type_of_chart = type(i[0])

            with st.spinner('Visualizing...'):
                if 'altair' in str(type_of_chart):
                    st.altair_chart(i[0], use_container_width=True)
                else:
                    pass

    return None
