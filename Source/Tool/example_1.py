import streamlit as st
import pandas as pd
import common


def get_data_set(omic_name):

    spinner_text_map =\
        {'genomics_mags': 'Embedding MAGs into vectors...',
         'genomics_mags_temporal': 'Embedding MAGs into vectors...',
         'genomics_mags_temporal_PCA': 'Embedding MAGs into vectors...',
         'genomics_kegg_temporal': 'Creating KO dataset...',
         'genomics_mags_annotated_temporal': 'Creating annotation dataset...',
         'genomics_mags_top_10_annotated_temporal': 'Creating annotation\
                                                     dataset...',
         'proteomics': 'Creating additional data set...',
         'metabolomics': 'Creating additional data set...',
         'phy_che': 'Creating additional data set...'}

    with st.spinner(spinner_text_map[omic_name]):
        df = common.get_cached_dataframe(common.EX_1, omic_name)

    # Fixing dataframe columns
    df = common.fix_dataframe_columns(df)

    return df


def example_1_genomics():
    st.header('Genomics')
    with st.spinner('Showing folder structure'):
        st.code('''rmags_filtered/
            ├── D03_O1.31.2.fa
            ├── D04_G2.13.fa
            ├── D04_G2.5.fa
            ├── D04_L6.fa
            ├── D04_O2.19.fa
            ├── D05_G3.14.1.fa
            ├── D05_G3.4.fa
            ├── D05_L3.17.fa
            ├── D08_G1.16.fa
            ├── D08_O6.fa
            ...
                ''')

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants
    data_set_list = ['W2V embedded MAGs', 'KEGG matrix',
                     'Product-annotated MAGs']
    choose_data_set = st.multiselect('Which data set do you want to see:',
                                     data_set_list)

    chosen_charts = []
    for i in choose_data_set:
        if i == 'W2V embedded MAGs':
            df_1 = get_data_set('genomics_mags_temporal')
            common.show_calculated_data_set(df_1, 'Embedded MAGs')
            labels_list = common.show_clustering_info(df_1, 'Genomics_1')

            # Traversing pairs in list
            for i in labels_list:
                temporal_feature, feature_list = common.find_temporal_feature(
                    df_1)
                feature_list = i[0]
                df_1[i[0]] = i[1]
                chosen_charts += common.visualize_data_set(
                        df_1, temporal_feature, feature_list,
                        'Cluster_Genomics_' + i[0])

        elif i == 'KEGG matrix':
            df_2 = get_data_set('genomics_kegg_temporal')
            common.show_calculated_data_set(df_2, 'KO matrix')
            labels_list = common.show_clustering_info(df_2, 'Genomics_2')

            # Traversing pairs in list
            for i in labels_list:
                temporal_feature, feature_list = common.find_temporal_feature(
                    df_2)
                feature_list = i[0]
                df_2[i[0]] = i[1]
                chosen_charts += common.visualize_data_set(
                        df_2, temporal_feature, feature_list,
                        'Cluster_Genomics_' + i[0])

        else:
            tmp_df_3 = get_data_set('genomics_mags_annotated_temporal')
            df_3 = get_data_set('genomics_mags_top_10_annotated_temporal')
            common.show_calculated_data_set(tmp_df_3, 'Product annotations')

            temporal_feature, feature_list = common.find_temporal_feature(df_3)
            chosen_charts += common.visualize_data_set(
                df_3, temporal_feature, feature_list, 'Genomics_3')

    # I should put cluster charts here, however I have to run it first
    # because I have rendered images and not altair charts
    # st.altair_chart()

    return chosen_charts


def example_1_metabolomics():
    st.header('Metabolomics')

    # Here I show the head() of the data set and some summary() and info()
    df = get_data_set('metabolomics')
    common.show_data_set(df)

    temporal_feature, feature_list = common.find_temporal_feature(df)

    chosen_charts = common.visualize_data_set(
        df, temporal_feature, feature_list, 'Metabolomics')

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return chosen_charts


def example_1_proteomics():
    st.header('Proteomics')

    # TODO: Implement the same data set creation as with the genomics data
    # Create a data set with w2v and then implement all data set visulizations
    # as with the genomics data
    with st.spinner('Showing folder structure'):
        st.code('''set_of_78/
            ├── D03_O1.31.2.faa
            ├── D04_G2.13.faa
            ├── D04_G2.5.faa
            ├── D04_L6.faa
            ├── D04_O2.19.faa
            ├── D05_G3.14.1.faa
            ├── D05_G3.4.faa
            ├── D05_L3.17.faa
            ├── D08_G1.16.faa
            ├── D08_O6.faa
            ...
                ''')

    # Here I show the head() of the data set and some summary() and info()
    df = get_data_set('proteomics')
    common.show_calculated_data_set(df, 'Protein properties')

    temporal_feature, feature_list = common.find_temporal_feature(df)

    chosen_charts = common.visualize_data_set(
        df, temporal_feature, feature_list, 'Proteomics')

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return chosen_charts


def example_1_phy_che():
    st.header('Physico-chemical')

    # Here I show the head() of the data set and some summary() and info()
    df = get_data_set('phy_che')
    common.show_data_set(df)

    temporal_feature, feature_list = common.find_temporal_feature(df)

    chosen_charts = common.visualize_data_set(
        df, temporal_feature, feature_list, 'Physico-Chemical')

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return chosen_charts


def example_1():

    col_1, col_2 = st.beta_columns(2)

    col_1.info('''
            This data set comes from the following paper:

            **Herold, M., Martínez Arbas, S., Narayanasamy, S. et al.\
            Integration of time-series meta-omics data reveals how microbial\
            ecosystems respond to disturbance. Nat Commun 11, 5281(2020).\
            https://doi.org/10.1038/s41467-020-19006-2**. Analyzed samples\
            were collected from a biological wastewater treatment plant in\
            Schifflange, Luxembourg. A precise location is shown on the map\
            located on the right.

            It contains **genomics**, **metabolomics**, **proteomics**, and\
            **physico-chemical** data. The code used to parse the data can be\
            found here: [GitLab]
            (https://git-r3lab.uni.lu/malte.herold/laots_niche_ecology_analysis)
            ''')

    col_2.map(pd.DataFrame({'lat': [49.513414], 'lon': [6.017925]}),
              zoom=13, use_container_width=True)

    example_1_omics_list = ['Genomics', 'Metabolomics', 'Proteomics',
                            'Physico-chemical']
    choose_omics = st.multiselect('Which omic do you want to see:',
                                  example_1_omics_list)

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
                        charts += example_1_genomics()
                elif i == 'Metabolomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        charts += example_1_metabolomics()
                elif i == 'Proteomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        charts += example_1_proteomics()
                else:
                    with column_list[curr_pos]:
                        curr_pos += 1
                        charts += example_1_phy_che()

        else:
            for i in choose_omics:
                if i == 'Genomics':
                    charts += example_1_genomics()
                elif i == 'Metabolomics':
                    charts += example_1_metabolomics()
                elif i == 'Proteomics':
                    charts += example_1_proteomics()
                else:
                    charts += example_1_phy_che()

    with st.beta_expander('Show/hide visualizations', expanded=True):
        for i in charts:
            type_of_chart = type(i[0])

            with st.spinner('Visualizing...'):
                if 'altair' in str(type_of_chart):
                    st.altair_chart(i[0], use_container_width=True)
                else:
                    pass

    return None
