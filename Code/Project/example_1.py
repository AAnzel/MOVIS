import os
import io
import random
import math
import gensim
import altair_saver
import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
import datetime as dt
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from gensim.models import Word2Vec
from sklearn.cluster import KMeans, OPTICS
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn import preprocessing, model_selection, metrics
from scipy.spatial.distance import jaccard, pdist, squareform

import omics_run

def create_main_shared():
    
    # TODO: Be careful when placing values in beta_columns
    # These are hardcoded values for column width
    tmp_col_1, tmp_col_2, tmp_col_3 = st.beta_columns([1.5, 2, 1])
    tmp_col_2.title('Tool title - Multi-omics time series')
    st.markdown(' ')
    st.markdown(' ')
    st.markdown('---')

    return None

def show_data_set(df):
    with st.beta_expander('Show the data set and related info', expanded = True):
        st.markdown('First 100 entries')
        st.dataframe(df.head(100))
        st.dataframe(df.describe())

    return None

def create_main_example_1_genomics():
    st.header ('Genomics')
    with st.beta_expander('Show folder structure', expanded = True):
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

    # I should put cluster charts here, however I have to run it first
    # because I have rendered images and not altair charts
    #st.altair_chart()
    
    return None


def create_main_example_1_metabolomics():
    st.header ('Metabolomics')

    # Here I show the head() of the data set and some summary() and info()
    df = omics_run.get_cached_dataframe(omics_run.EX_1, 'metabolomics')
    show_data_set(df)    

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return None

def create_main_example_1_proteomics():
    st.header ('Proteomics')

    # Here I show the head() of the data set and some summary() and info()

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return None

def create_main_example_1_phy_che():
    st.header ('Physico-chemical')

    # Here I show the head() of the data set and some summary() and info()
    df = omics_run.get_cached_dataframe(omics_run.EX_1, 'phy_che')
    show_data_set(df)
    
    with st.beta_expander('Show a correlation matrix', expanded = True):
        st.altair_chart(omics_run.phy_che.visualize_phy_che_heatmap(df),
                        use_container_width = True)



    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return None

def create_main_example_1():
    create_main_shared()
    
    st.info('''
            This data set comes from the following paper:

            **Herold, M., Martínez Arbas, S., Narayanasamy, S. et al. Integration\
            of time-series meta-omics data reveals how microbial ecosystems\
            respond to disturbance. Nat Commun 11, 5281 (2020).\
            https://doi.org/10.1038/s41467-020-19006-2**

            It contains **genomics**, **metabolomics**, **proteomics**, and\
            **physico-chemical** data. The code used to parse the data can be\
            found here: [GitHub](put_link)
            ''')
    
    example_1_omics_dict = {'Genomics':0, 'Metabolomics':0, 'Proteomics':0,
                            'Physico-chemical':0}
    choose_omics = st.multiselect('Which omic do you want to see:',
                                list(example_1_omics_dict.keys()))
    
    num_of_columns = 0
    for i in choose_omics:
        example_1_omics_dict[i] = 1
        num_of_columns += 1
    
    if num_of_columns >= 2:
        column_list = st.beta_columns(num_of_columns)
        curr_pos = 0

        for i in list(example_1_omics_dict.keys()):
            if example_1_omics_dict[i] != 0:
                if i == 'Genomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        create_main_example_1_genomics()
                elif i == 'Metabolomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        create_main_example_1_metabolomics()
                elif i == 'Proteomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        create_main_example_1_proteomics()
                else:
                    with column_list[curr_pos]:
                        curr_pos += 1
                        create_main_example_1_phy_che()
    
    else:
        for i in list(example_1_omics_dict.keys()):
            if example_1_omics_dict[i] != 0:
                if i == 'Genomics':
                    create_main_example_1_genomics()
                elif i == 'Metabolomics':
                    create_main_example_1_metabolomics()
                elif i == 'Proteomics':
                    create_main_example_1_proteomics()
                else:
                    create_main_example_1_phy_che()
    return None
