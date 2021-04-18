import os
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


def create_main_shared():
    
    # TODO: Be careful when placing values in beta_columns
    # These are hardcoded values for column width
    tmp_col_1, tmp_col_2, tmp_col_3 = st.beta_columns([1.8, 2, 1])
    tmp_col_2.title('Tool title - ML cruncher')
    st.markdown('---')

    return None
    
def create_main_example_1_genomics():
    st.subheader ('Genomics')
    with st.beta_expander('See folder structure', expanded = True):
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
    st.subheader ('Metabolomics')

    # Here I show the head() of the data set and some summary() and info()

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return None

def create_main_example_1_proteomics():
    st.subheader ('Proteomics')

    # Here I show the head() of the data set and some summary() and info()

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return None

def create_main_example_1_phy_che():
    st.subheader ('Physico-chemical')

    # Here I show the head() of the data set and some summary() and info()

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
    choose_omics = st.multiselect('Which omics do you want to see:',
                                list(example_1_omics_dict.keys()))
    
    num_of_columns = 0
    for i in choose_omics:
        print (i)
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
    

def create_main_example_2():
    create_main_shared()

    return None

def create_main_upload():
    create_main_shared()

    st.header('Dataset')

    imported_file = st.file_uploader('Upload your data set here. Maximum size\
                                     is 200MB.', type = ['csv', 'xlsx', 'xls'])
    delimiter = st.selectbox('Select the delimiter in your data set',
                            ['Comma (,)', 'Semicolon (;)', 'Tab (\\t)', 'Excel File'])
    
    return None


def create_sidebar():
    st.sidebar.markdown('Above should be a logo. Below should be a github repo\
                        logo that is clickable, and on the right should be\
                        the link to the paper.')
    
    # We need 5 columns so that we have nicely aligned images
    col_1, col_2, col_3, col_4, col_5 = st.sidebar.beta_columns([1, 2, 1, 2, 1])
    col_2.image(os.path.join('images', 'GitHub-Mark-120px-plus.png'), width = 52)
    col_2.markdown('[GitHub](https://github.com/AAnzel/Multi-omics_platform)')
    col_4.markdown('Paper doi with journal logo')
    
    st.sidebar.markdown('---')
    st.sidebar.markdown('Here put some info about the app.')
    st.sidebar.markdown('---')

    st.sidebar.markdown('**Choose your data set:**')
    choice_data_set = st.sidebar.radio('', ('Example 1', 'Example 2', 'Upload'),
                                       index = 0)

    if choice_data_set == 'Example 1':
        create_main_example_1()
    elif choice_data_set == 'Example 2':
        create_main_example_2()
    else:
        create_main_upload()
    
    return None



def create_main():
    
    create_sidebar()

    return None




def main():
    
    st.set_page_config(layout = 'wide')
    create_main()

    return None




if __name__ == '__main__':
    main()