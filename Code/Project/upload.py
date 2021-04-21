import os
import streamlit as st
import pandas as pd
import altair as alt

URANDOM_LENGTH = 5

def show_data_set(df):
    with st.beta_expander('Show the data set and related info', expanded = True):
        st.markdown('First 100 entries')
        st.dataframe(df.head(100))
        st.dataframe(df.describe())

    return None

def upload_data_set():

    # I should put some random number as a key, because I get errors
    imported_file = st.file_uploader('Upload your data set here. Maximum size\
                                     is 200MB.', type = ['csv', 'xlsx', 'xls'],
                                     key = os.urandom(URANDOM_LENGTH))
    delimiter = st.selectbox('Select the delimiter in your data set',
                            ['Comma (,)', 'Semicolon (;)', 'Tab (\\t)', 'Excel File'],
                            key = os.urandom(URANDOM_LENGTH))
    
    if imported_file is not None:
        return pd.read_csv(imported_file, delimiter = delimiter)
    else:
        return None

def upload_genomics():
    st.header ('Genomics')
    st.markdown ('')

    df = upload_data_set()
    
    if df is None:
        return None

    show_data_set(df)

    return None

def upload_proteomics():
    st.header ('Proteomics')
    st.markdown ('')

    df = upload_data_set()

    if df is None:
        return None

    show_data_set(df)

    return None

def upload_metabolomics():
    st.header ('Metabolomics')
    st.markdown ('')

    df = upload_data_set()
    
    if df is None:
        return None
    
    show_data_set(df)

    return None

def upload_transcriptomics():
    st.header ('Transcriptomics')
    st.markdown ('')

    df = upload_data_set()
    
    if df is None:
        return None
    
    show_data_set(df)

    return None

def upload_phy_che():
    st.header ('Physico-chemical')
    st.markdown ('')

    df = upload_data_set()
    
    if df is None:
        return None
    
    show_data_set(df)

    return None

def create_main_upload():

    st.header('Dataset')

    upload_omics_dict = {'Genomics':0, 'Metabolomics':0, 'Proteomics':0,
                            'Physico-chemical':0, 'Transcriptomics':0}
    choose_omics = st.multiselect('Which omic data do you want to upload:',
                                list(upload_omics_dict.keys()))
    
    num_of_columns = 0
    for i in choose_omics:
        upload_omics_dict[i] = 1
        num_of_columns += 1
    
    if num_of_columns >= 2:
        column_list = st.beta_columns(num_of_columns)
        curr_pos = 0

        for i in list(upload_omics_dict.keys()):
            if upload_omics_dict[i] != 0:
                if i == 'Genomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        upload_genomics()
                elif i == 'Metabolomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        upload_metabolomics()
                elif i == 'Proteomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        upload_proteomics()
                elif i == 'Transcriptomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        upload_transcriptomics()
                else:
                    with column_list[curr_pos]:
                        curr_pos += 1
                        upload_phy_che()
    
    else:
        for i in list(upload_omics_dict.keys()):
            if upload_omics_dict[i] != 0:
                if i == 'Genomics':
                    upload_genomics()
                elif i == 'Metabolomics':
                    upload_metabolomics()
                elif i == 'Proteomics':
                    upload_proteomics()
                elif i == 'Transcriptomics':
                    upload_transcriptomics()
                else:
                    upload_phy_che()
    
    return None
