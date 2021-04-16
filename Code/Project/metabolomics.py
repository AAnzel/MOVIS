import os
import random
import math
import streamlit
import gensim
import altair_saver
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


# Defining paths for each and every omic

path_root_data = os.path.join ('..', '..', 'Data', 'Extracted', 'First source', 'Databases')

path_all_fasta = os.path.join (path_root_data, 'fasta_files', 'AllBins')
path_genomics_78 = os.path.join (path_root_data, 'fasta_files', 'rmags_filtered')
path_genomics_kegg = os.path.join (path_root_data, 'Annotations', 'KEGG')
path_normalised_metabolomics = os.path.join (path_root_data, 'Metabolomics', 'Normalised_Tables')
path_proteomics_78 = os.path.join (path_root_data, 'Proteomics', 'set_of_78')
path_physico_chemical = os.path.join (path_root_data, 'PhysicoChemical')
path_second_source = os.path.join ('..', '..', 'Data', 'Extracted', 'Second source')

path_model_save_root = os.path.join ('..', 'Saved_models')
path_figures_save_root = os.path.join ('..', 'Output_figures')

num_of_mags = len([i for i in os.listdir(path_genomics_78) if i.endswith('fa')])
num_of_proteomics = len([i for i in os.listdir(path_proteomics_78) if i.endswith('faa')])
SEED = 42
END = num_of_mags
ALL_DAYS = 51
MAX_ROWS = 15000
EPOCHS = 10
NUM_OF_WORKERS = 8
START_DATE = dt.datetime.strptime ('2011-03-21', '%Y-%m-%d')
random.seed(SEED)
np.random.seed(SEED)
alt.data_transformers.enable('default', max_rows = MAX_ROWS) # Important if you want to visualize datasets with >5000 samples


# Functions below are shared among different omics
# Function that saves charts from list_of_charts with names from list_of_names
def save_charts (list_of_chart, list_of_names):
    
    for chart, name in zip(list_of_chart, list_of_names):
        print (chart, name)
        #altair_saver.save(chart, os.path.join (path_figures_save_root, name))
        chart.save(os.path.join (path_figures_save_root, name))

# This function creates new dataframe with column that represent season according to date
# It also concatenates important types with metabolite names
def season_data (data, temporal_column):
    new_df = data
    new_df['season'] = new_df[temporal_column].dt.month%12 // 3 + 1
    
    #important_types = [metabolite_column] + important_types
    #new_df['new_name'] = df[important_types].agg('\n'.join, axis=1)
    
    return new_df

def create_temporal_column (list_of_days, start_date, end):
    
    list_of_dates = []
    
    # This is specific to the metaomics data set I am using
    # Creating list of dates for every rMAG
    for i in list_of_days[:end]:
        
        tmp_datetime = start_date + dt.timedelta (weeks = int(i[1:3]))
        
        if tmp_datetime not in list_of_dates:
            list_of_dates.append (tmp_datetime)
        
        else:
            tmp_datetime = tmp_datetime.replace (day = tmp_datetime.day + 1)
            list_of_dates.append (tmp_datetime)
    
    return list_of_dates

# This functions are used for METABOLOMICS

def visualize_metabolites (data, temporal_column, metabolite_column, type_columns):
    
    data_seasoned = season_data (data, temporal_column)
    
    # Extract columns with float values
    float_columns = []
    
    for i in data_seasoned.columns:
        if data_seasoned[i].dtypes == 'float64' or data_seasoned[i].dtypes == 'float32':
            float_columns.append(i)
    
    # Create repeated chart with varying size encodings
    chart = alt.Chart(data_seasoned).mark_point(opacity = 1).encode(
        alt.X (temporal_column, type = 'temporal', scale = alt.Scale (nice = True)),
        alt.Y (metabolite_column, type = 'nominal'),
        alt.Size (alt.repeat("row"), type = 'quantitative'),
        alt.Color ('season:N', scale = alt.Scale (range = ['blue', 'green', 'orange', 'brown'])),
        alt.Tooltip (type_columns, type = 'nominal')
    ).properties(
        width = 1200
    ).repeat(
        row = float_columns
    ).resolve_scale(color = 'independent', size = 'independent')#.interactive()
    
    return chart
