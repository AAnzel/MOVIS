#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import altair as alt
import os, random, math
import gensim
import datetime as dt
import altair_saver
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from gensim.models import Word2Vec
from sklearn.cluster import KMeans, OPTICS
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn import preprocessing, model_selection, metrics
from scipy.spatial.distance import jaccard, pdist, squareform


# ## Defining paths for each and every omic

# In[2]:


path_root_data = os.path.join ('..', 'Data', 'Extracted', 'First source', 'Databases')

path_all_fasta = os.path.join (path_root_data, 'fasta_files', 'AllBins')
path_genomics_78 = os.path.join (path_root_data, 'fasta_files', 'rmags_filtered')
path_genomics_kegg = os.path.join (path_root_data, 'Annotations', 'KEGG')
path_normalised_metabolomics = os.path.join (path_root_data, 'Metabolomics', 'Normalised_Tables')
path_proteomics_78 = os.path.join (path_root_data, 'Proteomics', 'set_of_78')
path_physico_chemical = os.path.join (path_root_data, 'PhysicoChemical')
path_second_source = os.path.join ('..', 'Data', 'Extracted', 'Second source')

path_model_save_root = 'Saved_models'
path_figures_save_root = 'Output_figures'


# In[3]:


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


# ---
# # GENOMIC ANALYSIS
# ---
# ## MAG-related functions

# **Important**: I should review the way I look at MAGs. The names of all fasta files beggining with 'D_##' represent the days those MAGs were obtained. Therefore, I should look at this also as timeseries data. Also, maybe I should only consider 78 MAGs, and not all ~1300.
# After some consideration, I conclude that I should definetly use only 78 MAGs, because that way I wouldn't be tied to meta-omics data only. I also thinked about what should I visualize in that case. One idea is that I should also encode those MAGs with word2wec, and then make a 3D chart where one dimension is time, and other two dimensions would be PCA dimensions of those MAGs. I could also use this function to visualize proteomics data if I want.
# 
# Another important thing is that I should actually include FASTA headers and possibly use only them. That way, I could make figures like in a relevant paper where MAGs are groupped according to their taxonomy etc. I should look more into this.

# In[ ]:





# In[ ]:





# In[16]:


# Function that saves charts from list_of_charts with names from list_of_names
def save_charts (list_of_chart, list_of_names):
    
    for chart, name in zip(list_of_chart, list_of_names):
        print (chart, name)
        #altair_saver.save(chart, os.path.join (path_figures_save_root, name))
        chart.save(os.path.join (path_figures_save_root, name))
    


# In[4]:


# Function that splits each genome into k-mers thus creating even longer sentence (MAG)
# It returns tokenized genome i.e. [kmer, kmer,...]
def split_genome (genome, k = 5):
    new_genome = []
    n = len(genome)
    
    if n-k <=0:
        return genome
    else:
        for i in range(n-k):
            new_genome.append(genome[i:i+k])
        
        return new_genome


def vectorize_one_mag (one_mag, w2v_model):
    
    # We have to generate vectors for each word in one MAG and then create vector representation of that MAG
    # by averaging vectors of its words
    zero_vector = np.zeros(w2v_model.vector_size)
    word_vectors = []
    one_mag_vector = []
    
    for sentence in one_mag:
        for word in sentence:
            if word in w2v_model.wv:
                try:
                    word_vectors.append(w2v_model.wv[word])
                except KeyError:
                    print ('Key Error')
                    continue
    
    if word_vectors:
        word_vectors = np.asarray(word_vectors)
        one_mag_vector = word_vectors.mean (axis=0)
    
    else:
        one_mag_vector = zero_vector
    
    return one_mag_vector


# Function that vectorizes a MAG (document) with a pretrained word2vec model. It returns vector representation of a given MAG
# Vectorization is done by averaging word (k-mer) vectors for the whole document (MAG)
def vectorize_mags (w2v_model, path_fasta = path_genomics_78, end = 25):
    
    print ('Vectorizing MAGs')
    
    fasta_files = [i for i in os.listdir(path_fasta) if (i.endswith('fa') and i.startswith('D'))]
    list_of_mag_vectors = []
    
    # This was done so that I could work with first 'end' FASTA files only. Otherwise, I should just remove: i, and enumerate
    for i, fasta_file_name in enumerate(fasta_files):
        
        if i == end:
            break
        
        else:
            with open(os.path.join(path_fasta, fasta_file_name), 'r') as input_file:
                
                one_mag = []
                for fasta_string in SeqIO.parse(input_file, "fasta"):
                    
                    # Get kmers of a genome and create a sentence (list of words)
                    temp_kmers = split_genome (str(fasta_string.seq))
                    
                    # Create a document (list of sentences)
                    one_mag.append(temp_kmers)
                    
                # Vectorizing MAGs one by one
                list_of_mag_vectors.append (vectorize_one_mag (one_mag, w2v_model))
    
    print ('Finished vectorizing')
    
    return list_of_mag_vectors
    
    

# If one wants to import MAGs in order to vectorize them, one should use start argument in order to skip first 'start' MAGs
# If one wants to import MAGs to train word2vec model, one should use only end argument, so that first 'end' MAGs are used for training
# Todo: Implement randomisation in picking MAGs for training, and don't use first 'start' MAGs for training
# (create list of indexes from 0 to 1364, use sklearn split train test, then traverse directory and use only MAGs with indexes in train for training w2v)
def import_mags_and_build_model (end = 25, path_fasta = path_genomics_78):
    
    print ('Importing MAGs and building model')
    
    # There are 1364 MAGs enclosed in FASTA files
    # I have to traverse every FASTA file, and in each file every sequence
    
    fasta_files = [i for i in os.listdir(path_fasta) if (i.endswith('fa') and i.startswith('D'))]
    fasta_ids = []
    
    # This was done so that I could work with first 100 FASTA files only. Otherwise, I should just remove: i, and enumerate
    for i, fasta_file_name in enumerate(fasta_files):
        
        if i == end:
            break
        
        else:
            with open(os.path.join(path_fasta, fasta_file_name), 'r') as input_file:
                
                one_mag = []
                one_mag_ids = []
                for fasta_string in SeqIO.parse(input_file, "fasta"):
                    
                    # Get kmers of a genome and create a sentence (list of words)
                    temp_kmers = split_genome (str(fasta_string.seq))
                    
                    # Create a document (list of sentences)
                    one_mag.append(temp_kmers)
                    # Save FASTA ids for every MAG
                    one_mag_ids.append(str(fasta_string.id))
                    
                # Save list of ids for one MAG in global list
                fasta_ids.append(one_mag_ids)
                
                # If we do not have a model, we build one
                if i == 0:
                    print ('Building w2v model')
                    # We build our model on the first MAG
                    w2v_model = Word2Vec (sentences = one_mag, size = 100, workers = NUM_OF_WORKERS, seed=SEED)
                
                # Else we just expand its vocabulary
                else:
                    # Now we expand our vocabulary
                    w2v_model.build_vocab (one_mag, update = True)
                    
    print ('Finished building')
    
    return w2v_model, fasta_files, fasta_ids


def train_model (w2v_model, epochs, path_fasta = path_genomics_78, end = 25):
    
    print ('Starting model training')
    
    # There are 1364 MAGs enclosed in FASTA files
    # I have to traverse every FASTA file, and in each file every sequence
    
    fasta_files = [i for i in os.listdir(path_fasta) if (i.endswith('fa') and i.startswith('D'))]
    
    # This was done so that I could work with first 100 FASTA files only. Otherwise, I should just remove: i, and enumerate
    for i, fasta_file_name in enumerate(fasta_files):
        
        if i == end:
            break
        
        else:
            with open(os.path.join(path_fasta, fasta_file_name), 'r') as input_file:
                
                one_mag = []
                for fasta_string in SeqIO.parse(input_file, "fasta"):
                    
                    # Get kmers of a genome and create a sentence (list of words)
                    temp_kmers = split_genome (str(fasta_string.seq))
                    
                    # Create a document (list of sentences)
                    one_mag.append(temp_kmers)
                    
                
                w2v_model.train (one_mag, total_examples = w2v_model.corpus_count, epochs = epochs)
                    
    print ('Model training finished')
    
    return w2v_model


def visualize_with_pca (data, labels, centers):
    
    pca_model = PCA (n_components = 2, random_state = SEED)
    data_transformed = pca_model.fit_transform (data)
    
    data_transformed = pd.DataFrame (data_transformed)
    data_transformed.columns = ['PC_1', 'PC_2']
    data_transformed['Labels'] = labels
    
    chart_data = alt.Chart(data_transformed).mark_circle(opacity = 1).encode(
        alt.X ('PC_1:Q'),
        alt.Y ('PC_2:Q'),
        alt.Color ('Labels:N', legend = alt.Legend())
    )
    
    # This means we are visualising centroids from k_means (there are less centroids that data points)
    if labels.shape[0] != centers.shape[0]:
        
        centers_transformed = pca_model.fit_transform (centers)
        centers_transformed = pd.DataFrame (centers_transformed)
        centers_transformed.columns = ['PC_1', 'PC_2']
        
        chart_centers = alt.Chart(centers_transformed).mark_point(shape = 'diamond', color = 'black', size = 50, opacity = 0.7).encode(
            alt.X ('PC_1:Q'),
            alt.Y ('PC_2:Q'),
        )
        
        return chart_data + chart_centers
    
    # For DBSCAN there are no centroids
    else:
        return chart_data


# This function creates new dataframe with column that represent season according to date
# It also concatenate important types with metabolite names
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


def visualize_temporal_mags (data, list_of_days, start_date, end):
    
    list_of_dates = create_temporal_column (list_of_days, start_date, end)
    
    pca_model = PCA (n_components = 2, random_state = SEED)
    data_transformed = pca_model.fit_transform (data)
    
    data_transformed = np.hstack(((np.asarray(list_of_dates))[:, np.newaxis], data_transformed))
    data_transformed = pd.DataFrame (data_transformed, columns = ['DateTime', 'PCA_1', 'PCA_2'])
    
    data_transformed = season_data (data_transformed, 'DateTime')
    
    chart_data = alt.Chart(data_transformed).mark_circle(opacity = 1).encode(
        alt.X ('PCA_1:Q'),
        alt.Y ('PCA_2:Q'),
        alt.Color ('season:N', scale = alt.Scale (range = ['blue', 'green', 'orange', 'brown'])),
    ).properties(
        width = 1200
    )
    
    return chart_data


def import_kegg_and_create_df (end = 51, path_fasta = path_genomics_78, path_all_keggs = path_genomics_kegg):
    
    print ('Importing KEGG data')
    
    # There are 51 files for each day, in which there are KEGG IDs for each genome collected that day
    # I have to traverse every KEGG file, and create DataFrame for each and every one 
    
    kegg_files = [i for i in os.listdir(path_all_keggs) if (i.endswith('besthits') and i.startswith('D'))]
    rmags_78_names = [os.path.splitext(i)[0] for i in os.listdir(path_fasta) if (i.endswith('fa') and i.startswith('D'))]
    kegg_data_list = []
    
    # This was done so that I could work with first 100 files only. Otherwise, I should just remove: i, and enumerate
    for i, kegg_file_name in enumerate(kegg_files):
        
        if i == end:
            break
        
        else:
            # Now I create a DataFrame out of it and save it in the list of DataFrames
            tmp_df = pd.read_csv (os.path.join(path_genomics_kegg, kegg_file_name), delimiter = '\t')
            tmp_filter = tmp_df['Gene'].apply(lambda x: str(x).split('_')[0] + '_' + str(x).split('_')[1]).isin(rmags_78_names)
            tmp_df = tmp_df[tmp_filter]
                    
            tmp_df['Gene'] = tmp_df['Gene'].apply(lambda x: str(x).split('_')[0] + '_' + str(x).split('_')[1])
            tmp_df['ID'] = tmp_df['ID'].apply(lambda x: str(x).split(':')[1])
            tmp_df.drop (['maxScore', 'hitNumber'], axis = 1, inplace = True)
            tmp_df.reset_index(drop=True, inplace=True)
            
            kegg_data_list.append (tmp_df)
            
    
    print ('Finished importing')
    return create_kegg_matrix (kegg_data_list, path_fasta)


def create_kegg_matrix (list_data, path_fasta = path_genomics_78):
    
    print ('Creating KEGG matrix')
    
    rmags_78_names = [os.path.splitext(i)[0] for i in os.listdir(path_fasta) if (i.endswith('fa') and i.startswith('D'))]
    result_matrix_df = pd.DataFrame (columns = rmags_78_names)
    
    for i in list_data:
        tmp_df = i.value_counts().reset_index()
        
        for i, row in tmp_df.iterrows():
            result_matrix_df.at[row['ID'], row['Gene']] = row[0]
    
    result_matrix_df.fillna(0, inplace = True)
    
    print ('Finished creating')
    return result_matrix_df.T


def create_pairwise_jaccard (data):
    
    tmp_data = data.clip(0, 1)
    result = squareform(pdist(tmp_data.astype(bool), jaccard))
    
    return pd.DataFrame (result, index = data.index, columns = data.index)


def visualize_with_mds (data, start_date, end, path_fasta = path_genomics_78):
    
    mds_model = MDS(n_components = 2, random_state = SEED, dissimilarity = "precomputed", n_jobs = NUM_OF_WORKERS)
    mds_pos = mds_model.fit_transform(data)
    
    list_of_days = [i for i in os.listdir(path_fasta) if (i.endswith('fa') and i.startswith('D'))]
    temporal_column = create_temporal_column (list_of_days, start_date, end)
    
    data_transformed = pd.DataFrame (mds_pos)
    data_transformed.columns = ['MDS_1', 'MDS_2']
    data_transformed = np.hstack(((np.asarray(temporal_column))[:, np.newaxis], data_transformed))
    data_transformed = pd.DataFrame (data_transformed, columns = ['DateTime', 'MDS_1', 'MDS_2'])
    
    data_transformed = season_data (data_transformed, 'DateTime')
    
    chart_data = alt.Chart(data_transformed).mark_circle(opacity = 1).encode(
        alt.X ('MDS_1:Q'),
        alt.Y ('MDS_2:Q'),
        alt.Color ('season:N', scale = alt.Scale (range = ['blue', 'green', 'orange', 'brown'])),
    )
    
    return chart_data


# ## MAG examination

# ### KEGG examination

# In[6]:


kegg_matrix = import_kegg_and_create_df (end = ALL_DAYS, path_fasta = path_genomics_78, path_all_keggs = path_genomics_kegg)
kegg_matrix


# In[7]:


mag_scaler = preprocessing.StandardScaler()
scaled_keggs_df = mag_scaler.fit_transform(kegg_matrix)
#scaled_keggs_df = kegg_matrix.clip(0, 1)
scaled_keggs_df


# In[8]:


k_range_end = int(math.sqrt(num_of_mags)) # Usually it is sqrt(# of mags)

k_range = range(1, k_range_end)

k_mean_models = [KMeans (n_clusters = i, random_state = SEED) for i in k_range]
k_scores = [k_mean_model.fit(scaled_keggs_df).score(scaled_keggs_df) for k_mean_model in k_mean_models]
k_data = pd.DataFrame ({'k_range':k_range, 'k_scores':k_scores})


# In[9]:


k_num_chart = alt.Chart(data = k_data).mark_line().encode(
    alt.X ('k_range:Q'),
    alt.Y ('k_scores:Q')
)

k_num_chart


# In[10]:


# We can see from the chart above that 6 or 7 clusters are optimal for this task (where END = 25 MAGs)
num_of_clusters = 4

k_means_model = KMeans (n_clusters = num_of_clusters, random_state = SEED)
k_means_predicted = k_means_model.fit_predict(scaled_keggs_df)
k_means_predicted


# In[11]:


k_means_chart = visualize_with_pca (scaled_keggs_df, k_means_predicted, k_means_model.cluster_centers_)
k_means_chart


# ### KEGG examination but with pairwise Jaccard distance matrix (as seen in paper)

# In[43]:


kegg_pairwise = create_pairwise_jaccard (kegg_matrix)
kegg_pairwise


# In[46]:


kegg_mds_chart = visualize_with_mds(kegg_pairwise, START_DATE, END, path_genomics_78)
kegg_mds_chart


# ---
# # VAZNO:
# Sledece sto treba da se uradi je da se nadje transcriptomic data set i da se obradi i on u potpunosti. Nakon toga, treba da se sve podeli po skriptama i da se odluci o dizajnu. Posle ostaje jos da se napravi front end.
# 
# ---

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[12]:


# FOR CLUSTERING I SHOULD CREATE A DATAFRAME WITH MAGs INDEXES AND THEIR VECTOR REPRESENTATIONS
final_model, fasta_names, fasta_ids = import_mags_and_build_model (end = END, path_fasta = path_genomics_78)


# In[13]:


# Train model. It tooks ~10 minutes for END = 25 amount of MAGs
final_model = train_model (final_model, epochs = EPOCHS, end = END)


# In[14]:


final_model.wv.save_word2vec_format(os.path.join (path_model_save_root, 'model_78.bin'), binary=True) 


# Now I should vectorize documents with this model. For further use, I could save this model's weights, and use it to vectorize all mags. That would take a lot, but every MAG will have its vector representation
# > This could be done by importing one MAG at a time, then tokenizing it (like before), then getting vector representations of that MAG's sentences (genomes) and then finding the vector representation of the whole MAG (document). If I do that for one MAG at a time, There is no need to worry about memory
# 

# In[15]:


list_of_mag_vectors = vectorize_mags (final_model, path_fasta = path_genomics_78, end = END)
list_of_mag_vectors[0:2]


# In[16]:


mags_df = pd.DataFrame (list_of_mag_vectors)
mags_df


# ## Data preprocessing

# In[17]:


mag_scaler = preprocessing.StandardScaler()
scaled_mags_df = mag_scaler.fit_transform(mags_df)
scaled_mags_df


# ## Clustering

# ### 1. K-means

# In[18]:


k_range_end = int(math.sqrt(num_of_mags)) # Usually it is sqrt(# of mags)

k_range = range(1, k_range_end)

k_mean_models = [KMeans (n_clusters = i, random_state = SEED) for i in k_range]
k_scores = [k_mean_model.fit(scaled_mags_df).score(scaled_mags_df) for k_mean_model in k_mean_models]
k_data = pd.DataFrame ({'k_range':k_range, 'k_scores':k_scores})


# In[19]:


k_num_chart = alt.Chart(data = k_data).mark_line().encode(
    alt.X ('k_range:Q'),
    alt.Y ('k_scores:Q')
)

k_num_chart


# In[20]:


# We can see from the chart above that 6 or 7 clusters are optimal for this task (where END = 25 MAGs)
num_of_clusters = 4

k_means_model = KMeans (n_clusters = num_of_clusters, random_state = SEED)
k_means_predicted = k_means_model.fit_predict(scaled_mags_df)
k_means_predicted


# In[21]:


k_means_chart = visualize_with_pca (scaled_mags_df, k_means_predicted, k_means_model.cluster_centers_)
k_means_chart


# ### 2. OPTICS

# In[22]:


MIN_SAMPLES = 4

optics_model = OPTICS (min_samples = MIN_SAMPLES, n_jobs = NUM_OF_WORKERS)
optics_predicted = optics_model.fit_predict (scaled_mags_df)
optics_predicted


# In[23]:


# Visualize clusters, since there are no centroids, we are sending bogus array
optics_chart = visualize_with_pca (scaled_mags_df, optics_predicted, np.empty([optics_predicted.shape[0], 1], dtype=int))
optics_chart


# In[24]:


# Side by side comparison
cluster_comparison_chart = alt.hconcat (k_means_chart, optics_chart).resolve_scale(color='independent')
cluster_comparison_chart


# ## Evaluation

# In[25]:


eval_k_means = metrics.silhouette_score (scaled_mags_df, k_means_predicted)
eval_optics = metrics.silhouette_score (scaled_mags_df, optics_predicted)

print ('Silhouette scores: [best = 1, worst = -1]')
print ('\t1. K-means:', eval_k_means)
print ('\t2. OPTICS:', eval_optics)


# ## Visualizing rMAGs with time axis

# In[26]:


time_chart = visualize_temporal_mags (scaled_mags_df, fasta_names, START_DATE, END)
time_chart


# In[27]:


save_charts ([k_means_chart, optics_chart, cluster_comparison_chart, time_chart], ['genomics_k_means_chart.png', 'genomics_optics_chart.png', 'genomics_cluster_comparison_chart.png', 'genomics_time_chart.png'])


# ---
# # METABOLOMIC ANALYSIS
# ---
# ## Importing Metabolomic data

# In[5]:


metabolomics_file_name = os.path.join(path_normalised_metabolomics, os.listdir(path_normalised_metabolomics)[0])
metabolomics_df = pd.read_csv (metabolomics_file_name, delimiter = '\t')
metabolomics_df


# ## Data preprocessing

# In[6]:


metabolomics_df['date'] = pd.to_datetime(metabolomics_df['date'])
metabolomics_df.insert (0, 'date', metabolomics_df.pop('date'))
metabolomics_df.sort_values ('date', inplace = True, ignore_index = True)
metabolomics_df


# In[7]:


# Changing metabolite name if it is unknown
metabolomics_df.loc[metabolomics_df['known_type'].eq('unknown'), 'Metabolite'] = np.nan
metabolomics_df


# In[8]:


print ('Dataset uniqueness:')
print ('\t1. Timestamps:', len(metabolomics_df['date'].unique()))
print ('\t2. Metabolites:', len(metabolomics_df['Metabolite'].unique()))
print ('\t3. Types:', len(metabolomics_df['type'].unique()))
print ('\t4. Known types:', len(metabolomics_df['known_type'].unique()))
print ('\t5. Ns:', len(metabolomics_df['N'].unique()))
print ('\t6. Type 2s:', len(metabolomics_df['type2'].unique()))
print ('\t7. Measurements:', len(metabolomics_df['measurement'].unique()))


# In[9]:


# Saving the name column and removing unnecessairy columns
#metabolite_names = metabolomics_df['Metabolite']
#metabolomics_df.drop(labels = ['Metabolite', 'tp', 'KEGG.Compound.ID', 'Chebi.Name', 'Chebi.Name_combined'], axis = 1, inplace = True)
metabolomics_df.drop(labels = ['tp', 'KEGG.Compound.ID', 'Chebi.Name', 'Chebi.Name_combined'], axis = 1, inplace = True)
metabolomics_df


# In[10]:


# Dummy eencoding categorical data
scaled_metabolomics_df = pd.get_dummies(metabolomics_df, columns = ['type', 'known_type', 'N', 'type2', 'measurement'])
scaled_metabolomics_df


# In[11]:


# Standardizing data
metabolomics_scaler = preprocessing.StandardScaler()
scaled_metabolomics_df[['means', 'medians', 'sds', 'se', 'ci']] = metabolomics_scaler.fit_transform(metabolomics_df[['means', 'medians', 'sds', 'se', 'ci']])

scaled_metabolomics_df


# In[12]:


metabolomics_df.dropna(inplace = True)
metabolomics_df.reset_index(drop=True, inplace=True)
metabolomics_df


# ## Time series examination

# In[13]:


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


# In[14]:


metabolites_chart = visualize_metabolites(metabolomics_df, 'date', 'Metabolite', ['type', 'type2', 'measurement', 'N'])
metabolites_chart


# In[18]:


save_charts ([metabolites_chart], ['metabolomics_metabolites_chart.png'])


# ## Clustering

# In[19]:


# Deep learning temporal clustering

# Should I even do this? Previous visualizations are descriptive enough. It would be a lot of work for not much benefit


# ---
# # PROTEOMIC ANALYSIS
# ---
# ## Importing Proteomic data

# In[20]:


# I could create something similar to Fig. 5 of the original paper, where I would calculate mean of different proteomic feature values for each rMAG calculated by days
# So I would have a table: date | feature 1 | feature 2 | ...
# Where each feature is mean of all values for one day of each MAG in that rMAG


# In[21]:


def import_proteomics (end = 25, path_proteomics = path_proteomics_78):
    
    print ('Importing proteomics data')
    
    # There are 78 FASTA files
    # I have to traverse every FASTA file, and in each file every protein sequence
    
    fasta_files = [i for i in os.listdir(path_proteomics) if (i[-3:] == 'faa')]
    tmp_all = []
    
    # This was done so that I could work with first 100 FASTA files only. Otherwise, I should just remove: i, and enumerate
    for i, fasta_file_name in enumerate(fasta_files):
        
        if i == end:
            break
        
        else:
            with open(os.path.join(path_proteomics, fasta_file_name), 'r') as input_file:
                
                one_mag_list = []
                for fasta_string in SeqIO.parse(input_file, "fasta"):
                    
                    # Analyzing protein (peptide) and creating list of values for one MAG
                    sequence = str(fasta_string.seq)
                    
                    if '*' in sequence:
                        continue
                    
                    else:
                    
                        sequence_analysis = ProteinAnalysis (sequence)
                        
                        tmp_list = [sequence_analysis.molecular_weight(), sequence_analysis.gravy(), sequence_analysis.aromaticity(), sequence_analysis.instability_index(), sequence_analysis.isoelectric_point()]
                        
                        tmp_sec_str = sequence_analysis.secondary_structure_fraction()
                        tmp_list += [tmp_sec_str[0], tmp_sec_str[1], tmp_sec_str[2]]
                        tmp_list.append (sequence.count('K') + sequence.count('R') - sequence.count('D') - sequence.count('E')) # Electricity
                        
                        amino_acid_perc = sequence_analysis.get_amino_acids_percent()
                        
                        tmp_list.append (sum ([amino_acid_perc[aa] for aa in 'AGILPV']))
                        tmp_list.append (sum ([amino_acid_perc[aa] for aa in 'STNQ']))
                        tmp_list.append (sum ([amino_acid_perc[aa] for aa in 'QNHSTYCMW']))
                        tmp_list.append (sum ([amino_acid_perc[aa] for aa in 'AGILPVF']))
                        tmp_list.append (sum ([amino_acid_perc[aa] for aa in 'HKR']))
                        tmp_list.append (sum ([amino_acid_perc[aa] for aa in 'CM']))
                        tmp_list.append (sum ([amino_acid_perc[aa] for aa in 'DE']))
                        tmp_list.append (sum ([amino_acid_perc[aa] for aa in 'NQ']))
                        tmp_list.append (sum ([amino_acid_perc[aa] for aa in 'ST']))
                        
                        # Now I put all these values in one_mag_list as a numpy arrays
                        one_mag_list.append(np.asarray(tmp_list))
                        
                # Now I put one mag values, aggregated by mean, into the all mag list
                tmp_all.append (np.asarray(one_mag_list).mean (axis = 0))
    
    
    COLUMN_LIST = ['Molecular weight', 'Gravy', 'Aromaticity', 'Instability index', 'Isoelectric point', 'Secondary structure fraction 0', 'Secondary structure fraction 1', 'Secondary structure fraction 2', 'Electricity', 'Fraction aliphatic', 'Fraction uncharged polar', 'Fraction polar', 'Fraction hydrophobic', 'Fraction positive', 'Fraction sulfur', 'Fraction negative', 'Fraction amide', 'Fraction alcohol']
    all_mag_df = pd.DataFrame (tmp_all, columns = COLUMN_LIST)
    
    print ('Finished importing')
    
    return all_mag_df

def visualize_proteomics (data):
    
    # Adding another column that replaces temporal data for now
    if 'Index_tmp' not in data.columns:
        data.insert (0, 'Index_tmp', data.index.values)
    
    # Create repeated chart
    chart = alt.Chart(data).mark_area().encode(
        alt.X ('Index_tmp', type = 'quantitative'),
        alt.Y (alt.repeat('row'), type = 'quantitative'),
    ).properties(
        width = 1200
    ).repeat(
        row = data.columns.values
    )#.resolve_scale(color = 'independent', size = 'independent')#.interactive()
    
    return chart


# In[22]:


proteomics_data = import_proteomics (end = num_of_proteomics)
proteomics_data


# In[23]:


chart_proteomics = visualize_proteomics(proteomics_data)
chart_proteomics


# In[24]:


save_charts ([chart_proteomics], ['proteomics_chart_proteomics.png'])


# ---
# # PHYSICO-CHEMICAL ANALYSIS
# ---
# ## Importing Physico-chemical data

# In[25]:


phy_che_file_name = os.path.join(path_physico_chemical, [i for i in os.listdir(path_physico_chemical) if (i.endswith(('.tsv', '.csv')))][1])
phy_che_df = pd.read_csv (phy_che_file_name, decimal = ',')
phy_che_df


# ## Data preprocessing

# In[26]:


phy_che_df.drop(index = 0, axis = 1, inplace = True)
phy_che_df['Date'] = pd.to_datetime(phy_che_df['Date'])
phy_che_df['Time'] = pd.to_timedelta(phy_che_df["Time"], unit = 'h')
phy_che_df


# In[27]:


filtered_phy_che_df = phy_che_df[(phy_che_df['Date'] >= '2011-03-21') & (phy_che_df['Date'] <= '2012-05-03')]
tmp_column = pd.Series(filtered_phy_che_df['Date'] + filtered_phy_che_df['Time'])

filtered_phy_che_df.drop (['Date', 'Time'], axis = 1, inplace = True)
filtered_phy_che_df.reset_index(inplace = True, drop = True)
filtered_phy_che_df = filtered_phy_che_df.apply(lambda x: pd.to_numeric(x.astype(str).str.replace(',','.')))#, errors='coerce'))
filtered_phy_che_df.insert (0, 'DateTime', tmp_column.values)
filtered_phy_che_df


# In[28]:


# Visualize temperature, air_temperature, conductivity, inflow_pH, nitrate, oxygen, pH

def visualize_phy_che (data, temporal_column, list_of_columns):
    
    # Create repeated chart
    chart = alt.Chart(data).mark_line().encode(
        alt.X (temporal_column, type = 'temporal'),#, timeUnit = 'month'),
        alt.Y (alt.repeat('row'), type = 'quantitative'),
    ).properties(
        width = 1200
    ).repeat(
        row = list_of_columns
    )#.resolve_scale(color = 'independent', size = 'independent')#.interactive()
    
    return chart

def visualize_phy_che_heatmap (data):
    
    new_data = data.drop('DateTime', axis = 1)
    corr = new_data.corr().reset_index().melt('index')
    corr.columns = ['var_1', 'var_2', 'correlation']
    
    # Create correlation chart
    chart = alt.Chart(corr).mark_rect().encode(
        alt.X ('var_1', title = None, axis = alt.Axis(labelAngle = -45)),
        alt.Y ('var_2', title = None),
        alt.Color('correlation', legend=None, scale = alt.Scale(scheme='redblue', reverse = True)),
    ).properties(
        width = alt.Step(40),
        height = alt.Step(40)
    )
    
    chart += chart.mark_text(size = 12).encode(
        alt.Text ('correlation', format=".2f"),
        color = alt.condition("abs(datum.correlation) > 0.5", alt.value('white'), alt.value('black'))
    )
    
    return chart.transform_filter("datum.var_1 < datum.var_2") # This returns only lower triangle


# In[29]:


chart_phy_che = visualize_phy_che (filtered_phy_che_df, 'DateTime', filtered_phy_che_df.columns.values[4:])
chart_phy_che_corr = visualize_phy_che_heatmap (filtered_phy_che_df)
chart_phy_che_corr


# In[30]:


chart_phy_che


# In[31]:


save_charts ([chart_phy_che_corr, chart_phy_che], ['physico_chemical_chart_psy_che_corr.png', 'physico_chemical_chart_psy_che.png'])


# In[ ]:





# ---
# # Front-end
# ---

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




