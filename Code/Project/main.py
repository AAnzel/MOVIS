#!/usr/bin/env python
# coding: utf-8

import all_global

# ---
# # GENOMIC ANALYSIS
# ---
def calc_genomics:
    # **Important**: I should review the way I look at MAGs. The names of all fasta files beggining with 'D_##' represent the days those MAGs were obtained. Therefore, I should look at this also as timeseries data. Also, maybe I should only consider 78 MAGs, and not all ~1300.
    # After some consideration, I conclude that I should definetly use only 78 MAGs, because that way I wouldn't be tied to meta-omics data only. I also thinked about what should I visualize in that case. One idea is that I should also encode those MAGs with word2wec, and then make a 3D chart where one dimension is time, and other two dimensions would be PCA dimensions of those MAGs. I could also use this function to visualize proteomics data if I want.
    # 
    # Another important thing is that I should actually include FASTA headers and possibly use only them. That way, I could make figures like in a relevant paper where MAGs are groupped according to their taxonomy etc. I should look more into this.
    # ## MAG examination
    # ### KEGG examination

    kegg_matrix = import_kegg_and_create_df (end = ALL_DAYS, path_fasta = path_genomics_78, path_all_keggs = path_genomics_kegg)

    mag_scaler = preprocessing.StandardScaler()
    scaled_keggs_df = mag_scaler.fit_transform(kegg_matrix)
    #scaled_keggs_df = kegg_matrix.clip(0, 1)

    k_range_end = int(math.sqrt(num_of_mags)) # Usually it is sqrt(# of mags)
    k_range = range(1, k_range_end)

    k_mean_models = [KMeans (n_clusters = i, random_state = SEED) for i in k_range]
    k_scores = [k_mean_model.fit(scaled_keggs_df).score(scaled_keggs_df) for k_mean_model in k_mean_models]
    k_data = pd.DataFrame ({'k_range':k_range, 'k_scores':k_scores})
    
    k_num_chart = alt.Chart(data = k_data).mark_line().encode(
        alt.X ('k_range:Q'),
        alt.Y ('k_scores:Q')
    )

    # We can see from the chart above that 6 or 7 clusters are optimal for this task (where END = 25 MAGs)
    num_of_clusters = 4

    k_means_model = KMeans (n_clusters = num_of_clusters, random_state = SEED)
    k_means_predicted = k_means_model.fit_predict(scaled_keggs_df)

    k_means_chart = visualize_with_pca (scaled_keggs_df, k_means_predicted, k_means_model.cluster_centers_)

    # ### KEGG examination but with pairwise Jaccard distance matrix (as seen in paper)
    kegg_pairwise = create_pairwise_jaccard (kegg_matrix)
    kegg_mds_chart = visualize_with_mds(kegg_pairwise, START_DATE, END, path_genomics_78)

    # ---
    # # VAZNO:
    # Sledece sto treba da se uradi je da se nadje transcriptomic data set i da se obradi i on u potpunosti. Nakon toga, treba da se sve podeli po skriptama i da se odluci o dizajnu. Posle ostaje jos da se napravi front end.
    # 
    # ---
    
    # FOR CLUSTERING I SHOULD CREATE A DATAFRAME WITH MAGs INDEXES AND THEIR VECTOR REPRESENTATIONS
    final_model, fasta_names, fasta_ids = import_mags_and_build_model (end = END, path_fasta = path_genomics_78)
    
    # Train model. It tooks ~10 minutes for END = 25 amount of MAGs
    final_model = train_model (final_model, epochs = EPOCHS, end = END)
    
    final_model.wv.save_word2vec_format(os.path.join (path_model_save_root, 'model_78.bin'), binary=True) 

    # Now I should vectorize documents with this model. For further use, I could save this model's weights, and use it to vectorize all mags. That would take a lot, but every MAG will have its vector representation
    # > This could be done by importing one MAG at a time, then tokenizing it (like before), then getting vector representations of that MAG's sentences (genomes) and then finding the vector representation of the whole MAG (document). If I do that for one MAG at a time, There is no need to worry about memory
    # 
    
    list_of_mag_vectors = vectorize_mags (final_model, path_fasta = path_genomics_78, end = END)
    mags_df = pd.DataFrame (list_of_mag_vectors)
    
    # ## Data preprocessing
    mag_scaler = preprocessing.StandardScaler()
    scaled_mags_df = mag_scaler.fit_transform(mags_df)

    # ## Clustering

    # ### 1. K-means
    k_range_end = int(math.sqrt(num_of_mags)) # Usually it is sqrt(# of mags)
    k_range = range(1, k_range_end)

    k_mean_models = [KMeans (n_clusters = i, random_state = SEED) for i in k_range]
    k_scores = [k_mean_model.fit(scaled_mags_df).score(scaled_mags_df) for k_mean_model in k_mean_models]
    k_data = pd.DataFrame ({'k_range':k_range, 'k_scores':k_scores})

    k_num_chart = alt.Chart(data = k_data).mark_line().encode(
        alt.X ('k_range:Q'),
        alt.Y ('k_scores:Q')
    )

    # We can see from the chart above that 6 or 7 clusters are optimal for this task (where END = 25 MAGs)
    num_of_clusters = 4

    k_means_model = KMeans (n_clusters = num_of_clusters, random_state = SEED)
    k_means_predicted = k_means_model.fit_predict(scaled_mags_df)

    k_means_chart = visualize_with_pca (scaled_mags_df, k_means_predicted, k_means_model.cluster_centers_)

    # ### 2. OPTICS

    MIN_SAMPLES = 4

    optics_model = OPTICS (min_samples = MIN_SAMPLES, n_jobs = NUM_OF_WORKERS)
    optics_predicted = optics_model.fit_predict (scaled_mags_df)

    # Visualize clusters, since there are no centroids, we are sending bogus array
    optics_chart = visualize_with_pca (scaled_mags_df, optics_predicted, np.empty([optics_predicted.shape[0], 1], dtype=int))

    # Side by side comparison
    cluster_comparison_chart = alt.hconcat (k_means_chart, optics_chart).resolve_scale(color='independent')

    # ## Evaluation
    eval_k_means = metrics.silhouette_score (scaled_mags_df, k_means_predicted)
    eval_optics = metrics.silhouette_score (scaled_mags_df, optics_predicted)

    print ('Silhouette scores: [best = 1, worst = -1]')
    print ('\t1. K-means:', eval_k_means)
    print ('\t2. OPTICS:', eval_optics)

    # ## Visualizing rMAGs with time axis
    time_chart = visualize_temporal_mags (scaled_mags_df, fasta_names, START_DATE, END)
    
    save_charts ([k_means_chart, optics_chart, cluster_comparison_chart, time_chart], ['genomics_k_means_chart.png', 'genomics_optics_chart.png', 'genomics_cluster_comparison_chart.png', 'genomics_time_chart.png'])

# ---
# # METABOLOMIC ANALYSIS
# ---
def calc_metabolomics:
    # ## Importing Metabolomic data
    metabolomics_file_name = os.path.join(path_normalised_metabolomics, os.listdir(path_normalised_metabolomics)[0])
    metabolomics_df = pd.read_csv (metabolomics_file_name, delimiter = '\t')

    # ## Data preprocessing

    metabolomics_df['date'] = pd.to_datetime(metabolomics_df['date'])
    metabolomics_df.insert (0, 'date', metabolomics_df.pop('date'))
    metabolomics_df.sort_values ('date', inplace = True, ignore_index = True)

    # Changing metabolite name if it is unknown
    metabolomics_df.loc[metabolomics_df['known_type'].eq('unknown'), 'Metabolite'] = np.nan

    print ('Dataset uniqueness:')
    print ('\t1. Timestamps:', len(metabolomics_df['date'].unique()))
    print ('\t2. Metabolites:', len(metabolomics_df['Metabolite'].unique()))
    print ('\t3. Types:', len(metabolomics_df['type'].unique()))
    print ('\t4. Known types:', len(metabolomics_df['known_type'].unique()))
    print ('\t5. Ns:', len(metabolomics_df['N'].unique()))
    print ('\t6. Type 2s:', len(metabolomics_df['type2'].unique()))
    print ('\t7. Measurements:', len(metabolomics_df['measurement'].unique()))

    # Saving the name column and removing unnecessairy columns
    #metabolite_names = metabolomics_df['Metabolite']
    #metabolomics_df.drop(labels = ['Metabolite', 'tp', 'KEGG.Compound.ID', 'Chebi.Name', 'Chebi.Name_combined'], axis = 1, inplace = True)
    metabolomics_df.drop(labels = ['tp', 'KEGG.Compound.ID', 'Chebi.Name', 'Chebi.Name_combined'], axis = 1, inplace = True)

    # Dummy eencoding categorical data
    scaled_metabolomics_df = pd.get_dummies(metabolomics_df, columns = ['type', 'known_type', 'N', 'type2', 'measurement'])

    # Standardizing data
    metabolomics_scaler = preprocessing.StandardScaler()
    scaled_metabolomics_df[['means', 'medians', 'sds', 'se', 'ci']] = metabolomics_scaler.fit_transform(metabolomics_df[['means', 'medians', 'sds', 'se', 'ci']])

    metabolomics_df.dropna(inplace = True)
    metabolomics_df.reset_index(drop=True, inplace=True)

    # ## Time series examination
    metabolites_chart = visualize_metabolites(metabolomics_df, 'date', 'Metabolite', ['type', 'type2', 'measurement', 'N'])

    save_charts ([metabolites_chart], ['metabolomics_metabolites_chart.png'])
    
    # ## Clustering
    # Deep learning temporal clustering
    # Should I even do this? Previous visualizations are descriptive enough. It would be a lot of work for not much benefit


# ---
# # PROTEOMIC ANALYSIS
# ---
def calc_proteomics:
    # ## Importing Proteomic data
    
    # I could create something similar to Fig. 5 of the original paper, where I would calculate mean of different proteomic feature values for each rMAG calculated by days
    # So I would have a table: date | feature 1 | feature 2 | ...
    # Where each feature is mean of all values for one day of each MAG in that rMAG
    
    proteomics_data = import_proteomics (end = num_of_proteomics)
    chart_proteomics = visualize_proteomics(proteomics_data)

    save_charts ([chart_proteomics], ['proteomics_chart_proteomics.png'])

# ---
# # PHYSICO-CHEMICAL ANALYSIS
# ---
def calc_phy_che:
    # ## Importing Physico-chemical data
    phy_che_file_name = os.path.join(path_physico_chemical, [i for i in os.listdir(path_physico_chemical) if (i.endswith(('.tsv', '.csv')))][1])
    phy_che_df = pd.read_csv (phy_che_file_name, decimal = ',')
    
    # ## Data preprocessing
    phy_che_df.drop(index = 0, axis = 1, inplace = True)
    phy_che_df['Date'] = pd.to_datetime(phy_che_df['Date'])
    phy_che_df['Time'] = pd.to_timedelta(phy_che_df["Time"], unit = 'h')
    
    filtered_phy_che_df = phy_che_df[(phy_che_df['Date'] >= '2011-03-21') & (phy_che_df['Date'] <= '2012-05-03')]
    tmp_column = pd.Series(filtered_phy_che_df['Date'] + filtered_phy_che_df['Time'])

    filtered_phy_che_df.drop (['Date', 'Time'], axis = 1, inplace = True)
    filtered_phy_che_df.reset_index(inplace = True, drop = True)
    filtered_phy_che_df = filtered_phy_che_df.apply(lambda x: pd.to_numeric(x.astype(str).str.replace(',','.')))#, errors='coerce'))
    filtered_phy_che_df.insert (0, 'DateTime', tmp_column.values)

    # Visualize temperature, air_temperature, conductivity, inflow_pH, nitrate, oxygen, pH
    chart_phy_che = visualize_phy_che (filtered_phy_che_df, 'DateTime', filtered_phy_che_df.columns.values[4:])
    chart_phy_che_corr = visualize_phy_che_heatmap (filtered_phy_che_df)
    
    save_charts ([chart_phy_che_corr, chart_phy_che], ['physico_chemical_chart_psy_che_corr.png', 'physico_chemical_chart_psy_che.png'])




# FRONT-END
def main():




if __name__ == '__main__':
    main()
