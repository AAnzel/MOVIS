import all_global

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
