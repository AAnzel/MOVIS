import all_global

# This functions are used for PHYSICO-CHEMICAL

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
