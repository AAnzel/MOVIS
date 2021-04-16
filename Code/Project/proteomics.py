import all_global

# This functions are used for PROTEOMIC

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
