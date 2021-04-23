# Ideas for choosing visualization idioms depending on the used data set
## Time-series multi-omics data tool


***What?*** question:

Data set         | Data                                                 
---|---
Genomics         | 1 FASTA file for each timestamp and tabular with KOs 
Transcriptomics  | NaN                                                  
Metabolomics     | Tabular data (.csv, .xlsx...) with time column       
Proteomics       | 1 FASTA file for each timestamp and encoding tabular 
Physico-chemical | Tabular data (.csv, .xlsx...) with time column       

---
***Why?*** question:

Analyze                                        | Search                              | Query      
---|---|---
Consume: discover and present.                 | Target unknown: browse and explore. | Identify
Produce: record and derive. (**without** annotate) | Target known: **nothing** (**without** lookup and locate).      | Compare  
                   NaN                         |   NaN                               | **Without** summarize

---
***How?*** question:

Data set         | Encode                                               | Manipulate                        | Facet                     | Reduce
---|---|---|---|---
Genomics         | **FASTA**: W2V -> PCA -> Map; **KOs**: Map                              | Select and navigate (**without** change) | Juxtapose (**without** partition and superimpose) | Filter (**without** aggregate and embed)
Transcriptomics  | NaN                                                                     | Select and navigate (**without** change) | Juxtapose (**without** partition and superimpose) |Filter (**without** aggregate and embed)
Metabolomics     | Map                                                                     | Select and navigate (**without** change) | Juxtapose (**without** partition and superimpose) |Filter (**without** aggregate and embed)
Proteomics       | **FASTA**: W2V -> PCA -> Map; **Tabular**: BioPython -> Tabular -> Map  | Select and navigate (**without** change) | Juxtapose (**without** partition and superimpose) |Filter (**without** aggregate and embed)
Physico-chemical | Map                                                                     | Select and navigate (**without** change) | Juxtapose (**without** partition and superimpose) |Filter (**without** aggregate and embed)