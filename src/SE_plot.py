import pandas as pd
import matplotlib.pyplot as plt
import h5py
import seaborn as sns
import numpy as np

# stacked percentage bar plot for the number of SE and TE for each sample
df = pd.read_csv("samples_numberOf_SE_TE.txt", sep = "\t")
# convert to TCGA id
meta = pd.read_csv("meta/meta.txt", sep = "\t")
meta = meta.iloc[:,[0,2]]
meta['SampleID'] = meta['SampleID'].str.replace('_C', '')
data = pd.merge(df, meta, left_on='sample', right_on='SampleID', how='inner')
data.shape # (22, 5)
data['SE_percent'] = data['SE']/(data['SE'] + data['TE']) * 100
data['TE_percent'] = data['TE']/(data['SE'] + data['TE']) * 100
data = data.iloc[:, 4:7]
data = data.rename(columns={'SE_percent':'SE', 'TE_percent':'TE'})

# plot (data['SE'] ranges from 2.308278427681413 percent to 9.051652549410397 percent)
fig, ax = plt.subplots(figsize=(8, 6))
data.plot(x = 'TCGAID', kind = 'barh', stacked = True)
plt.ylabel('')
plt.xlabel('% of total enhancers')
#plt.savefig("percent_of_total_enhancers.png", bbox_inches='tight', dpi = 800)
plt.savefig("percent_of_total_enhancers.pdf", bbox_inches='tight', format = 'pdf')

# stacked percentage bar plot for SE and TE signal enrichment for each sample
def getStats(group, sample):
   data = pd.read_csv('superenhancer/SE_samplewise/samples_in_cluster' + str(group) + '/' + sample + '_superenhancer/OVCA_AllStitched.table.txt', skiprows = 5, sep = '\t')
   data['enrichment'] = data.iloc[:,6] - data.iloc[:,7]
   TE_enrichment = data.groupby('isSuper')['enrichment'].sum()[0]
   SE_enrichment = data.groupby('isSuper')['enrichment'].sum()[1]
   return {'sample': sample, 'TE_enrichment': TE_enrichment, 'SE_enrichment': SE_enrichment}

samples = pd.read_csv('samples_numberOf_SE_TE.txt', sep = "\t")
samples['group'] = pd.concat([pd.Series([1]*11), pd.Series([2]*11)], ignore_index=True)
samples = samples.iloc[:,[0,3]]

df = pd.DataFrame(columns=['sample', 'SE_enrichment', 'TE_enrichment'])
for index, row in samples.iterrows():
   g = row['group']
   s = row['sample']
   df = df.append(getStats(g, s), ignore_index = True)

df.shape # (22, 3)
df['SE_percent'] = df['SE_enrichment']/(df['SE_enrichment'] + df['TE_enrichment']) * 100
df['TE_percent'] = df['TE_enrichment']/(df['SE_enrichment'] + df['TE_enrichment']) * 100
df = df.iloc[:, [0,3,4]]
df = df.rename(columns={'SE_percent':'SE', 'TE_percent':'TE'})
# convert to TCGA id
meta = meta/meta.txt", sep = "\t")
meta = meta.iloc[:,[0,2]]
meta['SampleID'] = meta['SampleID'].str.replace('_C', '')
data = pd.merge(df, meta, left_on='sample', right_on='SampleID', how='inner')
data = data.iloc[:,[4,1,2]]

# plot (data['SE'] ranges from 23.90921055793199 percent to 46.97827174671048 percent with mean: 36.94677614766845 percent)
fig, ax = plt.subplots(figsize=(8, 6))
data.plot(x = 'TCGAID', kind = 'barh', stacked = True)
plt.ylabel('')
plt.xlabel('% of total H3K27ac enrichment')
plt.savefig("percent_of_total_H3K27ac_enrichment.pdf", bbox_inches='tight', format = 'pdf')

# unsupervised hierarchical clustering of primary samples using the genomic locations of all SEs
merged_df = pd.read_csv("samples_SE.bed", sep = "\t", header = None)
merged_df = merged_df.iloc[:,3]
merged_df.index = merged_df.values
samples = pd.read_csv("samples_numberOf_SE_TE.txt", sep = "\t")
samples = samples['sample'].values.tolist()

def getData(sample, merged_df):
   with h5py.File("SE_counts_per_sample/" + sample + "/counts.h5", "r") as f:
      region = list(f.keys())[2]
      data = pd.DataFrame(f[region][()])
      data.index = [s.decode('utf-8') for s in data.iloc[:,2]]
      data = data.iloc[:,6] # raw count
      data.name = sample
      merged_df = pd.merge(merged_df, data, left_index = True, right_index = True)
      return(merged_df)

for id in samples:
   merged_df = getData(id, merged_df)

merged_df = merged_df.drop(columns = [3])
merged_df.shape # (16629, 22)
correlation_matrix = np.corrcoef(merged_df, rowvar=False)

# assign color
def assign_color(cluster):
    if cluster == "Cluster_1":
        return "#d62728"
    elif cluster == "Cluster_2":
        return "#1f77b4"
    else:
        return None  # Handle other cases if needed


meta = pd.read_csv("Survival/survival.txt", sep = " ")
meta = meta.iloc[:,[0,6]]
meta['color'] = meta['NMF'].apply(assign_color)
meta['SampleID'] = meta['SampleID'].str.replace('_C', '')
cols = dict(zip(meta['SampleID'], meta['color']))
row_colors = merged_df.columns.map(cols)

obj = sns.clustermap(correlation_matrix, annot=True, row_colors = row_colors, cmap='coolwarm') # get the linkage matrix so that we could reorder the samples while preserving the clustering order
linkage_matrix = obj.dendrogram_row.linkage
linkage_matrix[18][[0,1]] = linkage_matrix[18][[1,0]]
linkage_matrix[19][[0,1]] = linkage_matrix[19][[1,0]]
linkage_matrix[20][[0,1]] = linkage_matrix[20][[1,0]]

# save as pdf, the only cluster2 sample clustered with cluster1 samples is OVCA.3274/"TCGA-23-2649"
fig, ax = plt.subplots(figsize=(8, 6))
g = sns.clustermap(correlation_matrix, row_linkage = linkage_matrix, col_linkage = linkage_matrix, annot = True, row_colors = row_colors, cmap='YlOrBr')
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_yticks([])
g.ax_heatmap.set_xticklabels([])
g.ax_heatmap.set_yticklabels([])
print(g)
plt.savefig("SE_counts.clustermap_color.pdf", format = 'pdf')
