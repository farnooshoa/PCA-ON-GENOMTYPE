import numpy as np
from pysam import VariantFile
from sklearn import decomposition
import pandas as pd

vcf_filename="ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf"
panel_file= "phase1_integrated_calls.20101123.ALL.panel"

genotypes =[]
samples =[]
Variant_ids=[]
with VariantFile (vcf_filename) as vcf_reader:
    counter=0
    for record in vcf_reader:
       counter += 1
       if counter % 100 ==0:
           alleles =[record.samples[x].allele_indices for x in record.samples]
           samples =[sample for sample in record.samples]
           genotypes.append(alleles)
           Variant_ids.append(record.id)
       if counter >= 100:
          break
       

with open(panel_file) as panel_file:
    labels = {}  # {sample_id: population_code}
    for line in panel_file:
        line = line.strip().split('\t')
        labels[line[0]] = line[1]       
             
genotypes = np.array(genotypes)
print(genotypes.shape)

matrix = np.count_nonzero(genotypes,axis=2)

matrix = matrix.T
print(matrix.shape)

pca = decomposition.PCA(n_components=1)
pca.fit(matrix)
print(pca.singular_values_)
to_plot = pca.transform(matrix)
print(to_plot.shape)


df = pd.DataFrame(matrix,columns=Variant_ids, index=samples)
df['Population code'] = df.index.map(labels)
df.to_csv("matrix.csv")