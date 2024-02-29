import numpy as np
from pysam import VariantFile
from sklearn import decomposition
import pandas as pd
import matplotlib.pyplot as plt

genotypes =[]
samples =[]
Variant_ids=[]
   
def read_vcf(vcf_filename):
 
   with VariantFile (vcf_filename) as vcf_reader:
      counter=0
      for record in vcf_reader:
           counter += 1
           if counter % 100 ==0:
            alleles =[record.samples[x].allele_indices for x in record.samples]
            samples =[sample for sample in record.samples]
            genotypes.append(alleles)
            Variant_ids.append(record.id)
           if counter >= 10000:
              break
       

   return genotypes, samples, Variant_ids 
 
def read_panel(panel_file):
 
 with open(panel_file) as panel_file:
    labels = {}  # {sample_id: population_code}
    for line in panel_file:
     line = line.strip().split('\t')
     labels[line[0]] = line[1]
 return labels

def Matrix (genotypes):
 
  genotypes= np.array(genotypes)
  matrix = np.count_nonzero(genotypes, axis=2).T

  return matrix

def perform_pca(matrix):
  
  pca = decomposition.PCA(n_components=2)
  transformed_matrix = pca.fit_transform(matrix)

  return transformed_matrix
  
def create_dataframe(matrix, Variant_ids, samples, labels):
  
  df = pd.DataFrame(matrix,columns=Variant_ids, index=samples)
  df['Population code'] = df.index.map(labels)

  return df



# Main part of the script
vcf_filename = "ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf"
panel_file = "phase1_integrated_calls.20101123.ALL.panel"

genotypes, samples, Variant_ids = read_vcf(vcf_filename)
labels = read_panel(panel_file)
matrix = Matrix (genotypes)
transformed_matrix= perform_pca(matrix)

# Create DataFrase and save to CSV

df = create_dataframe(matrix, Variant_ids, samples, labels)
df.rename(columns={df.columns[0]: 'Sample'}, inplace=True)
non_snp_colums = ['Population code', 'Sample']
df.to_csv("matrix.csv")

 