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

   return np.array(genotypes), samples, Variant_ids  

def read_panel(panel_file):
   with open(panel_file) as panel_file:
    labels = {}  # {sample_id: population_code}
    for line in panel_file:
        line = line.strip().split('\t')
        labels[line[0]] = line[1]       
    return labels         

def perform_pca(genotypes_matrix):
    matrix = np.count_nonzero(genotypes_matrix, axis=2).T

    pca = decomposition.PCA(n_components=2)
    transformed_matrix = pca.fit_transform(matrix)
    

    return transformed_matrix


def create_dataframe(matrix, Variant_ids, samples, labels):
   df = pd.DataFrame(matrix,columns=Variant_ids, index=samples)
   df['Population code'] = df.index.map(labels)
   return df

def save_to_csv(df, output_filename="matrix.csv"):
    df.to_csv(output_filename)


# Main part of the script
vcf_filename = "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf"
panel_file = "phase1_integrated_calls.20101123.ALL.panel"

genotypes_matrix, samples, variant_ids = read_vcf(vcf_filename)
labels = read_panel(panel_file)

transformed_matrix= perform_pca(genotypes_matrix)

# Create DataFrame and save to CSV
df = create_dataframe(transformed_matrix, variant_ids, samples, labels)
save_to_csv(df)

print(df[100])