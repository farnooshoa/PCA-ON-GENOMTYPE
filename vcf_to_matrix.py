import numpy as np
from pysam import VariantFile
from sklearn import decomposition
import pandas as pd
vcf_filename="ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf"
panel_file= "phase1_integrated_calls.20101123.ALL.panel"
def load_genomic_data(vcf_filename, panel_file, variant_sampling_rate=100):
    genotypes = []
    samples = []
    variant_ids = []

    with VariantFile(vcf_filename) as vcf_reader:
        for counter, record in enumerate(vcf_reader, start=1):
            if counter % variant_sampling_rate == 0:
                alleles = [record.samples[x].allele_indices for x in record.samples]
                samples = [sample for sample in record.samples]
                genotypes.append(alleles)
                variant_ids.append(record.id)
            if counter >= variant_sampling_rate:
                break

    with open(panel_file) as panel_file:
         labels = {}  # {sample_id: population_code}
         for line in panel_file:
             line = line.strip().split('\t')
             labels[line[0]] = line[1]       
             
    genotypes = np.array(genotypes)
    print(genotypes.shape)

def perform_pca(genotypes_matrix, n_components=1):
    matrix = np.count_nonzero(genotypes_matrix, axis=2)
    matrix = matrix.T

    pca = decomposition.PCA(n_components=n_components)
    pca.fit(matrix)

    to_plot = pca.transform(matrix)
    return matrix, to_plot, pca.singular_values_    

def create_dataframe(matrix, samples, variant_ids, labels):
    df = pd.DataFrame(matrix, columns=variant_ids, index=samples)
    df['Population code'] = df.index.map(labels)
    return df



def main():
    vcf_filename = "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf"
    panel_file = "phase1_integrated_calls.20101123.ALL.panel"
    
    # Load genomic data
    genotypes_matrix, samples, variant_ids, labels = load_genomic_data(vcf_filename, panel_file)

    # Perform PCA analysis
    matrix, to_plot, singular_values = perform_pca(genotypes_matrix)

    # Create DataFrame and export to CSV
    df = create_dataframe(matrix, samples, variant_ids, labels)
    df.to_csv("matrix.csv")

    # Print PCA results
    print("Original Genotype Matrix Shape:", genotypes_matrix.shape)
    print("Transformed Matrix Shape:", matrix.shape)
    print("PCA Singular Values:", singular_values)
    print("Transformed Data Shape:", to_plot.shape)

if __name__ == "__main__":
    main()