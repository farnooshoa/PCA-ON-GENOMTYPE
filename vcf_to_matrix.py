import numpy as np
from pysam import VariantFile

vcf_filename="ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
panel_file= "integrated_call_samples_v3.20130502.ALL.panel"

genotypes =[]
samples =[]
with VariantFile (vcf_filename) as vcf_reader:
    for record in vcf_reader:
        counter = 0
        alleles = [record.samples[x].allele_indices for x in record.samples]
        samples = [sample for sample in record.samples]
        genotypes.append(alleles)
        counter += 1
        if counter >= 100:
            break
np.array(genotypes)

