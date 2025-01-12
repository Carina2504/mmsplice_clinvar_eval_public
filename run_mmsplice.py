# Task: create run_mmsplice.py script: Apply MMSplice to the ClinVar variants. (This can take 5-10min.)

## Install dependencies:
# conda create -n mmsplice_env python=3.9 -y
# conda activate mmsplice_env
# conda install cyvcf2 cython -y
# pip install mmsplice


import warnings
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_save, predict_all_table
from mmsplice.utils import max_varEff


# optional: to suppress FutureWarning messages
warnings.filterwarnings("ignore", category=FutureWarning)

# paths to data
gtf = 'data/chr1.gtf.gz' # Ensembl reference annotation
fasta = 'data/chr1.fa' # chr1 reference sequence


def process_variants(variant_type):
    # process the ClinVar variants (benign or pathogenic on chromosome 1 in this case) and save outputs to file
    print(f"Processing {variant_type} variants")

    vcf = f"data/clinvar_chr1_{variant_type}.vcf.gz" # paths to ClinVar variants
    output_file = f"output/{variant_type}_predictions.tsv" # path to output folder to save the results in

    dl = SplicingVCFDataloader(gtf, fasta, vcf, tissue_specific=False)
    model = MMSplice()
    predictions = predict_all_table(model, dl, pathogenicity=True, splicing_efficiency=True)
    predictionsMax = max_varEff(predictions)
    predictionsMax.to_csv(output_file, sep="\t", index=False, header=True)

    print(f"Done processing {variant_type} variants.")


for variant in ["benign", "pathogenic"]:
    process_variants(variant)

print("Done.")