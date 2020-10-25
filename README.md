# DriverSignatures

Code to analyse the effects of mutation processes and selection on relative frequencies of cancer driver mutations.

## Usage

DriverSpectrums.R reads in whole genome and whole exome sequencing data from TCGA and ICGC, and associated metadata.

Preprocessing.R generates summary tables of mutation occurences and mutation signature scores per patient for use in downstream analysis steps.

DriverSignatureAssociations.R contains the code to analyse associations between driver mutation frequencies and mutational signatures across patients within each cancer type. 

DriverSelectionInference.R contains the code to infer relative fitness effects of driver within cancer types.
