# DriverSignatures

Code to analyse the effects of mutation processes and selection on relative frequencies of cancer driver mutations.

## Usage

### Requirements

In addition to cancer mutation data, the preprocessing scripts require mutation metadata to run:

- _Vogelstein.Cancer.Genome.Landscapes.Table.S2a.edit.2.Aug.2017.csv_: Driver genes from ref 16 in our paper, and is also given in Supp. Data 6
- _knownCanonical.5.June.2017.txt_: Table mapping gene ID's to the canonical transcripts for the genes. It can be obtained from the UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables)
- _knownToRefSeq.6.June.tsv_: Table giving the refseq ID for each canonical transcript; this can also be obtained from the UCSC Table Browser
- _trinucleotide.frequencies.tsv_: Table giving the frequencies of trinucleotide sequences in the human exome and human genome. This can be obtained using the get_context_freq of the R package SigsPack
- _singature.probabilities.txt_: Table containing the proportion of each of the 96 mutation types in each mutational signature - it can be obtained here: (https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt)
- _signatures.and.diseases.csv_: Table giving the cancer groups in which each signature has been found - this information is also available from: (https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt)
- _signature.annotations.csv_: Table containing reported associations for each signature, and identifying mutational signatures reported as 'clock-like' in reference 17 from our study.  
- _MasterSampleList.csv_: Table containing metadata for the different mutation datasets to be analyzed - this needs to be tailored to the specific datasets being studied (see _data_ for example)
- _Disease.Classifications.csv_: Table mapping cancer types to cancer groups with COSMIC mutational signature annotations
- _TCGA.Files.csv_: Table containing file paths for TCGA mutation data

In addition, please note line 396 of DriverSpectrums.R uses the results of running ANNOVAR (1) on the output from line 392 

### Preprocessing
DriverSpectrums.R loads mutation data into R, along with its associated metadata.

Preprocessing.R generates summary tables of mutation occurences and mutation signature scores per patient for use in downstream analysis steps.

### Analysis
DriverSignatureAssociations.R contains the code to analyse associations between driver mutation frequencies and mutational signatures across patients within each cancer type. 

DriverSelectionInference.R contains the code to infer relative fitness effects of driver within cancer types.

Plotting.R contains the code to generate the figures from the paper. Test

### References

1. Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research, 38:e164, 2010
