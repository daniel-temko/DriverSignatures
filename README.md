# DriverSignatures

Code to analyse the effects of mutation processes and selection on relative frequencies of cancer driver mutations.

## Usage

### Requirements

In addition to cancer mutation data, the preprocessing scripts require mutation metadata to run:

| File | Description |
| --- | --- |
| Vogelstein.Cancer.Genome.Landscapes.Table.S2a.edit.2.Aug.2017.csv | List of driver genes from ref 16 in our paper, and is also given in Supp. Data 6 |
| knownCanonical.5.June.2017.txt | Table mapping gene ID's to the canonical transcripts for the genes. It can be obtained from the UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables) |
| knownToRefSeq.6.June.tsv | Table giving the refseq ID for each canonical transcript; this can also be obtained from the UCSC Table Browser |
| trinucleotide.frequencies.tsv | Table giving the frequencies of trinucleotide sequences in the human exome and human genome. This can be obtained using the get_context_freq of the R package SigsPack |
| singature.probabilities.txt | Table containing the proportion of each of the 96 mutation types in each mutational signature - it can be obtained here: (https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt) |
| MasterSampleList.csv | Table containing metadata for the different mutation datasets to be analyzed - this needs to be tailored to the specific datasets being studied (see _data_ for example) |
| TCGA.Files.csv | Table contatining file paths for TCGA mutation data | 
| Disease.Classifications.csv | Table mapping cancer types to cancer groups with COSMIC mutational signature annotations |

### Preprocessing
DriverSpectrums.R loads mutation data into R, along with its associated metadata.

Preprocessing.R generates summary tables of mutation occurences and mutation signature scores per patient for use in downstream analysis steps.

### Analysis
DriverSignatureAssociations.R contains the code to analyse associations between driver mutation frequencies and mutational signatures across patients within each cancer type. 

DriverSelectionInference.R contains the code to infer relative fitness effects of driver within cancer types.
