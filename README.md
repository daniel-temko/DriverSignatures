# DriverSignatures

Code to analyse the effects of mutation processes and selection on relative frequencies of cancer driver mutations.

## Usage

### Preprocessing
DriverSpectrums.R loads mutation data into R, along with its associated metadata.

Preprocessing.R generates summary tables of mutation occurences and mutation signature scores per patient for use in downstream analysis steps.

### Analysis
DriverSignatureAssociations.R contains the code to analyse associations between driver mutation frequencies and mutational signatures across patients within each cancer type. 

DriverSelectionInference.R contains the code to infer relative fitness effects of driver within cancer types.
