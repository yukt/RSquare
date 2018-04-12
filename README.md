# R-Square Calculator - A Tool for Calculating Imputation Accuracy
## Installation
```
git clone https://github.com/yukt/RSquare.git
cd RSquare
bash install.sh
```

## Usage
```
RSquare -v [Validation.vcf.gz]         [Required] Input Validation File
        -i [Imputation.vcf.gz]         [Required] Input Imputation File
        -o [OutputPrefix]              [Required] Output Prefix
        --validationFormat [GT/DS/GP]  [Optional] Genotype info format (Default: GT)
        --imputationFormat [GT/DS/GP]  [Optional] Genotype info format (Default: DS)
```
