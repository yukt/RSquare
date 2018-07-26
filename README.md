# RSquare - A Tool for Calculating Imputation Accuracy
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
        --AF [AlleleFrequency File]    [Optional] See Note(a) below.
        --bins [Bins File]             [Optional] Default: See defaultBins.txt
Note: (a) AlleleFrequency file must contain two tab-delimited columns with header 'SNP\tAF'
          and SNPs should be same as the validation file.
```

The tool will compute SNP-wise imputation accuracy (.RSquare);

With `--AF` option, RSquare will automatically compute aggregated imputation accuracy (.aggRSquare) in addition; SNPs will be aggregated by Allele Frequency; if no Bins File is provided, it will use the default bins; `--bins` will be ignored when no `--AF` is detected.
