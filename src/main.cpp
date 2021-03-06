#include <iostream>
#include <getopt.h>
#include "FileTypeCheck.h"

void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}


void usage(FILE *fp)
{
    fprintf(fp, "\n");
    fprintf(fp, " ---------------------------------------------------------------------------------------------------- \n");
    fprintf(fp, "               R-Square Calculator - A Tool for Calculating Imputation Accuracy R-Square  \n");
    fprintf(fp, " ----------------------------------------------------------------------------------------------------\n");
    fprintf(fp, "\n (c) 2017 - Ketian Yu \n");
    fprintf(fp, "\n About:   This tool is used for evaluating the imputation accuracy of genotype data using R-Square.\n"
                "          Given Allele Frequency (AF), imputation accuracy aggregated by AF will also be generated.\n");
    fprintf(fp,  " URL  :   https://github.com/yukt/RSquare\n");


    fprintf(fp, "\n Usage:   RSquare -v [Validation.vcf.gz]         // [Required] Input Validation File\n");
    fprintf(fp, "                  -i [Imputation.vcf.gz]         // [Required] Input Imputation File\n");
    fprintf(fp, "                  -o [OutputPrefix]              // [Required] Output Prefix\n");
    fprintf(fp, "                  --validationFormat [GT/DS/GP]  // [Optional] Genotype info format (Default: GT)\n");
    fprintf(fp, "                  --imputationFormat [GT/DS/GP]  // [Optional] Genotype info format (Default: DS)\n");
    fprintf(fp, "                  --AF [AlleleFrequency File]    // [Optional] See Note (a)\n");
    fprintf(fp, "                  --bins [Bins File]             // [Optional] Default: See defaultBins.txt\n");
    fprintf(fp, " Note: (a) AlleleFrequency file must contain two tab-delimited columns with header 'SNP\tAF',\n"
                "           and SNPs should be same as the validation file. \n");
    fprintf(fp, "\n\n");
    exit(1);

}

void version(FILE *fp)
{
    fprintf(fp, "\n");
    fprintf(fp, "[INFO] Program version: 20180725\n");
    fprintf(fp, "\n\n");
}

void finish(FILE *fp)
{
    fprintf(fp, "\n");
    fprintf(fp, "RSquare finished successfully.\n");
    fprintf(fp, "\n\n");
}



int main(int argc, char **argv) {

    if (argc < 2) { usage(stderr); return 1; }
    else if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) { usage(stdout); return 0; }
    else if (strcmp(argv[1], "version") == 0 || strcmp(argv[1], "--version") == 0) { version(stdout); return 0; }

    int c;
    SummaryData R;

    static struct option loptions[] =
    {

        {"validation",       required_argument, NULL, 'v'},
        {"imputation",       required_argument, NULL, 'i'},
        {"output",           required_argument, NULL, 'o'},
        {"validationFormat", required_argument, NULL, 'f'},
        {"imputationFormat", required_argument, NULL, 'g'},
        {"AF",               required_argument, NULL, 'a'},
        {"bins",             required_argument, NULL, 'b'},
        {NULL,0,NULL,0}
    };


    while ((c = getopt_long(argc, argv, "v:i:o:f:g:a:b:",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'v': vcfCheck(c,optarg);       R.FileNameValidation    = optarg;   break;
            case 'i': vcfCheck(c,optarg);       R.FileNameImputation    = optarg;   break;
            case 'o': outputCheck(optarg);      R.OutputPrefix          = optarg;   break;
            case 'f': formatCheck(c, optarg);   R.validationFormat      = optarg;   break;
            case 'g': formatCheck(c, optarg);   R.imputationFormat      = optarg;   break;
            case 'a':                           R.FileNameAlleleFreq    = optarg;   break;
            case 'b':                           R.FileNameBins          = optarg;   break;
            case '?': usage(stderr);                                                break;
            default: error("[ERROR] Unknown argument: %s\n", optarg);
        }
    }

    if(R.FileNameValidation =="") { error("[ERROR] Missing Mandatory argument -v\n");}
    if(R.FileNameImputation =="") { error("[ERROR] Missing Mandatory argument -i\n");}
    if(R.OutputPrefix =="")       { error("[ERROR] Missing Mandatory argument -o\n");}

    version(stdout);
    R.analysis();
    finish(stdout);
    return 0;
}