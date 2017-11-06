#include <iostream>
#include <getopt.h>
#include "FileTypeCheck.h"
#include "SummaryData.h"



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
    fprintf(fp, " -------------------------------------------------------------------------------- \n");
    fprintf(fp, "     R-Square Calculator - A Tool for Calculating Imputation Accuracy R-Square  \n");
    fprintf(fp, " --------------------------------------------------------------------------------\n");
    fprintf(fp, "\n (c) 2017 - Ketian Yu \n");
    fprintf(fp, "\n About:   This tool is used for evaluating the imputation accuracy of genotype data.\n"
                "          [Optional] Given AlleleFrequency, aggregated RSquares will be generated.\n");
    fprintf(fp,  " URL  :   https://github.com/yukt/RSquare\n");


    fprintf(fp, "\n Usage:   RSquare -v [Validation.vcf.gz]         // Input Validation File\n");
    fprintf(fp, "                  -i [Imputation.vcf.gz]         // Input Imputation File\n");
    fprintf(fp, "                  -o [RSquareOutputPrefix]       // Output Prefix\n");
    fprintf(fp, "                  --validationFormat [GT/DS]     // [Optional] Default: GT\n");
    fprintf(fp, "                  --imputationFormat [GT/DS]     // [Optional] Default: DS\n");
    fprintf(fp, "                  --AF [AlleleFrequency File]    // [Optional] Input Allele Frequency File\n");
    fprintf(fp, " Note: AlleleFrequency File should contain exactly the same SNPs with imputation file, \n"
                "       and should begin with header 'CHROM POS REF ALT AF'.\n");
    fprintf(fp, "\n\n");
    exit(1);

}




int main(int argc, char **argv) {

    if (argc < 2) { usage(stderr); return 1; }
    else if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) { usage(stdout); return 0; }

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
        {NULL,0,NULL,0}
    };


    while ((c = getopt_long(argc, argv, "v:i:o:f:g:a:",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'v': vcfCheck(c,optarg);       R.FileNameValidation    = optarg;   break;
            case 'i': vcfCheck(c,optarg);       R.FileNameImputation    = optarg;   break;
            case 'o': outputCheck(optarg);      R.OutputPrefix          = optarg;   break;
            case 'f': formatCheck(c, optarg);   R.validationFormat      = optarg;   break;
            case 'g': formatCheck(c, optarg);   R.imputationFormat      = optarg;   break;
            case 'a': AFCheck(optarg);          R.FileAF                = optarg;   break;
            case '?': usage(stderr);                                                break;
            default: error("[ERROR:] Unknown argument: %s\n", optarg);
        }
    }

    if(R.FileNameValidation =="") { error("[ERROR:] Missing Mandatory argument -v: %s\n");}
    if(R.FileNameImputation =="") { error("[ERROR:] Missing Mandatory argument -i: %s\n");}


    R.analysis();
    return 0;
}