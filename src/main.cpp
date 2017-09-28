#include <iostream>
#include <getopt.h>
#include "RSquare.h"



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
    fprintf(fp, "\n (c) 2017 - Ketian \n");
    fprintf(fp, "\n About:   Write a Short Description ??? \n");
    fprintf(fp,  " URL    : Maybe Make a wiki Page for this ?\n");


    fprintf(fp, "\n Usage:   RSquare -v [Validation.vcf.gz] // Input Validation File\n");
    fprintf(fp, "                  -i [Imputed.vcf.gz]    // Input Imputed Dosage File\n");
    fprintf(fp, "                  -o [RsquareOutput]     // Output Prefix\n");
    fprintf(fp, "\n\n");
    exit(1);

}




int main(int argc, char **argv) {

    if (argc < 2) { usage(stderr); return 1; }
    else if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) { usage(stdout); return 0; }


        int c;
    String validation_file, imputation_file, output_file;

    static struct option loptions[] =
    {

        {"validation",  required_argument,  NULL, 'v'},
        {"imputation",  required_argument,  NULL, 'i'},
        {"output",  required_argument,  NULL, 'o'},
        {NULL,0,NULL,0}
    };

    RSquare R;
    while ((c = getopt_long(argc, argv, "v:i:o:",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'v': R.FileNameValidation = optarg; printf("-v with value %s\n", optarg); break;
            case 'i': R.FileNameImputation = optarg; printf("-i with value %s\n", optarg); break;
            case 'o': R.FileNameOutput = optarg;     printf("-o with value %s\n", optarg); break;
            case '?': usage(stderr); break;
            default: error("[ERROR:] Unknown argument: %s\n", optarg);
        }
    }

    if(R.FileNameValidation =="") { error("[ERROR:] Missing Mandatory argument -v: %s\n");}
    if(R.FileNameImputation =="") { error("[ERROR:] Missing Mandatory argument -i: %s\n");}


//    Dosage DosageData;
//    DosageData.read(validation_file);
//    cout << 7*1.0/2 << endl;

//    double n[] = {1,2,3,4,5}, m[] = {5,3,6,3,1};
//    vector<double> v1(n, n+5), v2(m,m+5);
//
//    RSquare R;
//    cout << "\n" << R.CalculateRSquare_VectorWise(v1,v2) << endl;

//
//    R.GetDosagefromVCFFile(validation_file, imputation_file);
//    R.CalculateRSquare();
//    R.printRSquare();




    R.Analysis();
    return 0;
}