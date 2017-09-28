#include <iostream>
#include <getopt.h>
#include "RSquare.h"

int main(int argc, char **argv) {
    int c;
    String validation_file, imputation_file, output_file;
    struct option longopts[] = {
            {"validation",  required_argument,  NULL, 'v'},
            {"imputation",  required_argument,  NULL, 'i'},
            {"output",  required_argument,  NULL, 'o'},
            {0,0,0,0}

    };

    if(argc < 7)
    {
        cout << "Error! Please check input. Paths for validation, imputation and output are all required." << endl;
        return 0;
    }


    RSquare R;

    while ((c = getopt_long(argc, argv, "v:i:o:", longopts, NULL)) != -1){
        switch (c) {
            case 'v':
                // codes needed here to check if optarg is a valid vcf file.

                R.FileNameValidation = optarg;
                printf("-v with value %s\n", optarg);
                break;

            case 'i':
                // codes needed here to check if optarg is a valid vcf file.

                R.FileNameImputation = optarg;
                printf("-i with value %s\n", optarg);
                break;

            case 'o':
                R.FileNameOutput = optarg;
                printf("-o with value %s\n", optarg);
                break;

            case ':':
                fprintf(stderr, "%s: option '-%c' requires an argument\n",
                        argv[0], optopt);
                break;
            case '?':
                fprintf(stderr, "%s: option '-%c' is invalid: ignored\n",
                        argv[0], optopt);
                break;
        }
    }


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