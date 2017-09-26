#ifndef RSQUARE_RSQUARE_H
#define RSQUARE_RSQUARE_H

#include <string>
#include <vector>
#include "VcfFileReader.h"
#include "Dosage.h"


using namespace std;



class RSquare{
    public:

        RSquare()
        {

        };

        ~RSquare()
        {

        };

        String FileNameValidation, FileNameImputation;
        int numMarkers, numSamples;

        Dosage DosageInfoValidation, DosageInfoImputation;
        vector<double> DosageValidation, DosageImputation, RSquareData;


        bool GetDosagefromVCFFile(String &VCFFileName_v, String &VCFFileName_i);
        double CalculateRSquare_VectorWise(vector<double> Validation, vector<double> Imputation);

        bool CalculateRSquare();
        void printRSquare();


};


#endif //RSQUARE_RSQUARE_H
