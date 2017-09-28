#ifndef RSQUARE_RSQUARE_H
#define RSQUARE_RSQUARE_H

#include <string>
#include <vector>
#include <fstream>
#include "VcfFileReader.h"
#include "Haplotype.h"



using namespace std;



class RSquare{

public:

    String FileNameValidation, FileNameImputation, FileNameOutput;

    int numMarkers, numSamples;
    vector<string> IndividualName, MarkerID;

    bool DS = true;
    Haplotype Validation, Imputation;
    vector<double> RSquareResult;

//    Dosage DosageInfoValidation, DosageInfoImputation;
//    vector<double> DosageValidation, DosageImputation, RSquareData;


    RSquare()
    {

    };

    ~RSquare()
    {

    };


    bool GetHaplotypeFromVCF();
    double CalculateRSquare_VectorWise(vector<double> Validation, vector<double> Imputation);
    bool CalculateRSquare();
    void printRSquare();
    bool outputRSquare();

    bool Analysis();




};


#endif //RSQUARE_RSQUARE_H
