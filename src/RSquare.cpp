#include "RSquare.h"

bool RSquare::GetDosagefromVCFFile(String &VCFFileName_v, String &VCFFileName_i)
{
    DosageInfoValidation.read(VCFFileName_v);
    DosageInfoImputation.read(VCFFileName_i);
    FileNameValidation = DosageInfoValidation.inFileName;
    FileNameImputation = DosageInfoImputation.inFileName;
    DosageValidation = DosageInfoValidation.DosageData;
    DosageImputation = DosageInfoImputation.DosageData;
    numMarkers = DosageInfoValidation.getNumMarkers();
    numSamples = DosageInfoValidation.getNumSamples();
    return 0;
}


double RSquare::CalculateRSquare_VectorWise(vector<double> Validation, vector<double> Imputation)
{
    int n = (int)Validation.size();
//    if (Imputation.size() != n)
//    {
//        cout << "\n Error! \n" << endl;
//        return false;
//    }

    double sumValidation = 0, sumValidation2 = 0, sumImputation = 0, sumImputation2 = 0, sumValidationImputation = 0;
    double EValidation, EImputation, varValidation, varImputation, covariance;

    for (int i = 0; i < n; i++){
        sumValidation += Validation[i];
        sumValidation2 += Validation[i]*Validation[i];
        sumImputation += Imputation[i];
        sumImputation2 += Imputation[i]*Imputation[i];
        sumValidationImputation += Validation[i]*Imputation[i];
    }

    EValidation = sumValidation * 1.0/n;
    EImputation = sumImputation * 1.0/n;
    covariance = sumValidationImputation * 1.0/n - EValidation * EImputation;
    varValidation = sumValidation2 * 1.0/n - EValidation*EValidation;
    varImputation = sumImputation2 * 1.0/n - EImputation*EImputation;

    return covariance/varValidation*covariance/varImputation;
}

bool RSquare::CalculateRSquare()
{

    vector<double>::iterator startV = DosageValidation.begin(), endV = startV + numSamples;
    vector<double>::iterator startI = DosageImputation.begin(), endI = startI + numSamples;

    cout << "\n Calculating R-Square ..." << endl;

    for (int i = 0; i < numMarkers; i++)
    {
        vector<double> tempValildation(startV, endV), tempImputation(startI, endI);


        RSquareData.push_back(CalculateRSquare_VectorWise(tempValildation, tempImputation));
        startV += numSamples; endV += numSamples;
        startI += numSamples; endI += numSamples;

    }
    return false;
}


void RSquare::printRSquare()
{
    for(int i = 0; i < RSquareData.size(); i++)
    {
        cout << RSquareData[i] << endl;
    }
}

