#include <Error.h>
#include "RSquare.h"

bool RSquare::GetHaplotypeFromVCF()
{
    Validation.DS = DS;
    Imputation.DS = DS;
    Validation.read(FileNameValidation);
    Imputation.read(FileNameImputation);

    // strcmp each markerID is need instead of simply compare numMarkers. Same for Samples.
    if (Validation.numMarkers == Imputation.numMarkers)
    {
        numMarkers = Validation.numMarkers;
        MarkerID = Validation.markerID;
//        cout << "\nValidation and Imputation the same Markers." << endl;
    }
    else
    {
        error("[Error:] Number of markers do not match in Validation and Imputation files.");

    }


    if (Validation.numSamples == Imputation.numSamples)
    {
        numSamples = Validation.numSamples;
        IndividualName = Validation.individualName;
//        cout << "\nValidation and Imputation have the same Samples." << endl;
    }
    else
    {
        error("[Error:] Number of samples do not match in Validation and Imputation files.");
    }

    return 0;
}


double RSquare::CalculateRSquare_VectorWise(vector<double> Validation, vector<double> Imputation)
{
    int n = int(Validation.size());

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
    RSquareResult.resize((unsigned long)numMarkers);

    vector<double>::iterator startV = Validation.HaplotypeData.begin(), endV = startV + numSamples;
    vector<double>::iterator startI = Imputation.HaplotypeData.begin(), endI = startI + numSamples;

    cout << "\nCalculating R-Square ..." << endl;

    for (int i = 0; i < numMarkers; i++)
    {
        vector<double> tempValildation(startV, endV), tempImputation(startI, endI);


        RSquareResult[i] = CalculateRSquare_VectorWise(tempValildation, tempImputation);
        startV += numSamples; endV += numSamples;
        startI += numSamples; endI += numSamples;

    }
    return false;
}


void RSquare::printRSquare()
{
    for(int i = 0; i < numMarkers; i++)
    {
        cout << RSquareResult[i] << endl;
    }
}

bool RSquare::outputRSquare()
{
    fstream fs;
    fs.open(OutputPrefix+"RSquareOutput", ios_base::out);


    fs << "MarkerID\tRSquare\n";

    for (int i = 0; i < numMarkers; i++)
    {
        fs << MarkerID[i] << "\t" << RSquareResult[i] << "\n";
    }

    fs.close();

    cout << "Success! Please check RSquare result:" << OutputPrefix+"RSquareOutput" << endl;

    return false;
}


bool RSquare::Analysis()
{
    GetHaplotypeFromVCF();
    CalculateRSquare();
    outputRSquare();
    return false;
}