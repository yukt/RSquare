//
// Created by Ketian Yu on 10/4/17.
//

#ifndef RSQUARE_SUMMARYDATA_H
#define RSQUARE_SUMMARYDATA_H

#include <string>
#include <vector>
#include <Error.h>
#include <fstream>
#include <sstream>
#include "VcfFileReader.h"



using namespace std;

class SummaryData
{
public:
    String FileNameValidation, FileNameImputation, OutputPrefix;
    int numRecords, numSamples;
    vector<string> SNP;
    vector<vector<double>> RSquareData;
    // [0]numObsGeno   [1]GoldFreq  [2]RSquare
    vector<vector<double>> SumDat;
    // [0]sumX [1]sumY [2]sumXY [3]sumX2 [4]sumY2 [5]n
    string validationFormat, imputationFormat;

    SummaryData()
    {
        validationFormat = "GT";
        imputationFormat = "DS";
    };

    ~SummaryData()
    {

    };

    bool   read();
    void   printData();
    bool   analysis();
    bool   RSquare();
    void   printRSquare();
    bool   output();

private:
    vector<double> vectorwiseRSquare( vector<int> index );


};

#endif //RSQUARE_SUMMARYDATA_H



