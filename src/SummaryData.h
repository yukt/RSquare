//
// Created by Ketian Yu on 10/4/17.
//

#ifndef RSQUARE_SUMMARYDATA_H
#define RSQUARE_SUMMARYDATA_H

#include <string>
#include <vector>
#include <Error.h>
#include <fstream>
#include "VcfFileReader.h"



using namespace std;

class SummaryData
{
public:
    String FileNameValidation, FileNameImputation, OutputPrefix;
    int numRecords, numSamples;
    vector<string> markerID, REF, ALT;
    vector<double> RSquareData;
    vector<vector<double>> SumDat;
    // [0]sumX [1]sumY [2]sumXY [3]sumX2 [4]sumY2

    SummaryData()
    {

    };

    ~SummaryData()
    {

    };

    bool   read();
    void   printData();
    bool   analysis();
    double vectorwiseRSquare( vector<int> index );
    bool   RSquare();
    bool   output();


};

#endif //RSQUARE_SUMMARYDATA_H



