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
    String FileNameValidation, FileNameImputation, OutputPrefix, FileAF;
    string validationFormat, imputationFormat;
    int numRecords, numSamples;
    int NumMax;

    // RSquare results:
    vector<string> SNP;
    vector<vector<double>> RSquareData;
    // [0]numObsGeno   [1]GoldFreq  [2]RSquare [3]ImputedFreq
    vector<vector<double>> SumDat;
    // [0]sumX [1]sumY [2]sumXY [3]sumX2 [4]sumY2 [5]n


    // for aggregate use:
    bool   makeAF_flag;
    vector <int> commonIndex;
    vector <vector<int>> aggregateIndex;
    vector <vector<double>> aggregateRSquare;

    SummaryData()
    {
        validationFormat = "GT";
        imputationFormat = "DS";
        makeAF_flag = false;
    };

    ~SummaryData()
    {

    };

    bool   analysis();


private:
    vector<double> vectorwiseRSquare( vector<int> index );
    void   printData();
    void   printRSquare();
    bool   sampleCheck();
    bool   loadNumMax();
    bool   loadAlleleFreq();
    bool   read();
    bool   RSquare();
    bool   output();
    bool   aggregate();
    bool   outputAggregate();
    bool   makeAF();


};

#endif //RSQUARE_SUMMARYDATA_H



