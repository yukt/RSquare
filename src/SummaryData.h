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
#include <StringBasics.h>


using namespace std;

class SummaryData
{
public:
    String FileNameValidation, FileNameImputation, OutputPrefix, FileNameAlleleFreq, FileNameBins;
    string validationFormat, imputationFormat;
    int numRecords, numSamples, numCommonSNPsAnalyzed;
    int NumMax;

    // RSquare results:
    vector<string> SNP;
    vector<vector<double>> RSquareData;
    // [0]numObsGeno   [1]MAF  [2]RSquare [3]GoldAltFreq [4]ImputedAltFreq
    vector<vector<double>> SumDat;
    // [0]sumX [1]sumY [2]sumXY [3]sumX2 [4]sumY2 [5]n


    // for aggregate use:
    vector <int> commonIndex;
    vector <double> bins;
    vector <string> pops;
    vector <vector<int>> aggregateIndex;
//    vector <vector<int>> aggregateIndexInGold;
    vector <double> aggregateAF;
    vector <vector<double>> aggregateRSquare;

    SummaryData()
    {
        validationFormat = "GT";
        imputationFormat = "DS";
        NumMax=99999999;
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
    bool   loadBins();
    void   initiateBins();
    bool   loadAlleleFreq();
    bool   read();
    bool   RSquare();
    bool   output();
    bool   aggregate();
    bool   outputAggregate();
    bool   outputIndex();


};

#endif //RSQUARE_SUMMARYDATA_H



