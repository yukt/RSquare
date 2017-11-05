//
// Created by Ketian Yu on 10/4/17.
//
#include "SummaryData.h"
#include "helperFunctions.h"
#include <iomanip>
#include <ctime>

bool SummaryData::read()
{
    VcfFileReader   inFileV, inFileI;
    VcfHeader       headerV, headerI;
    VcfRecord       recordV, recordI;
    int             posV,chrV,posI,chrI;

    inFileV.open(FileNameValidation, headerV);
    inFileI.open(FileNameImputation, headerI);

    // check samples.
    if (headerI.getNumSamples() != headerV.getNumSamples())
    {
        error("[Error: ] Number of samples do NOT match!");
    }

    double (*getDosageV) (VcfRecordGenotype&, int ) = &readGT;
    double (*getDosageI) (VcfRecordGenotype&, int ) = &readDS;
    if (validationFormat == "DS") { getDosageV = &readDS; }
    if (imputationFormat == "GT") { getDosageI = &readGT; }

    cout << "Reading VCFfiles ..." << endl;
    numSamples = headerV.getNumSamples();

    unsigned long maxNum = 0;
    while (inFileI.readRecord(recordI)) { maxNum++; }
    inFileI.close();
    inFileI.open(FileNameImputation, headerI);

    SumDat.resize(maxNum); SNP.resize(maxNum);
    numRecords = 0;

    bool EndRecord = !inFileV.readRecord(recordV); chrV = chr2int(recordV.getChromStr());
    while (inFileI.readRecord(recordI) & !EndRecord)
    {
        chrI = chr2int(recordI.getChromStr());
        while ((chrV<chrI) & !EndRecord) {EndRecord = !inFileV.readRecord(recordV); chrV = chr2int(recordV.getChromStr());}
        if    (chrV > chrI) continue;

        posI = recordI.get1BasedPosition(); posV = recordV.get1BasedPosition();
        while ((posV<posI) & !EndRecord) {EndRecord = !inFileV.readRecord(recordV); chrV = chr2int(recordV.getChromStr()); posV = recordV.get1BasedPosition();}
        if    (chrV > chrI or posV > posI) continue;

        stringstream ssV, ssI; ssV << posV; ssI << posI;
        string markerV = string(recordV.getChromStr()) + ":" + ssV.str()+":"+recordV.getRefStr()+":"+recordV.getAltStr();
        string markerI = string(recordI.getChromStr()) + ":" + ssI.str()+":"+recordV.getRefStr()+":"+recordV.getAltStr();
        if (markerV != markerI)
            continue;

        // Now the markers of records match.
        SNP[numRecords] = markerI;
        SumDat[numRecords].resize(6,0);
        // [0]sumX [1]sumY [2]sumXY [3]sumX2 [4]sumY2 [5]n

        VcfRecordGenotype& GenotypeV = recordV.getGenotypeInfo();
        VcfRecordGenotype& GenotypeI = recordI.getGenotypeInfo();
        vector<double> &temp = SumDat[numRecords];

        for (int i = 0; i < numSamples; i++)
        {
            double X = (*getDosageV)(GenotypeV, i);
            double Y = (*getDosageI)(GenotypeI, i);
            if (X>=0 and Y>=0){
                temp[0] += X; temp[1] += Y; temp[2] += X*Y; temp[3] += X*X; temp[4] += Y*Y; temp[5]++;
            }
        }

        numRecords++;
//        cout << numRecords << endl;

    }
    cout << "Finish reading " << numRecords << " Common Records." << endl;

    SumDat.resize(numRecords); SNP.resize(numRecords);
    inFileI.close(); inFileV.close();
    return false;
}

void SummaryData::printData()
{
    for (int i = 0; i < numRecords; i++)
    {
        cout << i << "\t" << SNP[i] << "\t";
        vector<double> &temp = SumDat[i];
        for (int j = 0; j < 5; j++)
            cout << temp[j] << "\t";
        cout << endl;
    }
}

bool SummaryData::analysis()
{
    time_t  startTime, endTime;
    struct tm * timeinfo;

    time(&startTime);

    read();
//    printData();
    RSquare();
//    printRSquare();
    output();

    time(&endTime);
    timeinfo = localtime (&endTime);
    double secondPassed = difftime(endTime,startTime);
    cout << "[INFO]	Analyzed [ " << numRecords << " ] SNPs" << endl;
    cout << "[INFO]	Analysis ends at: " <<  asctime(timeinfo);
    cout << "[INFO]	Analysis took " << secondPassed << " seconds" << endl;

    return false;
}

vector<double> SummaryData::vectorwiseRSquare(vector<int> index)
{
    vector<double> result;
    result.resize(3);

    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
    double EX, EY, varX, varY, cov;
    int n = 0;

    for (int i: index)
    {
        vector<double> &temp = SumDat[i];
        sumX  += temp[0];
        sumY  += temp[1];
        sumXY += temp[2];
        sumX2 += temp[3];
        sumY2 += temp[4];
        n     += temp[5];
    }

    result[0] = n;
    if (n == 0){
        result[1] = result[2] = 0;
        return result;
    }

    EX   = sumX *1.0/n;
    EY   = sumY *1.0/n;
    varX = sumX2*1.0/n - EX*EX;
    varY = sumY2*1.0/n - EY*EY;
    cov  = sumXY*1.0/n - EX*EY;

    result[1] = min(EX,EY)*0.5;
    if (varX == 0 or varY == 0){
        result[2] = 0;
        return result;
    }
    result[2] = 1.0*cov/varX*cov/varY;
    return result;

}

bool SummaryData::RSquare()
{
    cout << "Calculating RSquare ..." << endl;
    RSquareData.resize(numRecords);
    for (int i = 0; i < numRecords; i++)
    {
        RSquareData[i] = vectorwiseRSquare({i});
    }
    return false;
}

void SummaryData::printRSquare()
{
    cout << "SNP\tnumObsGeno\tGoldFreq\tRSquare\n";
    for (int i = 0; i < numRecords; i++){
        vector<double> &temp = RSquareData[i];
        cout << SNP[i] << "\t";
        cout << setprecision(6) << temp[0] << "\t" << temp[1] << "\t" << temp[2] << "\n";
    }
}

bool SummaryData::output()
{
    fstream fs;
    fs.open(OutputPrefix+".RSquareOutput", ios_base::out);
    fs << std::fixed << std::setprecision(6);
    fs << "SNP\tnumObsGeno\tGoldFreq\tRSquare\n";

    for (int i = 0; i < numRecords; i++)
    {
        vector<double> &temp = RSquareData[i];
        fs << SNP[i] << "\t";
        fs << (int)temp[0] << "\t" << temp[1] << "\t" << temp[2] << "\n";
    }

    fs.close();

    cout << "Success! Please check RSquare result: " << OutputPrefix+".RSquareOutput" << endl;

    return false;
}