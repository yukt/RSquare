//
// Created by Ketian Yu on 10/4/17.
//
#include "SummaryData.h"
#include "helperFunctions.h"
#include <iomanip>

bool SummaryData::read()
{
//    const int       MAXNUM = 10000000;
    VcfFileReader   inFileV, inFileI;
    VcfHeader       headerV, headerI;
    VcfRecord       recordV, recordI;

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

    SumDat.resize(maxNum); CHROM.resize(maxNum); POS.resize(maxNum); REF.resize(maxNum); ALT.resize(maxNum);
    numRecords = 0;

    while (inFileI.readRecord(recordI))
    {
        while (inFileV.readRecord(recordV) && recordV.get1BasedPosition() < recordI.get1BasedPosition());
        if    (recordI.get1BasedPosition() <  recordV.get1BasedPosition())
            continue;

        stringstream ssV, ssI;
        ssV << recordV.get1BasedPosition(); ssI << recordI.get1BasedPosition();

        string markerV = ssV.str() + "_" + recordV.getRefStr() + "_" + recordV.getAltStr();
        string markerI = ssI.str() + "_" + recordI.getRefStr() + "_" + recordI.getAltStr();
        if (markerV != markerI)
            continue;

        // Now the markers of records match.
        CHROM[numRecords] = recordI.getChromStr();
        POS[numRecords]   = recordI.get1BasedPosition();
        REF[numRecords]   = recordI.getRefStr();
        ALT[numRecords]   = recordI.getAltStr();
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

    SumDat.resize(numRecords); CHROM.resize(numRecords); POS.resize(numRecords); REF.resize(numRecords); ALT.resize(numRecords);


    return false;
}

void SummaryData::printData()
{
    for (int i = 0; i < numRecords; i++)
    {
        cout << i << "\t" << CHROM[i] << "\t" << POS[i] << "\t" << REF[i] << "\t" << ALT[i] << "\t";
        vector<double> &temp = SumDat[i];
        for (int j = 0; j < 5; j++)
            cout << temp[j] << "\t";
        cout << endl;
    }
}

bool SummaryData::analysis()
{
    read();
//    printData();
    RSquare();
//    printRSquare();
    output();

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

    EX   = sumX *1.0/n;
    EY   = sumY *1.0/n;
    varX = sumX2*1.0/n - EX*EX;
    varY = sumY2*1.0/n - EY*EY;
    cov  = sumXY*1.0/n - EX*EY;

    result[0] = n;
    result[1] = EX*0.5;
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
        cout << CHROM[i] << ":" << POS[i] << ":" << REF[i] << ":" << ALT[i] << "\t";
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
        fs << CHROM[i] << ":" << POS[i] << ":" << REF[i] << ":" << ALT[i] << "\t";
        fs << (int)temp[0] << "\t" << temp[1] << "\t" << temp[2] << "\n";
    }

    fs.close();

    cout << "Success! Please check RSquare result:" << OutputPrefix+".RSquareOutput" << endl;

    return false;
}