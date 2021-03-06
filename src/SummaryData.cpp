//
// Created by Ketian Yu on 10/4/17.
//
#include "SummaryData.h"
#include "helperFunctions.h"
#include <iomanip>
#include <ctime>
#include <cmath>
#include <VcfFileReader.h>

bool SummaryData::read()
{
    VcfFileReader   inFileV, inFileI;
    VcfHeader       headerV, headerI;
    VcfRecord       recordV, recordI;
    int             posV,chrV,posI,chrI;

    inFileV.open(FileNameValidation, headerV);
    inFileI.open(FileNameImputation, headerI);

    double (*getDosageV) (VcfRecordGenotype&, int ) = &readGT;
    double (*getDosageI) (VcfRecordGenotype&, int ) = &readDS;
    if      (validationFormat == "DS") { getDosageV = &readDS; }
    else if (validationFormat == "GP") { getDosageV = &readGP; }
    if      (imputationFormat == "GT") { getDosageI = &readGT; }
    else if (imputationFormat == "GP") { getDosageI = &readGP; }

    SumDat.resize(NumMax); SNP.resize(NumMax); commonIndex.resize(NumMax);
    numRecords = 0;

    bool EndRecord = !inFileV.readRecord(recordV); chrV = chr2int(recordV.getChromStr());
    int indexInImputed = -1, indexInValidation = 0;
    while (inFileI.readRecord(recordI) & !EndRecord)
    {
        indexInImputed++;
        chrI = chr2int(recordI.getChromStr());
        while ((chrV<chrI) & !EndRecord)
        {
            EndRecord = !inFileV.readRecord(recordV);
            chrV = chr2int(recordV.getChromStr());
            indexInValidation++;
        }
        if    (chrV > chrI) continue;

        posI = recordI.get1BasedPosition(); posV = recordV.get1BasedPosition();
        while ((posV<posI) & !EndRecord)
        {
            EndRecord = !inFileV.readRecord(recordV);
            chrV = chr2int(recordV.getChromStr());
            posV = recordV.get1BasedPosition();
            indexInValidation++;
        }
        if    (chrV > chrI or posV > posI) continue;

        stringstream ssV, ssI; ssV << posV; ssI << posI;
        string markerV = string(recordV.getChromStr()) + ":" + ssV.str()+":"+recordV.getRefStr()+":"+recordV.getAltStr();
        string markerI = string(recordI.getChromStr()) + ":" + ssI.str()+":"+recordI.getRefStr()+":"+recordI.getAltStr();
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
        commonIndex[numRecords] = indexInValidation;
        numRecords++;
    }

    SumDat.resize(numRecords); SNP.resize(numRecords); commonIndex.resize(numRecords);
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


vector<double> SummaryData::vectorwiseRSquare(vector<int> index)
{
    vector<double> result;
    result.resize(5);

    // X - validation Y - Imputed
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
    double EX, EY, varX, varY, cov;
    int numGeno = 0;

    for (int i: index)
    {
        vector<double> &temp = SumDat[i];
        sumX    += temp[0];
        sumY    += temp[1];
        sumXY   += temp[2];
        sumX2   += temp[3];
        sumY2   += temp[4];
        numGeno += temp[5];
    }

    result[0] = numGeno;
    if (numGeno == 0){
        result[1] = result[2] = result[3] = result[4] = 0;
        return result;
    }

    EX   = sumX *1.0/numGeno;
    EY   = sumY *1.0/numGeno;
    varX = sumX2*1.0/numGeno - EX*EX;
    varY = sumY2*1.0/numGeno - EY*EY;
    cov  = sumXY*1.0/numGeno - EX*EY;

    result[1] = min(EX,2-EX)*0.5;
    result[3] = EX*0.5;
    result[4] = EY*0.5;
    if (varX == 0 or varY == 0){
        result[2] = 0;
        return result;
    }
    result[2] = 1.0*cov/varX*cov/varY;
    return result;

}

bool SummaryData::RSquare()
{
    RSquareData.resize(numRecords);
    for (int i = 0; i < numRecords; i++)
    {
        RSquareData[i] = vectorwiseRSquare({i});
    }
    return false;
}

void SummaryData::printRSquare()
{
    cout << "SNP\tnumObsGeno\tGoldMAF\tGoldAltFreq\tImputedAltFreq\tRSquare\n";
    for (int i = 0; i < numRecords; i++){
        vector<double> &temp = RSquareData[i];
        cout << SNP[i] << "\t";
        cout << setprecision(6) << temp[0] << "\t" << temp[1] << "\t" << temp[3] << "\t" << temp[4] << "\t" << temp[2] << "\n";
    }
}

bool SummaryData::output()
{
    fstream fs;
    fs.open(OutputPrefix+".RSquare", ios_base::out);
    fs << std::fixed << std::setprecision(6);
    fs << "SNP\tnumObsGeno\tGoldMAF\tGoldAltFreq["<< validationFormat << "]\tImputedAltFreq["<< imputationFormat << "]\tRSquare\n";

    for (int i = 0; i < numRecords; i++)
    {
        vector<double> &temp = RSquareData[i];
        fs << SNP[i] << "\t";
        fs << (int)temp[0] << "\t" << temp[1] << "\t" << temp[3] << "\t" << temp[4] << "\t" << temp[2] << "\n";
    }

    fs.close();
    return false;
}

bool SummaryData::loadNumMax()
{
    VcfFileReader inFile;
    VcfHeader     header;
    VcfRecord     record;

    inFile.open(FileNameImputation, header);
    NumMax = 0;
    while (inFile.readRecord(record)) { NumMax++; }
    inFile.close();
    return false;
}

bool SummaryData::sampleCheck()
{
    VcfFileReader   inFileV, inFileI;
    VcfHeader       headerV, headerI;
    inFileV.open(FileNameValidation, headerV);
    inFileI.open(FileNameImputation, headerI);

    // check samples.
    if (headerI.getNumSamples() != headerV.getNumSamples())
    {
        error("[Error] Number of samples do NOT match!");
    }
    numSamples = headerV.getNumSamples();
    inFileI.close(); inFileV.close();
    return false;
}

void SummaryData::initiateBins()
{
    static const double arr[] = {0,0.0005,0.001,0.002,0.005,0.010,0.015,0.020,0.035,0.050};
    bins.assign(arr, arr + sizeof(arr) / sizeof(arr[0]));
    bins.resize(20);
    for(int i=1; i<11; i++) { bins[9+i]=i*0.1; }
}

bool SummaryData::loadBins()
{
    ifstream inFile;
    inFile.open(FileNameBins);

    bins.clear();
    bins.resize(99);

    double x, y=-1e-4;
    int count=0;
    while(inFile >> x)
    {
        if(x < y || x > 1){ error("[Error] Please check the bin file.");}
        bins[count] = x;
        y = x;
        count++;
    }
    bins.resize(count);
    inFile.close();
    return 0;
}

bool SummaryData::loadAlleleFreq()
{
    fstream fs;
    fs.open(FileNameAlleleFreq, ios_base::in);
    string header, line;

    double AF;
    vector<int> counter;

    int colCount=0, popCount=0;
    int col_AF=-1;

    char* pch, * end_str;
    string tempString;

    getline(fs, header);
    pch = strtok_r ((char*)header.c_str(),"\t", &end_str);

    while(pch!=NULL)
    {
        colCount++;
        if(colCount>1)
        {
            tempString=string(pch);
            if(tempString=="AF") { col_AF=colCount; }
            else { pops.push_back(tempString); popCount++; }
        }
        pch = strtok_r (NULL,"\t", &end_str);
    }



    int nBins = bins.size();
    aggregateIndex.resize(nBins);
//    aggregateIndexInGold.resize(nBins);
    aggregateAF.resize(nBins);
    counter.resize(nBins);

    for(int i=0; i<nBins; i++)
    {
        aggregateIndex[i].resize(commonIndex.size());
//        aggregateIndexInGold[i].resize(commonIndex.size());
    }

    int nRecordReadAF = -1, nRecordReadCommonIndex = -1;
    vector<string>::iterator iter = SNP.begin();
    numCommonSNPsAnalyzed = 0;

    for(int index : commonIndex)
    {
        nRecordReadCommonIndex++;
        while (nRecordReadAF < index) { getline(fs,line); nRecordReadAF++; }
        pch = strtok_r ((char*)line.c_str(),"\t", &end_str);
        string currentSNP = *iter;
        if(string(pch)!= currentSNP) { error("[Error] Allele Frequency File Does NOT Match Validation VCF File!\n"); }
        colCount = 1;

        while(pch!=NULL)
        {
            colCount++;
            pch = strtok_r (NULL,"\t", &end_str);
            if(colCount == col_AF)
            {
                if(pch==NULL)
                    break;
                AF = stod(pch);
                int group = nBins-1;

                if(AF > bins[0]) { for(int i=1; i<nBins; i++) { if( AF < bins[i] ) { group = i-1; break; }}}

                aggregateIndex[group][counter[group]] = nRecordReadCommonIndex;
//                aggregateIndexInGold[group][counter[group]] = nRecordReadAF;
                aggregateAF[group] += AF;
                counter[group]++;
                numCommonSNPsAnalyzed++;
                break;
            }
        }

        iter++;
    }

    for(int i=0; i<nBins; i++)
    {
        aggregateIndex[i].resize(counter[i]);
//        aggregateIndexInGold[i].resize(counter[i]);
        aggregateAF[i]/=counter[i];
    }
    return 0;
}

bool SummaryData::aggregate()
{
    unsigned long numSNP;
    vector<double> result;
    aggregateRSquare.resize(aggregateIndex.size()-1);
    for (int i = 0; i <aggregateRSquare.size(); i++)
    {
        numSNP = aggregateIndex[i].size();
        aggregateRSquare[i].resize(4);
        vector<double> &temp = aggregateRSquare[i];
        // aggregateRSquare [0]numSNP [1]GoldenFreq [2]ImputedFreq [3]RSquare

        temp[0] = numSNP;
        if (numSNP==0) { temp[1] = temp[2] = temp[3] = 0; }
        else
        {
            result = vectorwiseRSquare(aggregateIndex[i]);
            //result [0]numGeno [1]MAF [2]RSquare [3]GoldenAltFreq [4]ImputedAltFreq
            temp[1]=result[3]; temp[2]=result[4]; temp[3] = result[2];
        }
    }

    return false;
}

bool SummaryData::outputAggregate()
{
    fstream fs;
    fs.open(OutputPrefix+".aggRSquare", ios_base::out);
    fs << std::fixed << std::setprecision(6);
    fs << "AlleleFreq\tRSquare\t\tnumSNP\tGoldenFreq\tImputedFreq\n";

    for (int i = 0; i < aggregateRSquare.size(); i++)
    {
        vector<double> &temp = aggregateRSquare[i];
        if(temp[0]>0) { fs << aggregateAF[i] << "\t" << temp[3] << "\t" << (int)temp[0] << "\t" << temp[1] << "\t" << temp[2] << "\n"; }
    }

    fs.close();

    return false;
}

//bool SummaryData::outputIndex()
//{
//    fstream fs1, fs2;
//    fs1.open(OutputPrefix+".bin.0.common.index", ios_base::out);
//    fs2.open(OutputPrefix+".bin.0.golden.index", ios_base::out);
//    vector<int> &temp1 = aggregateIndex[0];
//    vector<int> &temp2 = aggregateIndexInGold[0];
//    for(int j = 0; j < temp1.size(); j++) { fs1 << temp1[j] << "\n"; fs2 << temp2[j] << "\n"; }
//    fs1.close();
//    fs2.close();
//
////    for(int i = 0; i < aggregateIndex.size(); i++)
////    {
////        fstream fs;
////        fs.open(OutputPrefix+".bin."+to_string(i).c_str()+".index", ios_base::out);
////        vector<int> &temp = aggregateIndex[i];
////        for(int j = 0; j < temp.size(); j++) { fs << temp[j] << "\n"; }
////        fs.close();
////    }
//
//    return 0;
//}


bool SummaryData::analysis()
{
    time_t  startTime, endTime;
    struct tm * timeinfo;
    double secondPassed;

    sampleCheck();

    if (FileNameAlleleFreq != ""){
        if (FileNameBins !="") { loadBins(); }
        else { initiateBins(); }
    }

    time(&startTime); timeinfo = localtime(&startTime);
    cout << "[INFO] Analysis started at: " << asctime(timeinfo);
    cout << "[INFO] Loading "<< validationFormat << " information from validation vcffile: " + FileNameValidation << " ..." << endl;
    read(); cout << "[INFO] Loaded [ " << numRecords << " ] common records\n" << "[INFO] Calculating RSquare ..." << endl;
    RSquare();
    output(); cout << "[INFO] Success! RSquare result: " + OutputPrefix + ".RSquare" << endl;
    if (FileNameAlleleFreq != "")
    {
        cout << "[INFO] Calculating aggregate RSquare ..." << endl;
        loadAlleleFreq();
        aggregate();
        outputAggregate();
//        outputIndex();
        cout << "[INFO] Success! Aggregate RSquare result: " + OutputPrefix+".aggRSquare" << endl;
    }

    time(&endTime); timeinfo = localtime(&endTime);
    secondPassed = difftime(endTime,startTime);
    cout << "[INFO] Analyzed [ " << numCommonSNPsAnalyzed << " ] SNPs" << endl;
    cout << "[INFO] Analysis ends at: " <<  asctime(timeinfo);
    cout << "[INFO] Analysis took " << secondPassed << " seconds" << endl;

    return false;
}
