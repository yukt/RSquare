//
// Created by Ketian Yu on 10/4/17.
//
#include "SummaryData.h"
#include "helperFunctions.h"

bool SummaryData::read()
{
    const int       MAXNUM = 100000, MAXNUM_SAME_POS = 100;
    VcfFileReader   inFileV, inFileI;
    VcfHeader       headerV, headerI;
    VcfRecord       recordV, recordI;
    VcfRecord       curV, curI;
    int             chrV, chrI, posV, posI;

    vector<String>        TempAltV;
    vector<SummaryRecord> TempRecordV;

    inFileV.open(FileNameValidation, headerV);
    inFileI.open(FileNameImputation, headerI);

    // check samples.
    if (headerI.getNumSamples() != headerV.getNumSamples())
    {
        error("[Error: ] Number of samples do NOT match!");
    }

    cout << "Reading VCFfiles ..." << endl;
    numSamples = headerV.getNumSamples();

    SumDat.resize(MAXNUM); markerID.resize(MAXNUM); REF.resize(MAXNUM); ALT.resize(MAXNUM);
    numRecords = 0;

    bool EndRecord = !(inFileV.readRecord(recordV) & inFileI.readRecord(recordI));
    chrV = chr2int(recordV.getChromStr()); chrI = chr2int(recordI.getChromStr());


    while (!EndRecord)
    {
        // read until CHROM of V and I are the same.
        bool ChrCheck = (chrV == chrI);
        while (!EndRecord and !ChrCheck)
        {
            while (chrV < chrI) { if(!(inFileV.readRecord(recordV))) {EndRecord = true; break;} chrV = chr2int(recordV.getChromStr());}
            while (chrV > chrI) { if(!(inFileI.readRecord(recordI))) {EndRecord = true; break;} chrI = chr2int(recordI.getChromStr());}
            if (chrV == chrI) { ChrCheck = true; }
        }
        if ( !ChrCheck ){ break; }

        // Under the same CHROM, read until POS of V and I are the same.
        int currentChr = chrV;
        posV = recordV.get1BasedPosition(); posI = recordI.get1BasedPosition();
        bool PosCheck = (posV == posI);
        while (!PosCheck & ChrCheck & !EndRecord)
        {
            while(posV < posI)
            {
                if(!(inFileV.readRecord(recordV)))     { EndRecord = true; break; }
                chrV = chr2int(recordV.getChromStr()); if(chrV > currentChr) { ChrCheck = false; break; }
                posV = recordV.get1BasedPosition();
            }
            while(posV > posI)
            {
                if(!(inFileI.readRecord(recordI)))     { EndRecord = true; break; }
                chrI = chr2int(recordI.getChromStr()); if(chrI > currentChr) { ChrCheck = false; break; }
                posI = recordI.get1BasedPosition();
            }
            if (ChrCheck and posV == posI) {PosCheck = true;}
        }

        if (!ChrCheck) { continue; }
        if (!PosCheck) { break;    }

        // start to compare ALT
        String altV = recordV.getAltStr(), altI = recordI.getAltStr();

        if (altV == altI)
        {
//            cout << numRecords << " " << chrV << " " << posV << " " << altV << endl;

            markerID[numRecords] = recordI.getIDStr();
            REF[numRecords]      = recordI.getRefStr();
            ALT[numRecords]      = recordI.getAltStr();
            SumDat[numRecords].resize(5,0);
            // [0]sumX [1]sumY [2]sumXY [3]sumX2 [4]sumY2

            VcfRecordGenotype& GenotypeV = recordV.getGenotypeInfo();
            VcfRecordGenotype& GenotypeI = recordI.getGenotypeInfo();
            vector<double> &temp = SumDat[numRecords];

            for (int i = 0; i < numSamples; i++)
            {
                double X = stod(*GenotypeV.getString("DS",i));
                double Y = stod(*GenotypeI.getString("DS",i));
                temp[0] += X; temp[1] += Y; temp[2] += X*Y; temp[3] += X*X; temp[4] += Y*Y;
            }

            numRecords++;


            EndRecord = !(inFileV.readRecord(recordV) and inFileI.readRecord(recordI));
            chrV = chr2int(recordV.getChromStr()); chrI = chr2int(recordI.getChromStr());
            continue;
        }

        else
        {
            TempAltV.clear();       TempAltV.resize(MAXNUM_SAME_POS);       TempAltV[0]    = altV;
            TempRecordV.clear();    TempRecordV.resize(MAXNUM_SAME_POS);    TempRecordV[0] = SummaryRecord(&recordV);
            int numTempV   = 1;
            int currentPos = posV;

            // read fileV until the next position or end of the fileV.
            while (true)
            {
                if (!inFileV.readRecord(recordV)) { EndRecord = true; break; }
                chrV = chr2int(recordV.getChromStr()); posV = recordV.get1BasedPosition();
                if (posV != currentPos or chrV != currentChr) { break; }
                TempAltV[numTempV] = recordV.getAltStr();
                TempRecordV[numTempV] = SummaryRecord(&recordV);
                numTempV++;
            }

            // read fileI and compare ALT until the next POS or END.
            while (true)
            {
                for (int j = 0; j < numTempV; j++)
                {
                    if (altI == TempAltV[j])
                    {
//                        cout << numRecords << " " << TempRecordV[j].ChromStr << " " << TempRecordV[j].BasedPosition <<  " " << altI << endl;

                        markerID[numRecords] = recordI.getIDStr();
                        REF[numRecords]      = recordI.getRefStr();
                        ALT[numRecords]      = recordI.getAltStr();
                        SumDat[numRecords].resize(5,0);
                        // [0]sumX [1]sumY [2]sumXY [3]sumX2 [4]sumY2

                        vector<double> &DosageV = TempRecordV[j].Dosage;
                        VcfRecordGenotype& GenotypeI = recordI.getGenotypeInfo();
                        vector<double> &temp = SumDat[numRecords];

                        for (int i = 0; i < numSamples; i++)
                        {
                            double X = DosageV[i];
                            double Y = stod(*GenotypeI.getString("DS",i));
                            temp[0] += X; temp[1] += Y; temp[2] += X*Y; temp[3] += X*X; temp[4] += Y*Y;
                        }

                        numRecords++;


                        break;
                    }
                }

                if (!inFileI.readRecord(recordI)) { EndRecord = true; break;}
                chrI = chr2int(recordI.getChromStr()); posI = recordI.get1BasedPosition();
                if (posI != currentPos or chrI != currentChr) { break; }
                altI = recordI.getAltStr();
            }

        }

    }
//
//    cout << endl;
//    cout << "Printing SumDat:" << endl;
//    for(int i=0; i < numRecords; i++)
//    {
//        cout << i << " ";
//        for (int j=0; j< 5; j++)
//            cout << SumDat[i][j] << " ";
//        cout << endl;
//    }
//

    return false;
}

void SummaryData::printData()
{
    for (int i = 0; i < numRecords; i++)
    {
        cout << i << "\t" << markerID[i] << "\t" << REF[i] << "\t" << ALT[i] << "\t";
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
    output();

    return false;
}

double SummaryData::vectorwiseRSquare(vector<int> index)
{
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
    double EX, EY, varX, varY, cov;

    int n = numSamples*int(index.size());

    for (int i: index)
    {
        vector<double> &temp = SumDat[i];
        sumX  += temp[0];
        sumY  += temp[1];
        sumXY += temp[2];
        sumX2 += temp[3];
        sumY2 += temp[4];
    }

    EX   = sumX *1.0/n;
    EY   = sumY *1.0/n;
    varX = sumX2*1.0/n - EX*EX;
    varY = sumY2*1.0/n - EY*EY;
    cov  = sumXY*1.0/n - EX*EY;

    return 1.0*cov/varX*cov/varY;

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

bool SummaryData::output()
{
    fstream fs;
    fs.open(OutputPrefix+"RSquareOutput", ios_base::out);


    fs << "MarkerID\tREF\tALT\tRSquare\n";

    for (int i = 0; i < numRecords; i++)
    {
        fs << markerID[i] << "\t" << REF[i] << "\t" << ALT[i] << "\t" << RSquareData[i] << "\n";
    }

    fs.close();

    cout << "Success! Please check RSquare result:" << OutputPrefix+"RSquareOutput" << endl;

    return false;
}