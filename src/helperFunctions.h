//
// Created by Ketian Yu on 10/4/17.
//

#ifndef RSQUARE_HELPERFUNCTIONS_H
#define RSQUARE_HELPERFUNCTIONS_H

#include <StringBasics.h>
#include <VcfRecordGenotype.h>
#include <VcfRecord.h>

double readGT(VcfRecordGenotype& Genotype, int i)
{
    string info = *Genotype.getString("GT", i);
    if (info == ".") return -1;
    double dosage = 0;
    for (char j : info) {
        if( j == '1' ) dosage++;
        else if (j == '.') return -1;
    }
    return dosage;
}

double readDS(VcfRecordGenotype& Genotype, int i)
{
    string info = *Genotype.getString("DS", i);
    if (info == ".") return -1;
    return stod(info);
}

double readGP(VcfRecordGenotype& Genotype, int i)
{
    string info = *Genotype.getString("GP", i);
    if (info == ".") return -1;
    string delimiter = ",";
    size_t pos = info.find(delimiter); info.erase(0, pos+1); pos = info.find(delimiter);
    string GP_01 = info.substr(0, pos); info.erase(0, pos+1);
    string GP_11 = info;
    if (GP_01 == "." or GP_11 == ".") return -1;
    return stod(GP_01) + 2*stod(GP_11);
}


int chr2int(String chr)
{
    if(chr == "X") { return 23; }
    if(chr == "Y") { return 24; }
    return chr;
}

class SummaryRecord{
public:
    String ChromStr;
    int BasedPosition;
    String RefStr, AltStr;
    int NumSamples;
    vector<double> Dosage;

    SummaryRecord(){
        ChromStr = "";
        BasedPosition = 0;
        RefStr = "";
        AltStr = "";
        NumSamples = 0;
    }
    ~SummaryRecord(){

    }

    SummaryRecord(VcfRecord *record)
    {
        ChromStr = record->getChromStr();
        BasedPosition = record->get1BasedPosition();
        RefStr = record->getRefStr();
        AltStr = record->getAltStr();
        NumSamples = record->getNumSamples();
        VcfRecordGenotype& Genotype = record->getGenotypeInfo();

        Dosage.resize((unsigned long)NumSamples);
        for (int i = 0; i < NumSamples; i++)
        {
            Dosage[i] = stod(*Genotype.getString("DS",i));
        }
    }

};

#endif //RSQUARE_HELPERFUNCTIONS_H
