//
// Created by Ketian Yu on 9/27/17.
//

#ifndef RSQUARE_HAPLOTYPE_H
#define RSQUARE_HAPLOTYPE_H

#include <string>
#include <vector>
#include "VcfFileReader.h"
#include "StringBasics.h"


using namespace std;


class Haplotype{
public:
    String inFileName;
    int numMarkers, numSamples;
    vector<string> individualName;
    vector<string> markerID;

    vector<double> HaplotypeData;

    bool DS;


    bool read(String &VCFFileName);
    int getNumMarkers();
    int getNumSamples();

};

#endif //RSQUARE_HAPLOTYPE_H
