#ifndef RSQUARE_DOSAGE_H
#define RSQUARE_DOSAGE_H

#include <string>
#include <vector>
#include "VcfFileReader.h"
#include "StringBasics.h"


using namespace std;


class Dosage{
    public:
        String inFileName;
        int numMarkers, numSamples;
        vector<string> individualName;
        vector<string> markerID;

        vector<double> DosageData;


        bool read(String &VCFFileName);
        int getNumMarkers();
        int getNumSamples();

};

#endif //RSQUARE_DOSAGE_H
