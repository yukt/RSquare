#include "Dosage.h"

bool Dosage::read(String &VCFFileName){
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;

    inFileName = VCFFileName;
    inFile.open(VCFFileName, header);
    numSamples = header.getNumSamples();


    for (int i = 0; i < numSamples; i++)
    {
        string tempName(header.getSampleName(i));
        individualName.push_back(tempName);
    }


    int numReadRecords = 0;

    cout << "Collecting Dosage Information ..." << endl;
    while (inFile.readRecord(record))
    {


        VcfRecordGenotype& Genotype = record.getGenotypeInfo();

        for(int j = 0; j < numSamples; j++)
        {
            DosageData.push_back(stod(*Genotype.getString("DS",j)));
            markerID.push_back(record.getIDStr());
        }

        numReadRecords++;

    }

    numMarkers = numReadRecords;

    cout << "Successfully Got Dosage for " << numMarkers << " Markers and " << numSamples << " Samples!" << endl;

    return 0;

}

int Dosage::getNumSamples()
{
    return numSamples;
}

int Dosage::getNumMarkers()
{
    return numMarkers;
}

