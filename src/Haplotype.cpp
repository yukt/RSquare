//
// Created by Ketian Yu on 9/27/17.
//
#include "Haplotype.h"

bool Haplotype::read(String &VCFFileName){

    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;

    inFileName = VCFFileName;
    inFile.open(VCFFileName, header);
    numSamples = header.getNumSamples();


    individualName.resize((unsigned long)numSamples);
    for (int i = 0; i < numSamples; i++)
    {
        string tempName(header.getSampleName((unsigned int)i));
        individualName[i] = tempName;
    }


    cout << "\nCollecting Haplotype Information ..." << endl;

    int numReadRecords = 0;

    while (inFile.readRecord(record))
    {
        numReadRecords++;
    }
    numMarkers = numReadRecords;
    inFile.close();


    markerID.resize((unsigned long)numMarkers);
    HaplotypeData.resize((unsigned long)numMarkers * numSamples);


    inFile.open(VCFFileName, header);

    if(DS){
        int i = 0;

        while (inFile.readRecord(record))
        {

            VcfRecordGenotype& Genotype = record.getGenotypeInfo();

            for(int j = 0; j < numSamples; j++)
            {
                HaplotypeData[i*numSamples+j] = stod(*Genotype.getString("DS",j));
                markerID[i] = record.getIDStr();
            }

            i++;

        }

        cout << "Successfully Got Dosage Data for " << i << " Markers and " << numSamples << " Samples!" << endl;
    }
    else
    {
        cout << "Error! This project is for DS format only. Still working on other formats :)" << endl;
        inFile.close();
        return false;
    }


    inFile.close();


    return true;

}

int Haplotype::getNumSamples()
{
    return numSamples;
}

int Haplotype::getNumMarkers()
{
    return numMarkers;
}
