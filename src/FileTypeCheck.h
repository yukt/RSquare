//
// Created by Ketian Yu on 10/3/17.
//

#include "SummaryData.h"

#ifndef RSQUARE_FILETYPECHECK_H
#define RSQUARE_FILETYPECHECK_H

#endif //RSQUARE_FILETYPECHECK_H

using namespace std;

bool vcfCheck(int option, String filename)
{
    // This is a slightly-modified version of the function provided by Sayantan
    // string HaplotypeSet::DetectFileType(String filename)

    IFILE fileStream = ifopen(filename, "r");
    string line;
    bool errorFlag = false;

    if(fileStream)
    {
        fileStream->readLine(line);
        if(line.length()<1)
        {
            ifclose(fileStream);
            errorFlag = true;
        }
        string tempString;
        tempString=(line.substr(0,17));

        char temp[tempString.length() + 1];
        std::strcpy(temp,tempString.c_str());
        for (char *iter = temp; *iter != '\0'; ++iter)
        {
            *iter = std::tolower(*iter);
        }

        if(((string)temp).compare("##fileformat=vcfv")==0)
        {
            ifclose(fileStream);
            return true;
        }
        else
        {
            ifclose(fileStream);
            errorFlag = true;
        }


    }
    else
    {
        errorFlag = true;
    }

    if(errorFlag)
    {
        error("[ERROR:] Invalid VCFfile -%c: %s\n", option, filename.c_str());
    }

    ifclose(fileStream);
    return false;

}

bool outputCheck(String filename)
{
    fstream fs(filename+".RSquareOutput", ios_base::out);
    if(fs.is_open())
    {
        fs.close();
    }
    else
    {
        error("[ERROR:] Invalid Prefix for Output -o: %s\n", filename.c_str());
    }
    return false;
}

bool formatCheck(int option, String format)
{
    if(format!="DS" and format!="GT" and format!="GP")
    {
        if(option == 'f')      error("[ERROR:] Invalid argument --validationFormat %s\n", format.c_str());
        else if(option == 'g') error("[ERROR:] Invalid argument --imputationFormat %s\n", format.c_str());
    }
    return false;
}

bool AFCheck(String filename)
{
    fstream fs(filename, ios_base::in);
    if(fs.is_open())
    {
        string SNP,AF;
        fs >> SNP >> AF;
        if (AF!="AF")
            error("[ERROR:] Please check the format of AlleleFrequency file --AF: %s\n", filename.c_str());
        fs.close();
    }
    else
    {
        error("[ERROR:] Allele Frequency file does not exist --AF: %s\n", filename.c_str());
    }
    return false;

}
