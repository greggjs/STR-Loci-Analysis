#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include "Sample.h"

typedef std::vector<Sample> SampleList;

SampleList populateLocusMapAndSampleList(LocusMap& locusMap, std::ifstream& inFile) {
    std::string line;
    SampleList sampleList;
    while(getline(inFile, line, '\r')) {
        Sample curr(line, locusMap);
        int i = 0;
        for (LocusPeak p : curr.getLocusPeaks()) {
            locusMap[i].addLocusPeaks(p.second);
            i++;
        }
        sampleList.push_back(curr);
    }
    return sampleList;
}

void calculateLocusProbs(LocusMap& locusMap, int sampleSize) {
    for (int i = 0; i < locusMap.size(); i++) {
        locusMap[i].calculateLocusProbs(sampleSize);
    }
}

int main(int argc, char* argv[]) {
    float pSignificance = atof(argv[2]);
    std::ifstream inFile(argv[1], std::ios::binary);

    std::string line;
    getline(inFile, line, '\r');
    LocusMap locusMap(line);
    SampleList sampleList = populateLocusMapAndSampleList(locusMap, inFile);
    calculateLocusProbs(locusMap, sampleList.size());
    //std::cout << locusMap << std::endl;


    return 0;
}
