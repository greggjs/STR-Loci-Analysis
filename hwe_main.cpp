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
    std::cout << "Populating Locus Map..." << std::endl;
    while (getline(inFile, line, '\r')) {
        Sample curr(line, locusMap);
        int i = 0;
        std::vector<LocusPeak> locusPeaks = curr.getLocusPeaks();
        for (int j = 0; j < locusPeaks.size(); j++) {
            locusMap[i].addLocusPeaks(locusPeaks[j].second);
            for (int k = j + 1; k < locusPeaks.size(); k++) {
                if (locusPeaks[j].first != locusPeaks[k].first) {
                    locusMap[i].addLocusPair(locusPeaks[j].second, locusPeaks[k]);
                }
            }
            i++;
        }

        sampleList.push_back(curr);
    }
    return sampleList;
}

void calculateLocusProbs(LocusMap& locusMap, int sampleSize, float psig) {
    for (int i = 0; i < locusMap.size(); i++) {
        locusMap[i].calculateLocusProbs(sampleSize);
        locusMap[i].calculateAlleleProbs(sampleSize);
        locusMap[i].calculateE(sampleSize);
    }
}

void calculateHWEandLinkage(LocusMap& locusMap, int sampleSize, float psig){
    for (int i = 0; i < locusMap.size() - 1; i++) {
        locusMap[i].isHWE(psig);
    }
    for (int i = 0; i < locusMap.size() - 1; i++) {
        for (int j = i + 1; j < locusMap.size(); j++) {
            locusMap[i].doLinkageCompares(locusMap[j], sampleSize, psig);
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Need 2 args (file and pSignificance)." << std::endl;
        return 1;
    }
    float pSignificance = atof(argv[2]);
    std::ifstream inFile(argv[1], std::ios::binary);

    std::string line;
    getline(inFile, line, '\r');
    LocusMap locusMap(line);
    SampleList sampleList = populateLocusMapAndSampleList(locusMap, inFile);
    calculateLocusProbs(locusMap, sampleList.size() - 1, pSignificance);
    calculateHWEandLinkage(locusMap, sampleList.size() - 1, pSignificance);
    //std::cout << sampleList.size() << std::endl;
    //std::cout << locusMap << std::endl;


    return 0;
}
