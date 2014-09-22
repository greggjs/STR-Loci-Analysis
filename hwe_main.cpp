#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include "Sample.h"

typedef std::vector<Sample> SampleList;

int main(int argc, char* argv[]) {
    float pSignificance = atof(argv[2]);
    std::ifstream inFile(argv[1], std::ios::binary);

    std::string line;
    getline(inFile, line, '\r');
    LocusMap locusMap(line);

    SampleList sampleList;
    while(getline(inFile, line, '\r')) {
        Sample curr(line, locusMap);
        int i = 0;
        for (LocusPeak p : curr.getLocusPeaks()) {
            locusMap[i].addLocusPeak(p.second.first);
            locusMap[i].addLocusPeak(p.second.second);
            i++;
        }
        sampleList.push_back(curr);
    }
    std::cout << sampleList[0] << std::endl;
    std::cout << locusMap << std::endl;

    return 0;
}
