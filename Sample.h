#include <cmath>
#include "LocusMap.h"

typedef std::pair<int, int> Peaks;
typedef std::pair<std::string, Key> LocusPeak;
typedef std::unordered_map<std::string, Key> LocusPeaks;

class Sample {
    friend std::ostream& operator<<(std::ostream& os, const Sample& sample) {
        os << sample.name << ": ";
        std::vector<LocusPeak> peaks(sample.allelePeaks.begin(), sample.allelePeaks.end());
        std::sort(peaks.begin(), peaks.end(),
            [](const LocusPeak& a, const LocusPeak& b){ return a.first < b.first; });
        for (LocusPeak peak : peaks) {
            os << "[" << peak.first << ": " << peak.second.first
                << ", " << peak.second.second << "], ";
        }
        os << "\n";
        return os;
    }

    public:
        Sample(const std::string& fileLine, LocusMap& locusMap) {
            std::string curr1, curr2;
            std::istringstream lineStream(fileLine);
            getline(lineStream, name, ',');
            Set locusSet = locusMap.getLocusSet();
            Set::iterator iter = locusSet.begin();
            while (lineStream && iter != locusSet.end()) {
                if (!getline(lineStream, curr1, ',')) {
                    break;
                }
                if (!getline(lineStream, curr2, ',')) {
                    break;
                }
                auto minAndMax = std::minmax({atoi(curr1.c_str()), atoi(curr2.c_str())});
                Key currPeaks = { minAndMax.first, minAndMax.second };
                allelePeaks[(*iter).getName()] = currPeaks;
                iter++;
            }
        }

        std::vector<LocusPeak> getLocusPeaks() const {
            std::vector<LocusPeak> peaks(allelePeaks.begin(), allelePeaks.end());
            std::sort(peaks.begin(), peaks.end(),
                [](const LocusPeak& a, const LocusPeak& b){ return a.first < b.first; });
            return peaks;
        }

    private:
        std::string name;
        LocusPeaks allelePeaks;
};
