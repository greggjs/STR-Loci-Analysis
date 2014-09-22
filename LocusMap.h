#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>

typedef std::unordered_map<int, int> LocusDist;

class Locus {
    friend std::ostream& operator<<(std::ostream& os, const Locus& a) {
        os << a.name << ": ";
        for (std::pair<int, int> p : a.locusDist) {
            os << "[" << p.first << ": " << p.second <<"],";
        }
        return os;
    }
    public:
        Locus(const std::string name) : name(name) {}
        bool operator<(const Locus& locus) const { return name < locus.name; }
        bool operator==(const Locus& locus) const { return name == locus.name; }

        std::string getName() const { return name; }
        void addLocusPeak(int peak) { locusDist[peak]++; }
    private:
        std::string name;
        LocusDist locusDist;
};

typedef std::vector<Locus> Set;

class LocusMap {
    friend std::ostream& operator<<(std::ostream& os, const LocusMap& lMap) {
        os << "LocusMap (size: " << lMap.locusSet.size() << "):\n";
        std::copy(lMap.locusSet.begin(), lMap.locusSet.end(),
            std::ostream_iterator<Locus>(os, "\n"));
        return os;
    }

    public:
        LocusMap(const std::string locusString) {
            std::string locus;
            std::istringstream locusStream(locusString);
            getline(locusStream, locus, ',');
            while(locusStream) {
                if (!getline(locusStream, locus, ',')) {
                    break;
                }
                locusSet.push_back(locus);
                getline(locusStream, locus, ',');
            }
        }

        Locus operator[](int index) const { return locusSet.at(index); }
        Locus& operator[](int index) { return locusSet.at(index); }

        Set& getLocusSet() { return locusSet; }

    private:
        Set locusSet;
};
