#include <string>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <unordered_map>

struct Key {
    int first;
    int second;
};

struct KeyHash {
    std::size_t operator()(const Key& k) const {
        return std::hash<int>()(k.first) ^
            (std::hash<int>()(k.second) << 1);
    }
};

struct KeyEqual {
    bool operator()(const Key& lhs, const Key& rhs) const {
        return lhs.first == rhs.first && lhs.second == rhs.second;
    }
};

typedef std::pair<Key, int> LocusDistPoint;
typedef std::pair<Key, double> LocusProbPoint;
typedef std::unordered_map<Key, int, KeyHash, KeyEqual> LocusDist;
typedef std::unordered_map<Key, double, KeyHash, KeyEqual> LocusProb;

class Locus {
    friend std::ostream& operator<<(std::ostream& os, const Locus& a) {
        os << "########## " << a.name << " ##########\nLocus Distribution:\n";
        for (LocusDistPoint p : a.locusDist) {
            os << "[" << p.first.first << ", " << p.first.second
                << ": " << p.second <<"], ";
        }
        os << "\nAlleles:\n";
        for (std::pair<int, int> allele : a.alleles) {
            os << "[" << allele.first << ": " << allele.second << "], ";
        }
        return os;
    }
    public:
        Locus(const std::string name) : name(name) {}
        bool operator<(const Locus& locus) const { return name < locus.name; }
        bool operator==(const Locus& locus) const { return name == locus.name; }

        std::string getName() const { return name; }
        void addLocusPeaks(Key peak) {
            locusDist[peak]++;
            alleles[peak.first]++;
            alleles[peak.second]++;
        }

        LocusProb calculateLocusProbs(double sampleSize) {
            for (LocusDistPoint point : locusDist) {
                locusProb[point.first] = (double)point.second / sampleSize;
            }
            return locusProb;
        }
    private:
        std::string name;
        LocusDist locusDist;
        std::unordered_map<int, int> alleles;
        LocusProb locusProb;

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
        int size() { return locusSet.size(); }

    private:
        Set locusSet;
};
