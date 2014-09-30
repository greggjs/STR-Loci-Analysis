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
typedef std::pair<int, int> AlleleDistPoint;
typedef std::pair<int, double> AlleleProbPoint;
typedef std::unordered_map<Key, int, KeyHash, KeyEqual> LocusDist;
typedef std::unordered_map<Key, double, KeyHash, KeyEqual> LocusProb;
typedef std::unordered_map<int, int> AlleleDist;
typedef std::unordered_map<int, double> AlleleProb;

class Locus {
    friend std::ostream& operator<<(std::ostream& os, const Locus& a) {
        os << "########## " << a.name << " ##########\nLocus Distribution:\n";
        for (LocusDistPoint p : a.locusDist) {
            os << "[" << p.first.first << ", " << p.first.second
                << ": " << p.second <<"], ";
        }
        os << "\nAlleles:\n";
        for (AlleleDistPoint allele : a.alleles) {
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

        AlleleProb calculateAlleleProbs(double sampleSize) {
            int highestAllelePeak = this->getHighestAllelePeak();
            for (AlleleDistPoint allele : alleles) {
                Key currPoint = { allele.first, allele.first };
                double sum = 0;
                for (int i = allele.first + 1; i <= highestAllelePeak; i++) {
                    Key allelePoint = {allele.first, i};
                    sum += locusProb[allelePoint];
                }
                alleleProb[allele.first] = locusProb[currPoint] + sum;
            }
            return alleleProb;
        }

        void printProbs() {
            std::cout << "########## " << name << " ##########\nLocus Distribution:\n";
            for (LocusProbPoint p : locusProb) {
                std::cout << "[" << p.first.first << ", " << p.first.second
                    << ": " << p.second <<"], ";
            }
            for (AlleleProbPoint allele : alleleProb) {
                std::cout << "[" << allele.first << ": " << allele.second << "], ";
            }
            std::cout << std::endl;
        }
    private:
        std::string name;
        LocusDist locusDist;
        AlleleDist alleles;
        LocusProb locusProb;
        AlleleProb alleleProb;

        int getHighestAllelePeak() {
            std::vector<AlleleDistPoint> alleleVector(alleles.begin(), alleles.end());
            std::sort(alleleVector.begin(), alleleVector.end(),
                [](AlleleDistPoint& a, AlleleDistPoint& b){ return a.first > b.first; });
            return alleleVector.at(0).first;
        }

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
