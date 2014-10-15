#include <string>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <boost/math/distributions/chi_squared.hpp>

struct MegaKey {
    double locus1;
    double locus2;
    std::string otherLocusName;
    double otherLocus1;
    double otherLocus2;
};

struct MegaKeyHash {
    std::size_t operator()(const MegaKey& k) const {
        return std::hash<double>() (k.locus1) ^
               (std::hash<double>() (k.locus2) << 1) ^
               (std::hash<std::string>() (k.otherLocusName) << 2) ^
               (std::hash<double>() (k.otherLocus1) << 3) ^
               (std::hash<double>() (k.otherLocus2) << 4);
    }
};

struct MegaKeyEqual {
    bool operator()(const MegaKey& lhs, const MegaKey& rhs) const {
        return lhs.locus1 == rhs.locus1 && lhs.locus2 == rhs.locus2
            && lhs.otherLocusName == rhs.otherLocusName
            && lhs.otherLocus1 == rhs.otherLocus1
            && lhs.otherLocus2 == rhs.otherLocus2;
    }
};

struct Key {
    double first;
    double second;
};

struct KeyHash {
    std::size_t operator()(const Key& k) const {
        return std::hash<double>() (k.first) ^
               (std::hash<double>() (k.second) << 1);
    }
};

struct KeyEqual {
    bool operator()(const Key& lhs, const Key& rhs) const {
        return lhs.first == rhs.first && lhs.second == rhs.second;
    }
};

typedef std::pair<Key, int> LocusDistPoint;
typedef std::pair<Key, double> LocusProbPoint;
typedef std::pair<Key, double> EVALSPoint;
typedef std::pair<double, int> AlleleDistPoint;
typedef std::pair<double, double> AlleleProbPoint;
typedef std::pair<MegaKey, int> LocusPairDistPoint;
typedef std::pair<MegaKey, double> LocusPairProbPoint;
typedef std::unordered_map<Key, int, KeyHash, KeyEqual> LocusDist;
typedef std::unordered_map<Key, double, KeyHash, KeyEqual> LocusProb;
typedef std::unordered_map<Key, double, KeyHash, KeyEqual> EVALS;
typedef std::unordered_map<MegaKey, int, MegaKeyHash, MegaKeyEqual> LocusPairDist;
typedef std::unordered_map<MegaKey, double, MegaKeyHash, MegaKeyEqual> LocusPairProb;
typedef std::unordered_map<double, int> AlleleDist;
typedef std::unordered_map<double, double> AlleleProb;


class Locus {
    friend std::ostream& operator<<(std::ostream& os, const Locus& a) {
        os << "########## " << a.name << " ##########\nLocus Distribution:\n";
        for (LocusDistPoint p : a.locusDist) {
            os << "[" << p.first.first << ", " << p.first.second
               << ": " << p.second << "], ";
        }
        os << "\nAlleles:\n";
        for (AlleleDistPoint allele : a.alleles) {
            os << "[" << allele.first << ": " << allele.second << "], ";
        }
        return os;
    }

    public:
        Locus(const std::string name) : name(name) {}

        bool operator<(const Locus& locus) const {
            return name < locus.name;
        }

        bool operator==(const Locus& locus) const {
            return name == locus.name;
        }

        bool operator!=(const Locus& locus) const {
            return name != locus.name;
        }

        std::string getName() const {
            return name;
        }

        void addLocusPeaks(Key peak) {
            locusDist[peak]++;
            alleles[peak.first]++;
            alleles[peak.second]++;
        }

        void addLocusPair(const Locus& other) {
            for (LocusDistPoint p : locusDist) {
                for (LocusDistPoint k : other.locusDist) {
                    MegaKey currKey = {
                        p.first.first,
                        p.first.second,
                        other.name,
                        k.first.first,
                        k.first.second
                    };
                    locusPairDist[currKey]++;
                }
            }
        }

        LocusProb calculateLocusProbs(double sampleSize) {
            for (LocusDistPoint point : locusDist) {
                locusProb[point.first] = (double)point.second / sampleSize;
            }
            return locusProb;
        }

        AlleleProb calculateAlleleProbs(double sampleSize) {
            for (AlleleDistPoint allele : alleles) {
                alleleProb[allele.first] = allele.second / (2 * sampleSize);
            }
            return alleleProb;
        }

        /** Finds the expected values of the locus distribution
         */
        EVALS calculateE(double sampleSize){
            for(LocusDistPoint n : locusDist) {
                double Eval;
                double population = sampleSize;
                double p = alleleProb[n.first.first];
                double q = alleleProb[n.first.second];
                if(n.first.first==n.first.second) {
                    Eval = p * q * population;
                } else {
                    Eval = 2 * p * q * population;
                }
                eVals[n.first] = Eval;
            }
            return eVals;
        }

        void isHWE(float psig) {
            double critVal = 0;
            double sub,ev,n;
            for(EVALSPoint e : eVals) {
                ev = e.second;
                n = locusDist[e.first];
                sub = n - ev;
                if(ev != 0) {
                    critVal += ((sub * sub) / ev);
                }
            }
            int alleleNum = this->getAlleleNumber();
            double df = 0.5 * (alleleNum) * (alleleNum - 1);
            boost::math::chi_squared mydist(df);
            double pval = boost::math::cdf(mydist, critVal);
            if((1 - pval) < psig) {
                std::cout << name << " is not likely to be in HWE" << std::endl;
            }
        }

        void doLinkageCompares(Locus& l, int sampleSize, float psig){
            this->calculateLocusProbs(sampleSize);
            l.calculateLocusProbs(sampleSize);
            double expression = 0;
            for(LocusProbPoint bigP1 : this->locusProb) {
                double big1 = bigP1.second;
                for(LocusProbPoint bigP2 : l.locusProb) {
                    double big2 = bigP2.second;
                    MegaKey actSeenKey = {
                        bigP1.first.first,
                        bigP1.first.second,
                        l.name,
                        bigP2.first.first,
                        bigP2.first.second
                    };
                    int actSeen = locusPairDist[actSeenKey];
                    double expectedSeen = big1 * big2 * sampleSize;
                    double square = (actSeen - expectedSeen);
                    std::cout << actSeen << " " << expectedSeen << std::endl;
                    expression += ((square * square) / expectedSeen);
                }
            }
            double alleleNum1 = this->getAlleleNumber();
            double alleleNum2 = l.getAlleleNumber();
            double sk1 = 0.5 * alleleNum1 * (alleleNum1 - 1);
            double sk2 = 0.5 * alleleNum2 * (alleleNum2 - 1);
            double df = (sk1 - 1) * (sk2 - 1);
            boost::math::chi_squared mydist(df);
            double pval = boost::math::cdf(mydist, expression);
            if ((1 - pval) < psig) {
                std::cout << this->name << " and " <<  l.name << "are likely linked " << expression << " " << df << std::endl;
            }
        }


        void printProbs() {

            std::cout << "########## " << name << " ##########\n";

            std::cout << "LocusDist" << std::endl;
            for(LocusDistPoint n : locusDist)
                std::cout << "[" << n.first.first << ", " << n.first.second
                          << ": " << n.second << "], ";

            std::cout << "\nLocusProb:" << std::endl;
            for (LocusProbPoint p : locusProb) {
                std::cout << "[" << p.first.first << ", " << p.first.second
                          << ": " << p.second << "], ";
            }

            std::cout << "\nAlleleDist:\n";
            for(AlleleDistPoint al : alleles) {
                std::cout << "[" << al.first << ", " << al.second << "], ";
            }

            std::cout << "\nAllelesProb:\n";
            double sum = 0;
            for (AlleleProbPoint allele : alleleProb) {
                std::cout << "[" << allele.first << ": " << allele.second << "], ";
                sum += allele.second;
            }
            std::cout << sum << std::endl;

            std::cout << "E:" << std::endl;
            for(EVALSPoint e : eVals) {
                std::cout << "[" << e.first.first << "," << e.first.second
                          << ":" << e.second << "], ";
            }
        }

        double getAlleleNumber(){
            double count = 0.0;
            for(AlleleProbPoint allele : alleleProb) {
                count+=1.0;
            }
            return count;
        }


    public:
        std::string name;
        LocusDist locusDist;
        AlleleDist alleles;
        LocusProb locusProb;
        AlleleProb alleleProb;
        EVALS eVals;
        EVALS eValsL;
        LocusPairDist locusPairDist;
        LocusPairProb locusPairProb;

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

        Locus operator[](int index) const {
            return locusSet.at(index);
        }

        Locus& operator[](int index) {
            return locusSet.at(index);
        }

        Set& getLocusSet() {
            return locusSet;
        }

        int size() {
            return locusSet.size();
        }

    private:
        Set locusSet;
};
