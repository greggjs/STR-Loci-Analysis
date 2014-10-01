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
typedef std::pair<Key,double> EVALSPoint;
typedef std::pair<double, int> AlleleDistPoint;
typedef std::pair<double, double> AlleleProbPoint;
typedef std::unordered_map<Key, int, KeyHash, KeyEqual> LocusDist;
typedef std::unordered_map<Key, double, KeyHash, KeyEqual> LocusProb;
typedef std::unordered_map<Key,double, KeyHash, KeyEqual> EVALS;
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

        std::string getName() const {
            return name;
        }

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

        void doLinkageCompares(Locus& l,int sampleSize, float psig){
            double expression = 0;
            this->calculateLocusProbs(sampleSize);
            l.calculateLocusProbs(sampleSize);
            for(LocusProbPoint bigP1 : this->locusProb) {
                int actualNumSeen1 = this->locusDist[bigP1.first];
                double big1 = bigP1.second;
                for(LocusProbPoint bigP2 : l.locusProb) {
                    double big2 = bigP2.second;
                    int actualNumSeen2 = l.locusDist[bigP2.first];
                    int actSeen = actualNumSeen1 + actualNumSeen2;
                    double expectedSeen = big1 * big2 * sampleSize;
                    double square = (actSeen - expectedSeen);
                    expression += ((square * square) / expectedSeen);
                }
            }
            int alleleNum1 = this->getAlleleNumber();
            int alleleNum2 = l.getAlleleNumber();
            int sk1 = 0.5 * alleleNum1 * (alleleNum1 - 1);
            int sk2 = 0.5 * alleleNum2 * (alleleNum2 - 1);
            int df = (sk1 - 1) * (sk2 - 1);
            boost::math::chi_squared mydist(df);
            double pval = boost::math::cdf(mydist, expression);
            if ((1 - pval) < psig) {
                std::cout << this->name << " and " <<  l.name << "are likely linked" << std::endl;
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

        int getAlleleNumber(){
            int count = 0;
            for(AlleleProbPoint allele : alleleProb) {
                count+=1;
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
