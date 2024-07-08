
#ifndef ELLIPSIS_SJCOUNTER_H
#define ELLIPSIS_SJCOUNTER_H

#include <mutex>
#include <map>

using namespace std;

class SJCounter {
private:
    mutex countMutex;
    map<string, map<string, double>> countTable;
public:
    SJCounter() = default;

    /**
     * Get SJCounter from file
     * @param file
     */
    explicit SJCounter(const string& file);

    void addCount(const string& spliceJunction, const string& cell, double count);

    map<string, map<string, double>> getCounts();

    /**
     * Write junction count table to file
     * @param outFile
     * @return true if count table not empty
     */
    bool writeToFile(const string& outFile);

    size_t getNumJunctions();

    /**
     * reset to empty count table
     */
    void clear();

    /**
     * merge 2 SJcounters (with no overlapping cells)
     * @param SJcounter2
     */
    void merge(const SJCounter& SJcounter2);
};


#endif //ELLIPSIS_SJCOUNTER_H
