
#include <fstream>
#include <set>
#include <vector>
#include <sstream>
#include <cassert>
#include "SJCounter.h"

SJCounter::SJCounter(const string& file){
    //open file
    fstream IFS = fstream(file.c_str());
    if (!IFS)
        throw ios_base::failure("Cannot open file " + file);

    string line;

    //read header (cells with novel junctions)
    vector<string> cells;
    string cell;
    getline(IFS, line);
    istringstream issHeader(line);
    while (issHeader >> cell){
        cells.push_back(cell);
    }

    //read junctionCounts
    string SJ; double count;
    while(getline(IFS, line)){
        istringstream iss(line);
        iss >> SJ;
        for (const auto& thisCell : cells){
            iss >> count;
            if (count > 0)
                addCount(SJ, thisCell, count);
        }
    }

    IFS.close();
}

void SJCounter::addCount(const string& spliceJunction, const string& cell, double count){
    lock_guard<mutex> lock (countMutex);
    assert(countTable.find(spliceJunction) == countTable.end() | countTable[spliceJunction].find(cell) == countTable[spliceJunction].end());
    countTable[spliceJunction][cell] = count;
}

map<string, map<string, double>>  SJCounter::getCounts(){
    lock_guard<mutex> lock(countMutex);
    return countTable;
}

bool SJCounter::writeToFile(const string& outFile) {

    //do not output anything if no novel junctions are observed
    if (countTable.empty()){
        return false;
    }

    // get all cells
    set<string> cells;
    for (auto const &jIt: countTable){
        for (auto const &cIt : jIt.second){
            cells.insert(cIt.first);
        }
    }

    //open file
    ofstream OFS(outFile);

    //write header
    for (auto const &cell : cells){
        OFS << "\t" << cell;
    }
    OFS << endl;

    //write counts for each junction
    for (auto const &jIt : countTable){
        string junction = jIt.first;
        OFS << junction;
        for (auto const &cell : cells){
            auto cIt = jIt.second.find(cell);
            if (cIt == jIt.second.end())
                OFS << "\t" << 0;
            else
                OFS << "\t" << cIt->second;
        }
        OFS << endl;
    }

    //close file
    OFS.close();

    return true;
}

size_t SJCounter::getNumJunctions() {
    return countTable.size();
}

void SJCounter::clear() {
    lock_guard<mutex> lock (countMutex);
    countTable.clear();
}

void SJCounter::merge(const SJCounter& SJcounter2){
    for (auto const &jIt : SJcounter2.countTable) {
        for (auto const &cIt: jIt.second) {
            this->addCount(jIt.first, cIt.first, cIt.second);
        }
    }
}
