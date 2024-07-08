
#include <iostream>
#include <sstream>

#include "geneInfo.h"

using namespace std;

GeneInfo::GeneInfo(const string& info){
    istringstream iss(info);
    iss >> geneID >> chrName >> start >> end >> strand;
}

string GeneInfo::getGeneID() const{
    return geneID;
}

string GeneInfo::getChrName(){
    return chrName;
}

size_t GeneInfo::getStart() const{
    return start;
}

size_t GeneInfo::getEnd() const{
    return end;
}

string GeneInfo::getStrand() {
    return strand;
}

bool GeneInfo::containsJunction(const string& j_chr, size_t j_start, size_t j_end, const string& j_strand) const{
    if (chrName != j_chr) return false;
    if (strand != j_strand) return false;
    if (j_start < start) return false;
    if (j_end > end) return false;
    return true;
}

bool GeneInfo::containsJunction(const string& j_chr, size_t j_start, size_t j_end) const{
    if (chrName != j_chr) return false;
    if (j_start < start) return false;
    if (j_end > end) return false;
    return true;
}