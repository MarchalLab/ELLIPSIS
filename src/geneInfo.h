
#ifndef ELLIPSIS_GENEINFO_H
#define ELLIPSIS_GENEINFO_H

#include <string>
#include <vector>

using namespace std;

class GeneInfo {
private:
    string geneID;
    string chrName;
    size_t start;
    size_t end;
    string strand;

public:
    explicit GeneInfo(const string& info);

    string getGeneID() const;

    string getChrName();

    size_t getStart() const;

    size_t getEnd() const;

    string getStrand();

    /**
     * Check if positions of junction are part of an annotated gene.
     * @param j_chr
     * @param j_start
     * @param j_end
     * @param j_strand
     * @return True if gene positions include junction positions
     */
    bool containsJunction(const string& j_chr, size_t j_start, size_t j_end, const string& j_strand) const;

    /**
     * Check if positions of junction are part of an annotated gene regardless of strand.
     * @param j_chr
     * @param j_start
     * @param j_end
     * @return True if gene positions include junction positions
     */
    bool containsJunction(const string& j_chr, size_t j_start, size_t j_end) const;

};


#endif //ELLIPSIS_GENEINFO_H
