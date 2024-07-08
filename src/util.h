
#ifndef ELLIPSIS_UTIL_H
#define ELLIPSIS_UTIL_H

#include "geneInfo.h"
#include "sgedge.h"
#include "spliceGraph.h"

#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <chrono>
#include <mutex>
#include <atomic>
#include <map>
#include <regex>


//-----------------------------//
// multi-threaded progress bar //
//-----------------------------//

class Progressbar{
private:
    size_t barwidth = 70;
    size_t numDone;
    size_t numTot;
    mutex m;
    string message;
    ostream &os;

public:

    /**
     * Constructor
     * @param numTot total number of jobs
     * @param message job description
     * @param os output stream where progress can be printed
     */
    Progressbar(size_t numTot, string message, ostream& os);

    /**
     * Set message of progressbar
     * @param message to be printed in front of progressbar
     */
    void setMessage(const string& message);

    /**
     * Add to progress
     * @param done number of jobs done
     */
    void addProgress(int done);

    /**
     * Get progress value
     * @return progress value between 0 and 1
     */
    float getProgress();

    /**
     * Print progressbar to output
     */
    void printProgress();

    /**
     * Update progress value and print progress
     * @param done number of jobs done
     */
    void update(int done);

};

//------------------------//
// Thread-safe log writer //
//------------------------//

class LogWriter{
private:
    mutex m;
    ofstream ofs;
public:
    LogWriter() = default;
    void init(const string& fileName);
    void writeLine(const string& line);
    void close();
};

//-------------------//
// Utility functions //
//-------------------//

#define MAX_TIMERS 16

class Util {
private:
    static double prevProgress;
    static int currentTimer;
    static std::chrono::time_point<std::chrono::system_clock> startTime[MAX_TIMERS];

public:
    /**
     * Check whether a file exists
     * @param filename
     * @return True of false
     */
    static bool fileExists(const std::string& filename);

    /**
     * Check whether a directory exists
     * @param path
     * @return True of false
     */
    static bool dirExists(const std::string& path);

    /**
     * Create a new directory
     * @param dir
     */
    static void createDir(const std::string& dir);

    /**
     * remove a directory and its contents
     * @param dir
     */
    static void removeDir(const std::string& dir);

    /**
     * Get all files in directory
     * @param dir
     * @return list of files in directory
     */
    static vector<string> getFileList(const string& dir);

    /**
     * Get all files in directory that correspond to pattern
     * @param dir
     * @param pattern
     * @return
     */
    static vector<string> getFileList(const string& dir, const regex& pattern);

    /**
     * Get number of files in directory
     * @param dir
     * @return
     */
    static int getNumFiles(const string& dir);

    /**
     * Get 3 non-overlapping sets
     * @param v1
     * @param v2
     * @return [v1-v2, v2-v1, v1 intersect v2]
     */
    static vector<vector<string>> getNonOverlappingSets(vector<string> v1, vector<string> v2);

    /**
     * Split a string by delimiter
     * @param s
     * @param delim
     * @return vector containing partial strings
     */
    static std::vector<std::string> splitString(const std::string &s, char delim);

    /**
     * Compare if 2 strings are equal (case insensitive)
     * @param s1
     * @param s2
     * @return
     */
    static bool compareCaseInsensitiveStrings(std::string s1, std::string s2);

    /**
     * Create a string with a human readable version of a time period
     * @param time Time period (expressed in s)
     * @return String with a human readable version of a time period
     */
    static std::string humanTime(double time);

    /**
     * Start a chronometer
     */
    static void startChrono();

    /**
     * Stop the chronometer
     * @return The time in
     */
    static double stopChrono();

    /**
     * Stop the chronometer and return a human readable string
     * @return A human readable string containg the elapsed time
     */
    static std::string stopChronoStr();

    /**
     * Get minimum of 3 integer values
     * @param a
     * @param b
     * @param c
     * @return minimum
     */
    static int min3(int a, int b, int c);

    /**
     * Transpose matrix A
     * @param A
     * @return transposed matrix At
     */
    static vector<vector<double>> transpose(vector<vector<double>>& A);

    /**
     * compute the matrix product AB
     * @param A
     * @param B
     * @return AB
     */
    static vector<vector<double>> matrixProd(vector<vector<double>>& A, vector<vector<double>>& B);

    /**
     * compute the matrix vector product Ab
     * @param A
     * @param b
     * @return Ab
     */
    static vector<double> matrixVectorProd(vector<vector<double>>& A, vector<double> b);

    /**
     * Compute the inverse of matrix A (using the Gaussian elimination method)
     * @param A
     * @return A^(-1)
     */
    static vector<vector<double>> invert(vector<vector<double>>& A);

    /**
     * solve ordinary least squares equation Xa = y
     * @param X row sorted matrix
     * @param y
     * @return a for which (Xa - y)² is minimized
     */
    static vector<double> solveOLS(vector<vector<double>>& X, vector<double> y);

    /**
     * Use a Cholesky decomposition to solve the linear system of equations Ax=b
     * A and b will be overwritten such that A contains the lower triangular L
     * and b contains the solution x
     * @param A an mxm square SPD matrix
     * @param b a 1xm vector (contains solution x after this routine)
     */
    static void choleskyLinSolve(std::vector<double>& A,
                                 std::vector<double>& b);

    /**
     * Compute the probability p(k) from a Poisson distribution with mean mu
     * @param k Number of observations
     * @param mu Average
     * @return The probability p(k)
     */
    static double poissonPDF(unsigned int k, double mu);

    /**
     * Compute the logprob p(k) from a Poisson distribution with mean mu
     * @param k Number of observations
     * @param mu Average
     * @return The probability p(k)
     */
    static double logPoissonPDF(unsigned int k, double mu);

    /**
     * Compute the probability density function f(x) of a normal distribution with mean mu and standard deviation var
     * @param x
     * @param mu
     * @param var sigma^2
     * @return
     */
    static double normalPDF(double x, double mu, double var);

    /**
     * Compute the logprob density function log[f(x)] of a normal distribution with mean mu and standard deviation var
     * @param x
     * @param mu
     * @param var sigma^2
     * @return
     */
    static double logNormalPDF(double x, double mu, double var);

    /**
     * Compute the cumulative distribution function of a normal distribution with mean mu and standard deviation var
     * @param x
     * @param mu
     * @param var sigma^2
     * @return
     */
    static double normalCDF(double x, double mu, double var);

    /**
     * The probability for the interval between x-0.5 and x for a normal distribution with mean mu and variance var
     * (with all negative values mapped to 0)
     * @param x
     * @param mu
     * @param var sigma ^ 2
     * @return
     */
    static double logNormalProb(double x, double mu, double var);

    /**
     * get true PSI values per gene and per cell
     * @param gtfFile GTF file containing only transcripts that are expressed in the ground truth
     * @param graphDir directory containing ground truth splicegraphs
     * @param clusterFlowDir directory containing PSI values per transcript
     * @param clusterFile cluster annotation file
     * @param outDir directory to write ground truth PSI values to
     */
    static void getTruePSI(const string& gtfFile, const string& graphDir, const string& clusterFlowDir, const string& clusterFile, const string& outDir);

    /**
     * Read in cluster annotation
     * @param clusterFile
     * @return map containing list of cells per cluster
     */
    static map<string, vector<string>> readClusterFile(const string& clusterFile);

    /**
     * read cluster flow file
     * @param clusterFlowFile
     * @return map : cluster -> (transcript -> inclusion)
     */
    static map<string, map<string, double>> readClusterFlowFile(const string& clusterFlowFile);

    /**
     * read in neighbor file
     * @param neighborFile
     * @param alignmentDir
     * @return
     */
    static map<string, vector<string>>readNeighborFile(const string& neighborFile, const string& alignmentDir);

    /**
     * read fractions of low quality mapping reads
     * @param lowQualFile
     * @return
     */
    static map<string, double> readLowQualFile(const string& lowQualFile);

    /**
     * convert string to unsigned int (with under/overflow detection)
     * @param s
     * @return
     */
    static size_t string2uint(const char *s);
};


#endif //ELLIPSIS_UTIL_H
