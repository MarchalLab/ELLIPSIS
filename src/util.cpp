
#include "util.h"
#include "spliceGraph.h"
#include "GTFReader.h"

//#include <math.h>
#include <cmath>
#include <sys/stat.h>
#include <iostream>
#include <set>
#include <algorithm>
#include <cassert>
#include <dirent.h>
#include <queue>
#include <map>
#include <experimental/filesystem>

using namespace std;

//--------------------------//
// thread-safe progress bar //
//--------------------------//

Progressbar::Progressbar(size_t numTot, string message, ostream& os) : os(os), message(move(message)), numDone(0),
                                                                              numTot(numTot){}

void Progressbar::setMessage(const string& mess) {
    message = mess;
}

void Progressbar::addProgress(int done) {
    // lock the mutex
    unique_lock<mutex> lock(m);

    //set progress value
    numDone += done;
}

float Progressbar::getProgress(){
    // lock the mutex
    unique_lock<mutex> lock(m);

    float value = float(numDone) / float(numTot);
    return value;
}

void Progressbar::printProgress(){
    // lock the mutex
    unique_lock<mutex> lock(m);

    float value = float(numDone) / float(numTot);

    //move cursor to beginning of line
    os << "\r";

    os << message << " \t[";
    int pos = int(float(barwidth) * value);
    for (int i = 0; i < barwidth; ++i) {
        if (i < pos) os << "=";
        else if (i == pos) os << ">";
        else os << " ";
    }
    os << "] " << int(value * 100.0) << " %";

    os.flush();
}

void Progressbar::update(int done){
    int prev = int(getProgress() * 100);

    addProgress(done);

    int newValue = int(getProgress() * 100);

    if (newValue/5 > prev/5){ //only print progress every 5% progress made
        printProgress();
    }
}

//------------------------//
// Thread-safe log writer //
//------------------------//

void LogWriter::init(const string& fileName){
    // lock the mutex
    unique_lock<mutex> lock(m);

    ofs = ofstream(fileName);
}

void LogWriter::writeLine(const string& line) {
    // lock the mutex
    unique_lock<mutex> lock(m);

    ofs << line << endl;
}

void LogWriter::close() {
    // lock the mutex
    unique_lock<mutex> lock(m);

    ofs.close();
}

//-------------------//
// Utility functions //
//-------------------//

int Util::currentTimer = 0;
chrono::time_point<chrono::system_clock> Util::startTime[MAX_TIMERS];

bool Util::fileExists(const std::string& filename) {
    std::ifstream file(filename.c_str(), std::ios::in);
    bool OK = file.good();
    file.close();
    return OK;
}

bool Util::dirExists(const std::string& path) {
    struct stat buffer;
    if (stat(path.c_str(), & buffer) == 0){
        return true;
    }
    return false;
}

void Util::createDir(const std::string& dir){
    if (mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0){
        if (errno == EEXIST){
            //throw runtime_error("Directory " + dir + " already exists");
            return;
        }
        if (errno == EACCES | errno == EROFS){
            throw runtime_error("No permission to write " + dir);
        }
        if (errno == ENOENT | errno == ENOTDIR){
            throw runtime_error("The parent directory of " + dir + " does not exist.");
        }
        if (errno == ENOSPC){
            throw runtime_error("No space left on device to create " + dir);
        }
        throw runtime_error(dir + " could not be created.");
    }
}

void Util::removeDir(const std::string& dir){
    std::experimental::filesystem::remove_all(dir);
}

vector<string> Util::getFileList(const string& dir){

    if (! dirExists(dir)){
        throw runtime_error(dir  +" does not exist.");
    }

    vector<string> files;

    DIR* d = opendir(dir.c_str());
    struct dirent * tp;
    while((tp = readdir(d)) != nullptr){
        string file(tp->d_name);
        if (file != "." and file != "..") {
            files.push_back(file);
        }
    }
    closedir(d);

    return files;
}

vector<string> Util::getFileList(const string& dir, const regex& pattern){

    if (! dirExists(dir)){
        throw runtime_error(dir  +" does not exist.");
    }

    vector<string> files;

    DIR* dirp = opendir(dir.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != nullptr) {
        string file(dp->d_name);
        if (regex_match(file, pattern)){
            files.push_back(file);
        }
    }
    closedir(dirp);

    return files;
}

int Util::getNumFiles(const string& dir){

    if (! dirExists(dir)){
        throw runtime_error(dir  +" does not exist.");
    }

    int numFiles = 0;

    DIR* dirp = opendir(dir.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != nullptr) {
        string file(dp->d_name);
        if (file != "." and file != "..") {
            numFiles ++;
        }
    }
    closedir(dirp);

    return numFiles;
}

vector<vector<string>> Util::getNonOverlappingSets(vector<string> v1, vector<string> v2){

    //remove double values in vectors & sort
    set<string> v1_set(v1.begin(), v1.end());
    set<string> v2_set(v2.begin(), v2.end());

    vector<string> unique_v1(v1_set.size());
    vector<string> unique_v2(v2_set.size());
    vector<string> intersect(min(v1.size(), v2.size()));

    vector<string>::iterator ls_v1 = set_difference(v1_set.begin(), v1_set.end(), v2_set.begin(), v2_set.end(), unique_v1.begin());
    vector<string>::iterator ls_v2 = set_difference(v2_set.begin(), v2_set.end(), v1_set.begin(), v1_set.end(), unique_v2.begin());
    vector<string>::iterator ls_intersect = set_intersection(v1_set.begin(), v1_set.end(), v2_set.begin(), v2_set.end(), intersect.begin());

    unique_v1.erase(ls_v1, unique_v1.end());
    unique_v2.erase(ls_v2, unique_v2.end());
    intersect.erase(ls_intersect, intersect.end());

    vector<vector<string>> ret{unique_v1, unique_v2, intersect};
    return ret;
}

std::vector<std::string> Util::splitString(const std::string &s, char delim){
    std::vector<std::string> elems;

    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

bool Util::compareCaseInsensitiveStrings(std::string s1, std::string s2){
    transform(s1.begin(), s1.end(), s1.begin(), ::tolower);
    transform(s2.begin(), s2.end(), s2.begin(), ::tolower);

    if(s1 == s2)
        return true; //s1 == s2
    return false; //s1 != s2

}

string Util::humanTime (double time)
{
    uint64_t timeInt = uint64_t(time);

    uint64_t days = timeInt / 86400;
    timeInt -= days * 86400;
    uint64_t hours = timeInt / 3600;
    timeInt -= hours * 3600;
    uint64_t min = timeInt / 60;
    timeInt -= min * 60;
    uint64_t sec = timeInt;
    uint64_t ms = uint64_t((time - timeInt) * 1000);

    ostringstream iss;
    if (days > 0) {
        iss << days << "d " << hours << "h " << min << "min";
    } else if (hours > 0) {
        iss << hours << "h " << min << "min " << sec << "s";
    } else if (min > 0) {
        iss << min << "min " << sec << "s";
    } else if (sec > 0) {
        iss << sec << "s " << ms << "ms";
    } else {
        iss << ms << "ms";
    }

    return iss.str();
}

void Util::startChrono()
{
    // make sure we don't use too many timers
    assert(currentTimer < MAX_TIMERS); //Too many timers set

    startTime[currentTimer] = chrono::system_clock::now();
    currentTimer++;
}

double Util::stopChrono()
{
    // make sure stopChrono isn't called too often
    currentTimer--;
    assert(currentTimer >= 0); //Stopped too many timers.

    chrono::duration<double> elapsed = chrono::system_clock::now() - startTime[currentTimer];
    return (elapsed.count());
}

std::string Util::stopChronoStr() {
    return humanTime(stopChrono());
}

int Util::min3(int a, int b, int c){
    return(min(a, b), c);
}

double Util::poissonPDF(unsigned int k, double mu)
{
    return exp(k*log(mu)-mu-lgamma(k+1));
}

double Util::logPoissonPDF(unsigned int k, double mu)
{
    return k*log(mu)-mu-lgamma(k+1);
}

double Util::normalPDF(double x, double mu, double var){
    return 1.0/(sqrt(2*var*M_PI)) * exp(-1.0/(2*var) * pow((x-mu),2));
}

double Util::logNormalPDF(double x, double mu, double var){
    if (x == 0.0){
        double ret  = log(normalPDF(x,mu,var) + normalCDF(x, mu, var));
        //make sure no -inf values
        if (isinf(ret)){
            ret = std::numeric_limits<int>::min();
        }
        return ret;
    }
    return -0.5*log(2*var*M_PI) - 0.5 * pow(x-mu,2)/var;
}

double Util::normalCDF(double x, double mu, double var){
    return 0.5 * (1 + erf( (x-mu)/sqrt(2*var) ));
}

double Util::logNormalProb(double x, double mu, double var){
    double prob;
    if (x-0.5 < 0){
        prob = normalCDF(ceil(x+0.5), mu, var);
    }else{
        prob = normalCDF(x+0.5, mu, var) - normalCDF(x-0.5, mu, var);
    }
    assert(0 <= prob & prob <= 1);

    //avoid -inf/nan values
    if (prob == 0){
        double lowLogProb = logNormalPDF(x, mu, var);
        assert(lowLogProb < 0 && ! isinf(lowLogProb));
        return(lowLogProb);
    }

    return log(prob);

}

vector<vector<double>> Util::transpose(vector<vector<double>>& A){

    size_t nRow = A.size();
    assert(nRow > 0);
    size_t nCol = A[0].size();
    assert(nCol > 0);

    //initialize vector
    vector<vector<double>> At(nCol, vector<double>(nRow));

    for (int i = 0; i < nRow; i++){
        for (int j = 0; j < nCol; j++){
            At[j][i] = A[i][j];
        }
    }

    return At;
}

vector<vector<double>> Util::matrixProd(vector<vector<double>>& A, vector<vector<double>>& B){

    size_t nRow_A = A.size(), nRow_B = B.size();
    assert(nRow_A > 0); assert(nRow_B > 0);
    size_t nCol_A = A[0].size(), nCol_B = B[0].size();
    assert(nCol_A > 0); assert(nCol_B > 0);

    assert(nCol_A == nRow_B);

    //initialize result
    vector<vector<double>> prod(nRow_A, vector<double>(nCol_B, 0.0));

    for (int i = 0; i < nRow_A; i++){
        for (int j = 0; j < nCol_B; j++){
            for (int k = 0; k < nCol_A; k++){
                prod[i][j]  += A[i][k] * B[k][j];
            }
        }
    }

    return prod;
}

vector<double> Util::matrixVectorProd(vector<vector<double>> &A, vector<double> b) {
    size_t nRow_A = A.size(), nRow_b = b.size();
    assert(nRow_A > 0); assert(nRow_b > 0);
    size_t nCol_A = A[0].size();
    assert(nCol_A > 0);

    assert(nCol_A == nRow_b);

    //initialize result
    vector<double> prod(nRow_A, 0.0);

    for (int i = 0; i < nRow_A; i++){
        for (int k = 0; k < nCol_A; k++){
            prod[i] += A[i][k] * b[k];
        }
    }

    return prod;
}

vector<vector<double>> Util::invert(vector<vector<double>> &A) {

    size_t nRow= A.size();
    assert(nRow > 0);

    //matrix is only invertible if it is square
    assert(nRow == A[0].size());

    //extend A with unit matrix
    vector<vector<double>> Aext(nRow, vector<double>(2*nRow, 0.0));
    for (int i = 0; i < nRow; i++) {
        //copy elements from A
        for (int j = 0; j < nRow; j++) {
            Aext[i][j] = A[i][j];
        }

        //set 1 on diagonal of unit matrix next to A
        Aext[i][i + nRow] = 1.0;
    }

    //use Gaussian elimination
    for (int i = 0; i < nRow; i++){

        double pivot = Aext[i][i];
        assert(pivot > 0); //make sure matrix is invertible

        //normalize row according to pivot
        for (int j = i; j < 2*nRow; j++){
            Aext[i][j] /= pivot;
        }

        //elimination
        for (int ii = 0; ii < nRow; ii++){
            if (ii != i){
                double p = Aext[ii][i];
                for (int j = 0; j < 2*nRow; j++){
                    Aext[ii][j] -= p * Aext[i][j];
                }
            }
        }

    }


    //get inverse
    vector<vector<double>> Ainv(nRow, vector<double>(nRow));
    for (int i = 0; i < nRow; i++){
        for (int j = 0; j < nRow; j++){
            Ainv[i][j] = Aext[i][j+nRow];
        }
    }

    return Ainv;
}

vector<double> Util::solveOLS(vector<vector<double>>& X, vector<double> y){

    //transpose X
    vector<vector<double>> Xt = Util::transpose(X);

    //compute XtX and Xty
    vector<vector<double>> XtX = Util::matrixProd(Xt, X);
    vector<double> Xty = Util::matrixVectorProd(Xt, y);

    //invert XtX
    vector<vector<double>> XtX_inv = Util::invert(XtX);

    //multiply inv(XtX) * Xty
    return Util::matrixVectorProd(XtX_inv, Xty);

}

void Util::choleskyLinSolve(std::vector<double>& A,
                            std::vector<double>& b)
{
    size_t m = b.size();
    assert(A.size() == m*m);

    // lower triangulate
    for(int i = 0; i < m; i++){
        for (int j = 0; j < i; j++){
            for(int k = 0; k < j; k++)
                A[m*i + j] -= (A[m*i+k]*A[m*j+k]);
            A[m*i+j] /= A[m*j+j];
        }
        for (int k = 0; k < i; k++)
            A[m*i+i] -= (A[m*i+k]*A[m*i+k]);
        assert(A[m*i+i] > 0);
        A[m*i+i] = sqrt(A[m*i+i]);
    }

    // forward substitution
    for(size_t i = 0; i < m; i++){
        for (size_t j = 0; j < i; j++)
            b[i] -= (A[m*i+j] * b[j]);
        b[i] /= A[m*i+i];
    }
    // backward substitution
    for(size_t i = m-1; i >= 0; --i){ //start at m-1
        for(size_t j = i+1; j < m; j++)
            b[i] -= (A[m*j+i] * b[j]);
        b[i] /= A[m*i+i];
    }
}


void Util::getTruePSI(const string& gtfFile, const string& graphDir, const string& clusterFlowDir, const string& clusterFile, const string& outDir) {

    //get mapping from clusters to cells
    map<string, vector<string>> clusterMap = readClusterFile(clusterFile);

    //mapping from node/edge to list of transcripts in which they are included
    map<SpliceGraph*, map<SGNode*, vector<string>>> nodeMapping;
    map<SpliceGraph*, map<SGEdge*, vector<string>>> edgeMapping;

    //get transcript mapping
    GTFReader gtfReader(gtfFile);
    gtfReader.getTranscriptMapping(graphDir, nodeMapping, edgeMapping);

    //compute clusterflow per gene
    for (const auto& it : nodeMapping){

        SpliceGraph* graph = it.first;
        map<SGNode*, vector<string>> nodeMap = it.second;
        map<SGEdge*, vector<string>> edgeMap = edgeMapping[it.first];
        string geneID = graph->getGeneID();

        //read clusterFlow for this gene
        map<string, map<string, double>> clusterFlow = readClusterFlowFile(clusterFlowDir + "/" + geneID + ".tsv");

        //write PSI per cell to file
        ofstream out(outDir + "/" + geneID);

        //write header
        vector<SGNode*> orderedNodes; vector<SGEdge*> orderedEdges;
        out << "cell" << graph->getVarNames(orderedNodes, orderedEdges) << endl;

        for (const auto& clusterIt : clusterMap){

            string cluster = clusterIt.first;

            for (const auto& cell : clusterIt.second){

                //write cellName to file
                out << cell;

                //compute cell PSI for each node
                for (auto node : orderedNodes){
                    vector<string> transcripts = nodeMapping[graph][node];
                    double nodePSI = 0;
                    for (const auto& t : transcripts)
                        nodePSI += clusterFlow[cluster][t];
                    out << "\t" << round(nodePSI* 100) * 0.01;
                }

                //compute cell PSI for each edge
                for (auto edge: orderedEdges){
                    vector<string> transcripts = edgeMapping[graph][edge];
                    double edgePSI = 0;
                    for (const auto& t : transcripts)
                        edgePSI += clusterFlow[cluster][t];
                    out << "\t" << round(edgePSI*100) * 0.01;
                }

                out << endl;
            }
        }

        out.close();

        graph->deleteGraph();
        delete(graph);
    }
}

map<string, vector<string>> Util::readClusterFile(const string& clusterFile){

    ifstream IFS(clusterFile);
    if (!IFS){
        throw ios_base::failure("Cannot open file " + clusterFile);
    }

    string line, cell, cluster;
    map<string, vector<string>> clusterMap;

    while (getline(IFS, line)){
        istringstream iss(line);
        iss >> cell >> cluster;
        clusterMap[cluster].push_back(cell);
    }
    IFS.close();

    return clusterMap;
}

map<string, map<string, double>> Util::readClusterFlowFile(const string& clusterFlowFile){

    map<string, map<string, double>> clusterFlow;

    fstream clusterFlowStream(clusterFlowFile);
    string line, cluster, transcript;
    double flow;
    while(getline(clusterFlowStream, line)){
        istringstream iss(line);
        iss >> cluster >> transcript >> flow;
        clusterFlow[cluster][transcript] = flow;
    }
    clusterFlowStream.close();

    return clusterFlow;
}

map<string, vector<string>> Util::readNeighborFile(const string& neighborFile, const string& alignmentDir){

    ifstream IFS(neighborFile);
    if (!IFS){
        throw ios_base::failure("Cannot open file " + neighborFile);
    }

    string line, cell, neighbor;
    map<string, vector<string>> neighborMap;

    while (getline(IFS, line)){
        istringstream iss(line);
        iss >> cell;

        //check if cell exists in alignmentDir
        if (not Util::fileExists(alignmentDir + "/" + cell))
            throw runtime_error( "The cells in " + neighborFile + " do not correspond to an existing bam file in " +
                                 alignmentDir + "please check your if the header corresponds with the bam filenames (including .bam)");

        while ((iss >> neighbor)){

            //check if cell exists in alignmentDir
            if (not Util::fileExists(alignmentDir + "/" + cell))
                throw runtime_error( "The cells in " + neighborFile + " do not correspond to an existing bam file in " +
                                     alignmentDir + "please check your if the header corresponds with the bam filenames (including .bam)");

            //add neighbor
            neighborMap[cell].push_back(neighbor);
        }

    }

    return neighborMap;

}

map<string, double> Util::readLowQualFile(const string& lowQualFile){

    map<string, double> lowQualFrac;

    fstream IFS(lowQualFile);
    string line, cell;
    double frac;

    //skip first line (header)
    getline(IFS, line);

    while(getline(IFS, line)){
        istringstream iss(line);
        iss >> cell >> frac;
        lowQualFrac[cell] = frac;
    }

    IFS.close();

    return lowQualFrac;

}

size_t Util::string2uint(const char *s){

    char* err_ptr;

    long ret = strtol(s, &err_ptr, 10);

    if (ret < 0)
        throw invalid_argument("value cannot be negative");

    if (*err_ptr != '\0')
        throw invalid_argument( "received non integer value" );

    if (ret == numeric_limits<size_t>::max() | ret == numeric_limits<long>::max())
        throw invalid_argument("value too big");

    return size_t(ret);

}