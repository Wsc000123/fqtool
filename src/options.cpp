#include "options.h"

Options::Options(){
    version = "0.0.0";
    compile = std::string(__TIME__) + " " + std::string(__DATE__);
    in1 = "";
    in2 = "";
    out1 = "";
    out2 = "";
    reportTitle = "Fastq Report";
    thread = 1;
    compression = 2;
    phred64 = false;
    donotOverwrite = false;
    inputFromSTDIN = false;
    outputToSTDOUT = false;
    readsToProcess = 0;
    interleavedInput = false;
    insertSizeMax = 512;
    overlapDiffLimit = 5;
    verbose = false;
    jsonFile = "report.json";
    htmlFile = "report.html";
}

void Options::update(){
    adapter.adapterSeqR1Provided = adapter.inputAdapterSeqR1.empty() ? false : true;
    adapter.adapterSeqR2Provided = adapter.inputAdapterSeqR2.empty() ? false : true;
    if(indexFilter.enabled){
        initIndexFilter(indexFilter.index1File, indexFilter.index2File, indexFilter.threshold);
    }
}

bool Options::isPaired(){
    return in2.length() > 0 || interleavedInput;
}

bool Options::adapterCutEnabled(){
    if(adapter.enableTriming){
        if(isPaired() || !adapter.inputAdapterSeqR1.empty()){
            return true;
        }
    }
    return false;
}

bool Options::validate(){
    return true;
}

bool Options::shallDetectAdapter(bool isR2){
    if(!adapter.enableTriming){
        return false;
    }
    if(isR2){
        return isPaired() && adapter.enableDetectForPE && adapter.inputAdapterSeqR2 == "auto";
    }else{
        if(isPaired()){
            return adapter.enableDetectForPE && adapter.inputAdapterSeqR1 == "auto";
        }else{
            return adapter.inputAdapterSeqR1 == "auto";
        }
    }
}

void Options::initIndexFilter(const std::string& blacklistFile1, const std::string& blacklistFile2, int threshold){
    if(blacklistFile1.empty() && blacklistFile2.empty()){
        return;
    }
    if(!blacklistFile1.empty()){
        util::valid_file(blacklistFile1);
        indexFilter.blacklist1 = makeListFromFileByLine(blacklistFile1);
    }
    if(!blacklistFile2.empty()){
        util::valid_file(blacklistFile2);
        indexFilter.blacklist2 = makeListFromFileByLine(blacklistFile2);
    }
    if(indexFilter.blacklist1.empty() && indexFilter.blacklist2.empty()){
        return;
    }
    indexFilter.enabled = true;
    indexFilter.threshold = threshold;
}

std::vector<std::string> Options::makeListFromFileByLine(const std::string& filename){
    std::vector<std::string> ret;
    std::ifstream fr(filename);
    std::string line;
    while(std::getline(fr, line)){
        util::strip(line);
        if(line.find_first_not_of("ATCG") != std::string::npos){
            util::error_exit("processing " + filename + ", each line should be one index, which can only contain A/T/C/G");
        }
        ret.push_back(line);
    }
    return ret;
}

std::string Options::getAdapter1(){
    if(adapter.inputAdapterSeqR1 == "" || adapter.inputAdapterSeqR1 == "auto"){
        return "unspecified";
    }else{
        return adapter.inputAdapterSeqR1;
    }
}

std::string Options::getAdapter2(){
    if(adapter.inputAdapterSeqR2 == "" || adapter.inputAdapterSeqR2 == "auto"){
        return "unspecified";
    }else{
        return adapter.inputAdapterSeqR2;
    }
}
