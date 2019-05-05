#include "jsonreporter.h"

JsonReporter::JsonReporter(Options* opt){
    mOptions = opt;
    mDupHist = NULL;
    mDupRate = 0;
}

JsonReporter::~JsonReporter(){
}

void JsonReporter::setDupHist(size_t* dupHist, double* dupMeanGC, double dupRate){
    mDupHist = dupHist;
    mDupMeanGC = dupMeanGC;
    mDupRate = dupRate;
}

void JsonReporter::setInsertHist(long* insertHist, int insertSizePeak){
    mInsertHist = insertHist;
    mInsertSizePeak = insertSizePeak;
}

void JsonReporter::report(FilterResult* fresult, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2){
    std::ofstream ofs(mOptions->jsonFile);
    ofs << "{" << std::endl;
    long preTotalReads = preStats1->getReads();
    long preTotalBases = preStats1->getBases();
    long preQ20Bases = preStats1->getQ20();
    long preQ30Bases = preStats1->getQ30();
    long preTotalGC = preStats1->getGCNumber();
    int preRead1Length = preStats1->getMeanLength();
    int preRead2Length = 0;
    long postTotalReads = preStats1->getReads();
    long postTotalBases = preStats1->getBases();
    long postQ20Bases = preStats1->getQ20();
    long postQ30Bases = preStats1->getQ30();
    long postTotalGC = preStats1->getGCNumber();
    int postRead1Length = postStats1->getMeanLength();
    int postRead2Length = 0;
    
    if(preStats2 && postStats2){
        preTotalReads += postStats2->getReads();
        preTotalBases += postStats2->getBases();
        preQ20Bases += postStats2->getQ20();
        preQ30Bases += postStats2->getQ30();
        preTotalGC += postStats2->getGCNumber();
        postTotalReads += postStats2->getReads();
        postTotalBases += postStats2->getBases();
        postQ20Bases += postStats2->getQ20();
        postQ30Bases += postStats2->getQ30();
        postTotalGC += postStats2->getGCNumber();
        preRead2Length = preStats2->getMeanLength();
        postRead2Length = postStats2->getMeanLength();
    }

    double preQ20Rate = preTotalBases == 0 ? 0.0 : (double)preQ20Bases / preTotalBases;
    double preQ30Rate = preTotalBases == 0 ? 0.0 : (double)preQ30Bases / preTotalBases;
    double preGCRate = preTotalBases == 0 ? 0.0 : (double)preTotalGC / preTotalBases;
    double postQ20Rate = postTotalBases == 0 ? 0.0 : (double)postQ20Bases / postTotalBases;
    double postQ30Rate = postTotalBases == 0 ? 0.0 : (double)postQ30Bases / postTotalBases;
    double postGCRate = postTotalBases == 0 ? 0.0 : (double)postTotalGC / postTotalBases;

    ofs << "\t" << "\"summary\": {" << std::endl;
    ofs << "\t\t" << "\"before_filtering\": {" << std::endl;
    jsonutil::writeRecord(ofs, "\t\t", "total_reads", preTotalReads);
    jsonutil::writeRecord(ofs, "\t\t", "total_bases", preTotalBases);
    jsonutil::writeRecord(ofs, "\t\t", "Q20_bases", preQ20Bases);
    jsonutil::writeRecord(ofs, "\t\t", "Q30_bases", preQ30Bases);
    jsonutil::writeRecord(ofs, "\t\t", "Q20_rate", preQ20Rate);
    jsonutil::writeRecord(ofs, "\t\t", "Q30_rate", preQ30Rate);
    jsonutil::writeRecord(ofs, "\t\t", "read1_mean_length", preRead1Length);
    if(mOptions->isPaired()){
         jsonutil::writeRecord(ofs, "\t\t", "read2_mean_length", preRead2Length);
    }
    jsonutil::writeRecord(ofs, "\t\t", "gc_content", preGCRate);
    ofs << "\t\t" << "}," << std::endl;

    ofs << "\t\t" << "\"after_filtering\": {" << std::endl;
    jsonutil::writeRecord(ofs, "\t\t", "total_reads", postTotalReads);
    jsonutil::writeRecord(ofs, "\t\t", "total_bases", postTotalBases);
    jsonutil::writeRecord(ofs, "\t\t", "Q20_bases", postQ20Bases);
    jsonutil::writeRecord(ofs, "\t\t", "Q30_bases", postQ30Bases);
    jsonutil::writeRecord(ofs, "\t\t", "Q20_rate", postQ20Rate);
    jsonutil::writeRecord(ofs, "\t\t", "Q30_rate", postQ30Rate);
    jsonutil::writeRecord(ofs, "\t\t", "read1_mean_length", postRead1Length);
    if(mOptions->isPaired()){
         jsonutil::writeRecord(ofs, "\t\t", "read2_mean_length", postRead2Length);
    }
    jsonutil::writeRecord(ofs, "\t\t", "gc_content", postGCRate);
    ofs << "\t\t" << "}," << std::endl;

    if(fresult){
        ofs << "\t\t" << "\"filtering_result\": ";
        fresult->reportJsonBasic(ofs, "\t\t");
    }
    ofs << "\t" << "}," << std::endl;

    if(mOptions->duplicate.enabled){
        ofs << "\t" << "\"duplication\": {" << std::endl;
        ofs << "\t\t\"rate\": " << mDupRate << "," << std::endl;
        ofs << "\t\t\"histogram\": [";
        for(int d = 1; d < mOptions->duplicate.histSize; ++d){
            ofs << mDupHist[d];
            if(d != mOptions->duplicate.histSize - 1){
                ofs << ",";
            }
        }
        ofs << "]," << std::endl;
        ofs << "\t\t\"mean_gc\": [";
        for(int d = 1; d < mOptions->duplicate.histSize; ++d){
            ofs << mDupMeanGC[d];
            if(d != mOptions->duplicate.histSize - 1){
                ofs << ",";
            }
        }
        ofs << "]" << std::endl;
        ofs << "\t" << "}";
        ofs << "," << std::endl;
    }

    if(mOptions->isPaired()){
        ofs << "\t" << "\"insert_size\": {" << std::endl;
        ofs << "\t\t\"peak\": " << mInsertSizePeak << "," << std::endl;
        ofs << "\t\t\"unknown\": " << mInsertHist[mOptions->insertSizeMax] << "," << std::endl;
        ofs << "\t\t\"histogram\": [";
        for(int d = 0; d < mOptions->insertSizeMax; ++d){
            ofs << mInsertHist[d];
            if(d!=mOptions->insertSizeMax-1){
                ofs << ",";
            }
        }
        ofs << "]" << std::endl;
        ofs << "\t" << "}";
        ofs << "," << std::endl;
    }

    if(mOptions->adapter.enableTriming){
        ofs << "\t" << "\"adapter_trim\": ";
        fresult->reportAdaptersJsonSummary(ofs, "\t");
    }

    if(fresult && (mOptions->polyXTrim.enabled || mOptions->polyGTrim.enabled)){
        ofs << "\t" << "\"polyx_trimming\": ";
        fresult->reportPolyXTrimJson(ofs, "\t");
    }

    ofs << "}" << std::endl;
}
