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
    long preTotalReads = preStats1->getReads();
    long preTotalBases = preStats1->getBases();
    long preQ20Bases = preStats1->getQ20();
    long preQ30Bases = preStats1->getQ30();
    long preTotalGC = preStats1->getGCNumber();
    int preRead1Length = preStats1->getMeanLength();
    int preRead2Length = 0;
    long postTotalReads = postStats1->getReads();
    long postTotalBases = postStats1->getBases();
    long postQ20Bases = postStats1->getQ20();
    long postQ30Bases = postStats1->getQ30();
    long postTotalGC = postStats1->getGCNumber();
    int postRead1Length = postStats1->getMeanLength();
    int postRead2Length = 0;
    
    if(preStats2 && postStats2){
        preTotalReads += preStats2->getReads();
        preTotalBases += preStats2->getBases();
        preQ20Bases += preStats2->getQ20();
        preQ30Bases += preStats2->getQ30();
        preTotalGC += preStats2->getGCNumber();
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

    // whole report
    jsn::json jReport;
    // pre filtering qc summary
    jsn::json jPreFilterQC;
    jPreFilterQC["total_reads"] = preTotalReads;
    jPreFilterQC["total_bases"] = preTotalBases;
    jPreFilterQC["Q20_bases"] = preQ20Bases;
    jPreFilterQC["Q30_bases"] = preQ30Bases;
    jPreFilterQC["Q20_rate"] = preQ20Rate;
    jPreFilterQC["Q30_rate"] = preQ30Rate;
    jPreFilterQC["read1_mean_length"] = preRead1Length;
    if(mOptions->isPaired()){
        jPreFilterQC["read2_mean_length"] = preRead2Length;
    }
    jPreFilterQC["gc_content"] = preGCRate;
    jReport["summary"]["before_filtering"] = jPreFilterQC;

    // after filter qc summary
    jsn::json jPostFilterQC;
    jPostFilterQC["total_reads"] = postTotalReads;
    jPostFilterQC["total_bases"] = postTotalBases;
    jPostFilterQC["Q20_bases"] = postQ20Bases;
    jPostFilterQC["Q30_bases"] = postQ30Bases;
    jPostFilterQC["Q20_rate"] = postQ20Rate;
    jPostFilterQC["Q30_rate"] = postQ30Rate;
    jPostFilterQC["read1_mean_length"] = postRead1Length;
    if(mOptions->isPaired()){
        jPostFilterQC["read2_mean_length"] = postRead2Length;
    }
    jPostFilterQC["gc_content"] = postGCRate;
    jReport["summary"]["after_filtering"] = jPostFilterQC;

    // filter result
    jsn::json jFilterResult;
    fresult->reportJsonBasic(jFilterResult);
    jReport["filtering_result"] = jFilterResult;

    // duplication result
    if(mOptions->duplicate.enabled){
        jsn::json jDupResult;
        jDupResult["rate"] = mDupRate;
        std::vector<int32_t> dupVec(mDupHist, mDupHist + mOptions->duplicate.histSize);
        jDupResult["histogram"] = dupVec;
        std::vector<double> gcVec(mDupMeanGC, mDupMeanGC + mOptions->duplicate.histSize);
        jDupResult["mean_gc"] = gcVec;
        jReport["duplication"] = jDupResult;
    }

    // inset size result
    if(mOptions->isPaired()){
        jsn::json jInsert;
        jInsert["peak"] = mInsertSizePeak;
        jInsert["unknown"] = mInsertHist[mOptions->insertSizeMax];
        std::vector<int32_t> insVec(mInsertHist, mInsertHist + mOptions->insertSizeMax);
        jInsert["histogram"] = insVec;
        jReport["insert_size"] = jInsert;
    }

    // adapter trimming result
    if(mOptions->adapter.enableTriming){
        jsn::json jAdapterTrim;
        fresult->reportAdaptersJsonSummary(jAdapterTrim);
        jReport["adapter_trim"] = jAdapterTrim;
    }

    // polyx trimming result
    if(fresult && (mOptions->polyXTrim.enabled || mOptions->polyGTrim.enabled)){
        jsn::json jPolyXTrim;
        fresult->reportPolyXTrimJson(jPolyXTrim);
        jReport["polyx_trimming"] = jPolyXTrim;
    }
    // read1 before filtering
    if(preStats1){
        jReport["read1_before_filtering"] = preStats1->reportJson();
    }
    // read2 before filtering
    if(preStats2){
        jReport["read2_before_filtering"] = preStats2->reportJson();
    }
    // read1 after filtering
    if(postStats1){
        std::string name = "read1_after_filtering";
        if(mOptions->mergePE.enabled){
            name = "merged_and_filtered";
        }
        jReport[name] = postStats1->reportJson();
    }
    // read2 after filtering
    if(postStats2 && !mOptions->mergePE.enabled){
        jReport["read2_after_filtering"] = postStats2->reportJson();
    }
    ofs << jReport;
    ofs.close();
}
