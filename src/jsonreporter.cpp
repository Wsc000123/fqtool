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
    // pre filtering qc Summary
    jsn::json jPreFilterQC;
    jPreFilterQC["TotalReads"] = preTotalReads;
    jPreFilterQC["TotalBases"] = preTotalBases;
    jPreFilterQC["Q20Bases"] = preQ20Bases;
    jPreFilterQC["Q30Bases"] = preQ30Bases;
    jPreFilterQC["Q20BaseRate"] = preQ20Rate;
    jPreFilterQC["Q30BaseRate"] = preQ30Rate;
    jPreFilterQC["Read1Length"] = preRead1Length;
    if(mOptions->isPaired()){
        jPreFilterQC["Read2Length"] = preRead2Length;
    }
    jPreFilterQC["GCRate"] = preGCRate;
    jReport["Summary"]["BeforeFiltering"] = jPreFilterQC;

    // after filter qc Summary
    jsn::json jPostFilterQC;
    jPostFilterQC["TotalReads"] = postTotalReads;
    jPostFilterQC["TotalBases"] = postTotalBases;
    jPostFilterQC["Q20Bases"] = postQ20Bases;
    jPostFilterQC["Q30Bases"] = postQ30Bases;
    jPostFilterQC["Q20BaseRate"] = postQ20Rate;
    jPostFilterQC["Q30BaseRate"] = postQ30Rate;
    jPostFilterQC["Read1Length"] = postRead1Length;
    if(mOptions->isPaired()){
        jPostFilterQC["Read2Length"] = postRead2Length;
    }
    jPostFilterQC["GCRate"] = postGCRate;
    jReport["Summary"]["AfterFiltering"] = jPostFilterQC;

    // filter result
    jsn::json jFilterResult;
    fresult->reportJsonBasic(jFilterResult);
    jReport["FilterResult"] = jFilterResult;

    // duplication result
    if(mOptions->duplicate.enabled){
        jsn::json jDupResult;
        jDupResult["Rate"] = mDupRate;
        std::vector<int32_t> dupVec(mDupHist, mDupHist + mOptions->duplicate.histSize);
        jDupResult["Histogram"] = dupVec;
        std::vector<double> gcVec(mDupMeanGC, mDupMeanGC + mOptions->duplicate.histSize);
        jDupResult["MeanGC"] = gcVec;
        jReport["Duplication"] = jDupResult;
    }

    // inset size result
    if(mOptions->isPaired()){
        jsn::json jInsert;
        jInsert["Peak"] = mInsertSizePeak;
        jInsert["Unknown"] = mInsertHist[mOptions->insertSizeMax];
        std::vector<int32_t> insVec(mInsertHist, mInsertHist + mOptions->insertSizeMax);
        jInsert["Histogram"] = insVec;
        jReport["InsertSize"] = jInsert;
    }

    // adapter trimming result
    if(mOptions->adapter.enableTriming){
        jsn::json jAdapterTrim;
        fresult->reportAdaptersJsonSummary(jAdapterTrim);
        jReport["AdapterTrim"] = jAdapterTrim;
    }

    // polyx trimming result
    if(fresult && (mOptions->polyXTrim.enabled || mOptions->polyGTrim.enabled)){
        jsn::json jPolyXTrim;
        fresult->reportPolyXTrimJson(jPolyXTrim);
        jReport["PolyxTrimming"] = jPolyXTrim;
    }
    // read1 before filtering
    if(preStats1){
        jReport["Read1BeforeFiltering"] = preStats1->reportJson();
    }
    // read2 before filtering
    if(preStats2){
        jReport["Read2BeforeFiltering"] = preStats2->reportJson();
    }
    // read1 after filtering
    if(postStats1){
        std::string name = "Read1AfterFiltering";
        if(mOptions->mergePE.enabled){
            name = "MergedAndFiltered";
        }
        jReport[name] = postStats1->reportJson();
    }
    // read2 after filtering
    if(postStats2 && !mOptions->mergePE.enabled){
        jReport["Read2AfterFiltering"] = postStats2->reportJson();
    }
    ofs << jReport.dump(4);
    ofs.close();
}
