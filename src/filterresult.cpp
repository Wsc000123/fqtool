#include "filterresult.h"

FilterResult::FilterResult(Options* opt, bool paired){
    mOptions = opt;
    mPaired = paired;
    mTrimmedAdapterBases = 0;
    mTrimmedAdapterReads = 0;
    mTrimmedPolyXBases = (size_t*)std::calloc(4, sizeof(size_t));
    mTrimmedPolyXReads = (size_t*)std::calloc(4, sizeof(size_t));
    for(int i = 0; i < COMMONCONST::FILTER_RESULT_TYPES; ++i){
        mFilterReadStats[i] = 0;
    }
    mCorrectionMatrix = (size_t*)std::calloc(64, sizeof(size_t));
    mCorrectedReads = 0;
    mMergedPairs = 0;
    mSummarized = false;
}

FilterResult::~FilterResult(){
    free(mCorrectionMatrix);
    free(mTrimmedPolyXBases);
    free(mTrimmedPolyXReads);
}

void FilterResult::addFilterResult(int result){
    if(result < COMMONCONST::PASS_FILTER || result >= COMMONCONST::FILTER_RESULT_TYPES){
        return;
    }
    if(mPaired){
        mFilterReadStats[result] += 2;
    }else{
        mFilterReadStats[result] += 1;
    }
}

void FilterResult::addFilterResult(int result, int n){
    if(result < COMMONCONST::PASS_FILTER || result >= COMMONCONST::FILTER_RESULT_TYPES){
        return;
    }
    mFilterReadStats[result] += n;
}

void FilterResult::addPolyXTrimmed(int base, int length){
      mTrimmedPolyXReads[base] += 1;
      mTrimmedPolyXBases[base] += length;
}

void FilterResult::addMergedPairs(int n){
    mMergedPairs += n;
}

FilterResult* FilterResult::merge(std::vector<FilterResult*>& list){
    if(list.size() == 0){
        return NULL;
    }
    FilterResult* result = new FilterResult(list[0]->mOptions, list[0]->mPaired);
    size_t* target = result->getFilterReadStats();
    size_t* correction = result->getCorrectionMatrix();
    for(size_t i = 0; i < list.size(); ++i){
        // update mFilterReadStats 
        size_t* curStats = list[i]->getFilterReadStats();
        for(int j = 0; j < COMMONCONST::FILTER_RESULT_TYPES; ++j){
            target[j] += curStats[j];
        }
        // update mCorrectionMatrix
        size_t* curCorr = list[i]->getCorrectionMatrix();
        for(int p = 0; p < 64; ++p){
            correction[p] += curCorr[p];
        }
        // update mTrimmedAdapterReads/Bases
        result->mTrimmedAdapterReads += list[i]->mTrimmedAdapterReads;
        result->mTrimmedAdapterBases += list[i]->mTrimmedAdapterBases;
        result->mMergedPairs += list[i]->mMergedPairs;

        // update mTrimmedPolyXBases/Reads
        for(int b = 0; b < 4; ++b){
            result->mTrimmedPolyXBases[b] += list[i]->mTrimmedPolyXBases[b];
            result->mTrimmedPolyXReads[b] += list[i]->mTrimmedPolyXReads[b];
        }

        // update read counting
        result->mCorrectedReads += list[i]->mCorrectedReads;

        // merge adapter stats
        for(auto& e: list[i]->mAdapter1Count){
            if(result->mAdapter1Count.count(e.first) > 0){
                result->mAdapter1Count[e.first] += e.second;
            }else{
                result->mAdapter1Count[e.first] = e.second;
            }
        }

        for(auto& e: list[i]->mAdapter2Count){
            if(result->mAdapter2Count.count(e.first) > 0){
                result->mAdapter2Count[e.first] += e.second;
            }else{
                result->mAdapter2Count[e.first] = e.second;
            }
        }
    }
    return result;
}

void FilterResult::summary(bool force){
    if(mSummarized && !force){
        return;
    }
    mCorrectedBases = 0;
    for(int p = 0; p < 64; ++p){
        mCorrectedBases += mCorrectionMatrix[p];
    }
    mSummarized = true;
}

size_t FilterResult::getTotalCorrectedBases(){
    if(!mSummarized){
        summary();
    }
    return mCorrectedBases;
}

void FilterResult::addCorrection(char from, char to){
    int f = from & 0x07;
    int t = to & 0x07;
    ++mCorrectionMatrix[f * 8 +t];
}

void FilterResult::incCorrectedReads(int count){
    mCorrectedReads += count;
}

size_t FilterResult::getCorrectionNum(char from, char to){
    int f = from & 0x07;
    int t = to & 0x07;
    return mCorrectionMatrix[f * 8 + t];
}

void FilterResult::addAdapterTrimmed(const std::string& adapter, bool isR2){
    if(adapter.empty()){
        return;
    }
    ++mTrimmedAdapterReads;
    mTrimmedAdapterBases += adapter.length();
    if(!isR2){
        if(mAdapter1Count.count(adapter) > 0){
            ++mAdapter1Count[adapter];
        }else{
            mAdapter1Count[adapter] = 1;
        }
    }else{
        if(mAdapter2Count.count(adapter) > 0){
            ++mAdapter2Count[adapter];
        }else{
            mAdapter2Count[adapter] = 1;
        }
    }
}

void FilterResult::addAdapterTrimmed(const std::string& adapterR1, const std::string& adapterR2){
    mTrimmedAdapterReads += 2;
    mTrimmedAdapterBases += adapterR1.length() + adapterR2.length();
    if(!adapterR1.empty()){
        if(mAdapter1Count.count(adapterR1) > 0){
            ++mAdapter1Count[adapterR1];
        }else{
            mAdapter1Count[adapterR1] = 1;
        }
    }

    if(!adapterR2.empty()){
        if(mAdapter2Count.count(adapterR2) > 0){
            ++mAdapter2Count[adapterR2];
        }else{
            mAdapter2Count[adapterR2] = 1;
        }
    }
}

std::ostream& operator<<(std::ostream& os, FilterResult* re){
    Options* mOptions = re->mOptions;
    os << "reads passed filter: " << re->mFilterReadStats[COMMONCONST::PASS_FILTER] << "\n";
    os << "reads failed due to low quality: " << re->mFilterReadStats[COMMONCONST::FAIL_QUALITY] << "\n";
    os << "reads failed due to too many N: " << re->mFilterReadStats[COMMONCONST::FAIL_N_BASE] << "\n";
    if(mOptions->lengthFilter.enabled){
        os << "reads failed due to too short: " << re->mFilterReadStats[COMMONCONST::FAIL_LENGTH] << "\n";
        if(mOptions->lengthFilter.maxReadLength > 0){
            os << "reads failed due to too long: " << re->mFilterReadStats[COMMONCONST::FAIL_TOO_LONG] << "\n";
        }
    }
    if(mOptions->complexityFilter.enabled){
        os << "reads failed due to low complexity: " << re->mFilterReadStats[COMMONCONST::FAIL_COMPLEXITY] << "\n";
    }
    if(mOptions->adapter.enableTriming){
        os << "reads with adapter trimmed: " << re->mTrimmedAdapterReads << "\n";
        os << "bases trimmed due to adapters: " << re->mTrimmedAdapterBases << "\n";
    }
    if(mOptions->correction.enabled){
        os << "reads corrected by overlap analysis: " << re->mCorrectedReads << "\n";
        os << "bases corrected by overlap analysis: " << re->getTotalCorrectedBases() << "\n";
    }
    return os;
}

void FilterResult::reportJsonBasic(jsn::json& j){
    j["passed_filter_reads"] = mFilterReadStats[COMMONCONST::PASS_FILTER];
    j["low_quality_reads"] = mFilterReadStats[COMMONCONST::FAIL_QUALITY];
    j["too_many_N_reads"] = mFilterReadStats[COMMONCONST::FAIL_N_BASE];
    if(mOptions->correction.enabled){
        j["corrected_reads"] = mCorrectedReads;
        j["corrected_bases"] = getTotalCorrectedBases();
    }
    if(mOptions->complexityFilter.enabled){
        j["low_complexity_reads"] = mFilterReadStats[COMMONCONST::FAIL_COMPLEXITY];
    }
    if(mOptions->lengthFilter.enabled){
        j["too_short_reads"] = mFilterReadStats[COMMONCONST::FAIL_LENGTH];
        if(mOptions->lengthFilter.maxReadLength > 0){
            j["too_long_reads"] = mFilterReadStats[COMMONCONST::FAIL_TOO_LONG];
        }
    }
}

void FilterResult::reportHtmlBasic(std::ofstream& ofs, size_t totalReads, size_t totalBases){
    ofs << "<table class='summary_table'>\n";
    htmlutil::outputTableRow(ofs, "reads passed filters:", htmlutil::formatNumber(mFilterReadStats[COMMONCONST::PASS_FILTER]) + " (" + std::to_string(mFilterReadStats[COMMONCONST::PASS_FILTER] * 100.0 / totalBases) + "%)");
    htmlutil::outputTableRow(ofs, "low_quality_reads", htmlutil::formatNumber(mFilterReadStats[COMMONCONST::FAIL_QUALITY]) + " (" + std::to_string(mFilterReadStats[COMMONCONST::FAIL_QUALITY] * 100.0 / totalBases) + "%)");
    htmlutil::outputTableRow(ofs, "too_many_N_reads", htmlutil::formatNumber(mFilterReadStats[COMMONCONST::FAIL_N_BASE]) + " (" + std::to_string(mFilterReadStats[COMMONCONST::FAIL_N_BASE] * 100.0 / totalBases) + "%)");
    if(mOptions->correction.enabled){
        htmlutil::outputTableRow(ofs, "corrected_reads", htmlutil::formatNumber(mCorrectedReads) + " (" + std::to_string(mCorrectedReads * 100.0 / totalReads) + "%)");
        htmlutil::outputTableRow(ofs, "corrected_bases", htmlutil::formatNumber(getTotalCorrectedBases()) + " (" + std::to_string(getTotalCorrectedBases()) + " (" + std::to_string(getTotalCorrectedBases() * 100.0 / totalBases) + "%)");
    }
    if(mOptions->complexityFilter.enabled){
        htmlutil::outputTableRow(ofs, "low_complexity_reads", htmlutil::formatNumber(mFilterReadStats[COMMONCONST::FAIL_COMPLEXITY]) + " (" + std::to_string(mFilterReadStats[COMMONCONST::FAIL_COMPLEXITY] * 100.0 / totalReads) + "%)");
    }
    if(mOptions->lengthFilter.enabled){
        htmlutil::outputTableRow(ofs, "too_short_reads", htmlutil::formatNumber(mFilterReadStats[COMMONCONST::FAIL_LENGTH]) + " (" + std::to_string(mFilterReadStats[COMMONCONST::FAIL_LENGTH] * 100.0 / totalReads) + "%)");
        if(mOptions->lengthFilter.maxReadLength > 0){
            htmlutil::outputTableRow(ofs, "too_long_reads", htmlutil::formatNumber(mFilterReadStats[COMMONCONST::FAIL_TOO_LONG]) + " (" + std::to_string(mFilterReadStats[COMMONCONST::FAIL_TOO_LONG] * 100.0 /totalReads) + "%)");
        }
    }
    ofs << "</table>\n";
}

void FilterResult::reportAdaptersJsonDetails(jsn::json& j, std::map<std::string, size_t>& adapterCounts){
    size_t totalAdapters = 0;
    for(auto& e: adapterCounts){
        totalAdapters += e.second;
    }
    if(totalAdapters == 0){
        return;
    }
    const double dTotalAdapters = (double)totalAdapters;
    size_t reported = 0;
    for(auto& e: adapterCounts){
        if(e.second / dTotalAdapters < mOptions->adapter.reportThreshold){
            continue;
        }
        j[e.first] = e.second;
        reported += e.second;
    }
    size_t unreported = totalAdapters - reported;
    if(unreported > 0){
        j["others"] = unreported;
    }
}

void FilterResult::reportAdaptersHtmlDetails(std::ofstream& ofs, std::map<std::string, size_t>& adapterCounts, size_t totalBases){
    size_t totalAdapters = 0;
    size_t totalAdaptersBases = 0;
    for(auto& e: adapterCounts){
        totalAdapters += e.second;
        totalAdaptersBases += e.first.length();
    }
    double frac = (double)totalAdaptersBases / (double)totalBases;
    if(mPaired){
        frac *= 2.0;
    }
    if(frac < 0.01){
        ofs << "<div class='sub_section_tips'>The input has little adapter percentage (~" << std::to_string(frac * 100.0) << "%), probably it's trimmed before.</div>\n";
    }
    if(totalAdapters == 0){
        return;
    }

    ofs << "<table class='summary_table'>\n";
    ofs << "<tr><td class='adapter_col' style='font-size:14px;color:#ffffff;background:#556699'>" << "Sequence" << "</td><td class='col2' style='font-size:14px;color:#ffffff;background:#556699'>" << "Occurrences" << "</td    ></tr>\n";    
    const double dTotalAdapters = (double)totalAdapters;
    size_t reported = 0;
    for(auto& e: adapterCounts){
        if(e.second / dTotalAdapters < mOptions->adapter.reportThreshold){
            continue;
        }
        htmlutil::outputTableRow(ofs, e.first, e.second);
        reported += e.second;
    }
    size_t unreported = totalAdapters - reported;
    if(unreported > 0){
        std::string tag = "other adapter sequences";
        if(reported == 0){
            tag = "all adapter sequences";
        }
        htmlutil::outputTableRow(ofs, tag, unreported);
    }
    ofs << "</table>\n";
}

void FilterResult::reportAdaptersJsonSummary(jsn::json& j){
    j["adapter_trimmed_reads"] = mTrimmedAdapterReads;
    j["adapter_trimmed_bases"] = mTrimmedAdapterBases;
    j["read1_adapter_sequence"] = mOptions->adapter.adapterSeqR1Provided ? mOptions->adapter.inputAdapterSeqR1 : mOptions->adapter.detectedAdapterSeqR1;
    if(mPaired){
        j["read2_adapter_sequence"] = mOptions->adapter.adapterSeqR2Provided ? mOptions->adapter.inputAdapterSeqR2 : mOptions->adapter.detectedAdapterSeqR2;
    }
    jsn::json jR1AdapterCount;
    reportAdaptersJsonDetails(jR1AdapterCount, mAdapter1Count);
    j["read1_adapter_counts"] = jR1AdapterCount;
    if(mPaired){
        jsn::json jR2AdapterCount;
        reportAdaptersJsonDetails(jR2AdapterCount, mAdapter2Count);
        j["read2_adapter_counts"] = jR2AdapterCount;
    }
}

void FilterResult::reportAdaptersHtmlSummary(std::ofstream& ofs, size_t totalBases){
    ofs << "<div class='subsection_title' onclick=showOrHide('read1_adapters')>Adapter or bad ligation of read1</div>\n";
    ofs << "<div id='read1_adapters'>\n";
    reportAdaptersHtmlDetails(ofs, mAdapter1Count, totalBases);
    ofs << "</div>\n";
    if(mPaired){
        ofs << "<div class='subsection_title' onclick=showOrHide('read2_adapters')>Adapter or bad ligation of read2</div>\n";
        ofs << "<div id='read2_adapters'>\n";
        reportAdaptersHtmlDetails(ofs, mAdapter2Count, totalBases);
        ofs << "</div>\n";
    }
}

void FilterResult::reportPolyXTrimJson(jsn::json& j){
    const char atcg[4] = {'A', 'T', 'C', 'G'};
    j["total_polyx_trimmed_reads"] = std::accumulate(mTrimmedPolyXReads, mTrimmedPolyXReads + 4, 0);
    jsn::json jPolyReads;
    for(int b = 0; b < 4; ++b){
        jPolyReads[atcg[b]] = mTrimmedPolyXReads[b];
    }
    j["polyx_trimmed_reads"] = jPolyReads;
    j["total_polyx_trimmed_bases"] = std::accumulate(mTrimmedPolyXBases, mTrimmedPolyXBases + 4, 0);
    jsn::json jPolyBases;
    for(int b = 0; b < 4; ++b){
        jPolyBases[atcg[b]] = mTrimmedPolyXBases[b];
    }
    j["polyx_trimmed_bases"] = jPolyBases;
}
