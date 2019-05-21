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
    j["PassedFilterReads"] = mFilterReadStats[COMMONCONST::PASS_FILTER];
    j["LowQualityReads"] = mFilterReadStats[COMMONCONST::FAIL_QUALITY];
    j["TooManyNReads"] = mFilterReadStats[COMMONCONST::FAIL_N_BASE];
    if(mOptions->correction.enabled){
        j["CorrectedReads"] = mCorrectedReads;
        j["CorrectedBases"] = getTotalCorrectedBases();
    }
    if(mOptions->complexityFilter.enabled){
        j["LowComplexityReads"] = mFilterReadStats[COMMONCONST::FAIL_COMPLEXITY];
    }
    if(mOptions->lengthFilter.enabled){
        j["TooShortReads"] = mFilterReadStats[COMMONCONST::FAIL_LENGTH];
        if(mOptions->lengthFilter.maxReadLength > 0){
            j["TooLongReads"] = mFilterReadStats[COMMONCONST::FAIL_TOO_LONG];
        }
    }
}

CTML::Node FilterResult::reportHtmlBasic(size_t totalReads, size_t totalBases){
    CTML::Node table("table.summary_table");
    table.AppendChild(htmlutil::make2ColRowNode("Reads Passed Filters", std::to_string(mFilterReadStats[COMMONCONST::PASS_FILTER]) + "(" + std::to_string(mFilterReadStats[COMMONCONST::PASS_FILTER] * 100.0 / totalBases) + "%)")); 
    table.AppendChild(htmlutil::make2ColRowNode("Low Quality Reads", std::to_string(mFilterReadStats[COMMONCONST::FAIL_QUALITY]) + "(" + std::to_string(mFilterReadStats[COMMONCONST::FAIL_QUALITY] * 100.0 /totalBases) + "%)"));
    table.AppendChild(htmlutil::make2ColRowNode("Too Many N Reads", std::to_string(mFilterReadStats[COMMONCONST::FAIL_N_BASE]) + "(" + std::to_string(mFilterReadStats[COMMONCONST::FAIL_N_BASE] * 100.0 / totalBases) + "%)"));
    if(mOptions->correction.enabled){
        table.AppendChild(htmlutil::make2ColRowNode("Corrected Reads", std::to_string(mCorrectedReads) + "(" + std::to_string(mCorrectedReads * 100.0 / totalReads) + "%)"));
        table.AppendChild(htmlutil::make2ColRowNode("Corrected Bases", std::to_string(getTotalCorrectedBases()) + "(" + std::to_string(getTotalCorrectedBases() * 100.0 / totalBases)  + "%)"));
    }
    if(mOptions->complexityFilter.enabled){
        table.AppendChild(htmlutil::make2ColRowNode("Low Complexity Reads", std::to_string(mFilterReadStats[COMMONCONST::FAIL_COMPLEXITY]) + "(" + std::to_string(mFilterReadStats[COMMONCONST::FAIL_COMPLEXITY] * 100.0 / totalReads) + "%)"));
    }
    if(mOptions->lengthFilter.enabled){
        table.AppendChild(htmlutil::make2ColRowNode("Too Short Reads", std::to_string(mFilterReadStats[COMMONCONST::FAIL_LENGTH]) + "(" + std::to_string(mFilterReadStats[COMMONCONST::FAIL_LENGTH] * 100.0 / totalReads) + "%)"));
        if(mOptions->lengthFilter.maxReadLength > 0){
            table.AppendChild(htmlutil::make2ColRowNode("Too Long Reads", std::to_string(mFilterReadStats[COMMONCONST::FAIL_TOO_LONG]) + "(" + std::to_string(mFilterReadStats[COMMONCONST::FAIL_TOO_LONG] * 100.0 /totalReads) + "%)"));
        }
    }
    return table;
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
        j["Others"] = unreported;
    }
}

CTML::Node FilterResult::reportAdaptersHtmlDetails(std::map<std::string, size_t>& adapterCounts, size_t totalBases){
    CTML::Node table("table.summary_table");
    CTML::Node tableHeadRow("tr");
    CTML::Node tableHeadCol1("td.adapter_col", "Sequence");
    tableHeadCol1.SetAttribute("style", "font-size:14px;color:#ffffff;background:#556699");
    CTML::Node tableHeadCol2("td.col2", "Occurences");
    tableHeadCol2.SetAttribute("style", "font-size:14px;color:#ffffff;background:#556699");
    tableHeadRow.AppendChild(tableHeadCol1).AppendChild(tableHeadCol2);
    table.AppendChild(tableHeadRow);
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
    if(totalAdapters == 0){
        return table;
    }
    const double dTotalAdapters = (double)totalAdapters;
    size_t reported = 0;
    for(auto& e: adapterCounts){
        if(e.second / dTotalAdapters < mOptions->adapter.reportThreshold){
            continue;
        }
        CTML::Node row("tr");
        row.AppendChild(CTML::Node("td.adapter_col", e.first));
        row.AppendChild(CTML::Node("td.col2", std::to_string(e.second) + "(" + std::to_string(e.second*100.0/dTotalAdapters) + "%)"));
        table.AppendChild(row);
        reported += e.second;
    }
    size_t unreported = totalAdapters - reported;
    if(unreported > 0){
        std::string tag = "other adapter sequences";
        if(reported == 0){
            tag = "all adapter sequences";
        }
        table.AppendChild(htmlutil::make2ColRowNode(tag, std::to_string(unreported) + "(" + std::to_string(unreported*100.0/dTotalAdapters) + "%)"));
    }
    return table;
}

void FilterResult::reportAdaptersJsonSummary(jsn::json& j){
    j["AdapterTrimmedReads"] = mTrimmedAdapterReads;
    j["AdapterTrimmedBases"] = mTrimmedAdapterBases;
    j["Read1AdapterSequence"] = mOptions->adapter.adapterSeqR1Provided ? mOptions->adapter.inputAdapterSeqR1 : mOptions->adapter.detectedAdapterSeqR1;
    if(mPaired){
        j["Read2AdapterSequence"] = mOptions->adapter.adapterSeqR2Provided ? mOptions->adapter.inputAdapterSeqR2 : mOptions->adapter.detectedAdapterSeqR2;
    }
    jsn::json jR1AdapterCount;
    reportAdaptersJsonDetails(jR1AdapterCount, mAdapter1Count);
    j["Read1AdapterCounts"] = jR1AdapterCount;
    if(mPaired){
        jsn::json jR2AdapterCount;
        reportAdaptersJsonDetails(jR2AdapterCount, mAdapter2Count);
        j["Read2AdapterCounts"] = jR2AdapterCount;
    }
}

CTML::Node FilterResult::reportAdaptersHtmlSummary(size_t totalBases){
    // adapters
    CTML::Node adapterSection("div.section_div");
    CTML::Node adapterSectionTitle("div.section_title", "Adapters");
    adapterSectionTitle.SetAttribute("onclick", "showOrHide('adapters')");
    CTML::Node adapterSectionTitleLink("a");
    adapterSectionTitleLink.SetAttribute("name", "summary");
    adapterSectionTitle.AppendChild(adapterSectionTitleLink);
    adapterSection.AppendChild(adapterSectionTitle);
    CTML::Node adaptersID("div#adapters");
    // adapters->read1
    CTML::Node adapter1Section("div.subsection_title", "Adapter or bad ligation of read1");
    adapter1Section.SetAttribute("onclick", "showOrHide('read1_adapters')");
    CTML::Node adapter1ID("div#read1_adapters");
    adapter1ID.AppendChild(reportAdaptersHtmlDetails(mAdapter1Count, totalBases));
    adaptersID.AppendChild(adapter1Section);
    adaptersID.AppendChild(adapter1ID);
    // adapters->raed2
    if(mPaired){
        CTML::Node adapter2Section("div.subsection_title", "Adapter or bad ligation of read2");
        adapter2Section.SetAttribute("onclick", "showOrHide('read2_adapters')");
        CTML::Node adapter2ID("div#read2_adapters");
        adapter2ID.AppendChild(reportAdaptersHtmlDetails(mAdapter2Count, totalBases));
        adaptersID.AppendChild(adapter2Section);
        adaptersID.AppendChild(adapter2ID);
    }
    adapterSection.AppendChild(adaptersID);
    return adapterSection;
}

CTML::Node FilterResult::reportPolyXTrimHtml(){
    CTML::Node polyXSection("div.section_div");
    CTML::Node polyXSectionTitle("div.section_title", "PolyX Trimming");
    polyXSectionTitle.SetAttribute("onclick", "showOrHide('polyx')");
    CTML::Node polyXSectionTitleLink("a");
    polyXSectionTitleLink.SetAttribute("name", "summary");
    polyXSectionTitle.AppendChild(polyXSectionTitleLink);
    polyXSection.AppendChild(polyXSectionTitle);
    CTML::Node polyXID("div#polyx");
    CTML::Node polyXTable("table.summary_table");
    polyXTable.AppendChild(htmlutil::make2ColRowNode("TotalPolyXTrimmedReads", std::accumulate(mTrimmedPolyXReads, mTrimmedPolyXReads + 4, 0)));
    polyXTable.AppendChild(htmlutil::make2ColRowNode("TotalPolyXTrimmedBases", std::accumulate(mTrimmedPolyXBases, mTrimmedPolyXBases + 4, 0)));
    std::string nc = "ATCG";
    for(int b = 0; b < 4; ++b){
        polyXTable.AppendChild(htmlutil::make2ColRowNode("ReadsTrimmedByPoly" + std::string(1, nc[b]), mTrimmedPolyXReads[b]));
    }
    for(int b = 0; b < 4; ++b){
        polyXTable.AppendChild(htmlutil::make2ColRowNode("BasesTrimmedByPoly" + std::string(1, nc[b]), mTrimmedPolyXBases[b]));
    }
    polyXID.AppendChild(polyXTable);
    polyXSection.AppendChild(polyXID);
    return polyXSection;
}

void FilterResult::reportPolyXTrimJson(jsn::json& j){
    const char atcg[4] = {'A', 'T', 'C', 'G'};
    j["TotalPolyxTrimmedReads"] = std::accumulate(mTrimmedPolyXReads, mTrimmedPolyXReads + 4, 0);
    jsn::json jPolyReads;
    for(int b = 0; b < 4; ++b){
        jPolyReads[std::string(1, atcg[b])] = mTrimmedPolyXReads[b];
    }
    j["PolyxTrimmedReads"] = jPolyReads;
    j["TotalPolyxTrimmedBases"] = std::accumulate(mTrimmedPolyXBases, mTrimmedPolyXBases + 4, 0);
    jsn::json jPolyBases;
    for(int b = 0; b < 4; ++b){
        jPolyBases[std::string(1, atcg[b])] = mTrimmedPolyXBases[b];
    }
    j["PolyxTrimmedBases"] = jPolyBases;
}
