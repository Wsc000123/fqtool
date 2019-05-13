#include "stats.h"

Stats::Stats(Options* opt, bool isRead2, int bufferMargin){
    mOptions = opt;
    mReads = 0;
    mBases = 0;
    mIsRead2 = isRead2;
    mMinReadLen = 0;
    mMaxReadLen = 0;
    mMinQual = 127;
    mMaxQual = 33;
    mEvaluatedSeqLen = opt->est.seqLen1;
    if(mIsRead2){
        mEvaluatedSeqLen = opt->est.seqLen2;
    }
    mCycles = mEvaluatedSeqLen;
    mBufLen = mEvaluatedSeqLen + bufferMargin;
    mQ20Total = 0;
    mQ30Total = 0;
    mSummarized = false;
    mKmerMax = 0;
    mKmerMin = 0;
    mKmerLen = 0;
    if(opt->kmer.enabled){
        mKmerLen = opt->kmer.kmerLen;
    }
    mKmerBufLen = 1 << (mKmerLen * 2);
    mLengthSum = 0;
    mOverRepSampling = 0;
    if(opt->overRepAna.enabled){
        mOverRepSampling = opt->overRepAna.sampling;
    }
    
    allocateRes();
    
    if(opt->overRepAna.enabled){
        initOverRepSeq();
    }
}

void Stats::allocateRes(){
    for(int i = 0; i < 8; ++i){
        mQ20Bases[i] = 0;
        mQ30Bases[i] = 0;
        mBaseContents[i] = 0;
        mCycleQ20Bases[i] = new size_t[mBufLen];
        std::memset(mCycleQ20Bases[i], 0, sizeof(size_t) * mBufLen);
        mCycleQ30Bases[i] = new size_t[mBufLen];
        std::memset(mCycleQ30Bases[i], 0, sizeof(size_t) * mBufLen);
        mCycleBaseContents[i] = new size_t[mBufLen];
        std::memset(mCycleBaseContents[i], 0, sizeof(size_t) * mBufLen);
        mCycleBaseQuality[i] = new size_t[mBufLen];
        std::memset(mCycleBaseQuality[i], 0, sizeof(size_t) * mBufLen);
    }

    mCycleTotalBase = new size_t[mBufLen];
    std::memset(mCycleTotalBase, 0, sizeof(size_t) * mBufLen);
    mCycleTotalQuality = new size_t[mBufLen];
    std::memset(mCycleTotalQuality, 0, sizeof(size_t) * mBufLen);
    
    if(mKmerLen){
        mKmer = new size_t[mKmerBufLen];
        std::memset(mKmer, 0, sizeof(size_t) * mKmerBufLen);
    }
}

void Stats::extendBuffer(int newBufLen){
    if(newBufLen <= mBufLen){
        return;
    }
    size_t* newBuf = NULL;

    for(int i = 0; i < 8; ++i){
        newBuf = new size_t[newBufLen];
        std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
        std::memcpy(newBuf, mCycleQ20Bases[i], sizeof(size_t) * mBufLen);
        delete mCycleQ20Bases[i];
        mCycleQ20Bases[i] = newBuf;

        newBuf = new size_t[newBufLen];
        std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
        std::memcpy(newBuf, mCycleQ30Bases[i], sizeof(size_t) * mBufLen);
        delete mCycleQ30Bases[i];
        mCycleQ30Bases[i] = newBuf;

        newBuf = new size_t[newBufLen];
        std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
        std::memcpy(newBuf, mCycleBaseContents[i], sizeof(size_t) * mBufLen);
        delete mCycleBaseContents[i];
        mCycleBaseContents[i] = newBuf;

        newBuf = new size_t[newBufLen];
        std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
        std::memcpy(newBuf, mCycleBaseQuality[i], sizeof(size_t) * mBufLen);
        delete mCycleBaseQuality[i];
        mCycleBaseQuality[i] = newBuf;
    }

    newBuf = new size_t[newBufLen];
    std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
    memcpy(newBuf, mCycleTotalBase, sizeof(size_t) * mBufLen);
    delete mCycleTotalBase;
    mCycleTotalBase = newBuf;

    newBuf = new size_t[newBufLen];
    std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
    memcpy(newBuf, mCycleTotalQuality, sizeof(size_t) * mBufLen);
    delete mCycleTotalQuality;
    mCycleTotalQuality = newBuf;

    mBufLen = newBufLen;
}

Stats::~Stats(){
    for(int i = 0; i < 8; ++i){
        delete mCycleQ20Bases[i];
        mCycleQ20Bases[i] = NULL;
        
        delete mCycleQ30Bases[i];
        mCycleQ30Bases[i] = NULL;
        
        delete mCycleBaseContents[i];
        mCycleBaseContents[i] = NULL;
        
        delete mCycleBaseQuality[i];
        mCycleBaseQuality[i] = NULL;
    }
    
    delete mCycleTotalBase;
    delete mCycleTotalQuality;
    
    for(auto& e : mQualityCurves){
        delete e.second;
    }
    
    for(auto& e : mContentCurves){
        delete e.second;
    }
    
    if(mKmerLen){
        delete mKmer;
    }

    deleteOverRepSeqDist();
}

void Stats::summarize(bool forced){
    if(mSummarized && !forced){
        return;
    }

    // cycle and total bases
    bool getMinReadLen = false;
    int c = 0;
    for(c = 0; c < mBufLen; ++c){
        mBases += mCycleTotalBase[c];
        if(!getMinReadLen && c > 1 && mCycleTotalBase[c] < mCycleTotalBase[c-1]){
            mMinReadLen = c;
            getMinReadLen = true;
        }
        if(mCycleTotalBase[c] == 0){
            break;
        }
    }
    mCycles = c;
    mMaxReadLen = c;

    // Q20, Q30, base content
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < mCycles; ++j){
            mQ20Bases[i] += mCycleQ20Bases[i][j];
            mQ30Bases[i] += mCycleQ30Bases[i][j];
            mBaseContents[i] += mCycleBaseContents[i][j];
        }
        mQ20Total += mQ20Bases[i];
        mQ30Total += mQ30Bases[i];
    }

    // quality curve for mean qual
    double* meanQualCurve = new double[mCycles];
    std::memset(meanQualCurve, 0, sizeof(double) * mCycles);
    for(int i = 0; i < mCycles; ++i){
        meanQualCurve[i] = (double)mCycleTotalQuality[i] / (double)mCycleTotalBase[i];
    }
    mQualityCurves["Mean"] = meanQualCurve;

    // quality curves and base contents curves for different nucleotides
    char nucleotides[5] = {'A', 'T', 'C', 'G', 'N'};
    for(int i = 0; i < 5; ++i){
        int b = nucleotides[i] & 0x07;
        double* qualCurve = new double[mCycles];
        std::memset(qualCurve, 0, sizeof(double) * mCycles);
        double* contentCurve = new double[mCycles];
        std::memset(contentCurve, 0, sizeof(double) * mCycles);
        for(int j = 0; j < mCycles; ++j){
            if(mCycleBaseContents[b][j] == 0){
                qualCurve[j] = meanQualCurve[j];
            }else{
                qualCurve[j] = (double)mCycleBaseQuality[b][j] / (double)mCycleBaseContents[b][j];
            }
            contentCurve[j] = (double)mCycleBaseContents[b][j] / (double)mCycleTotalBase[j];
        }
        mQualityCurves[std::string(1, nucleotides[i])] = qualCurve;
        mContentCurves[std::string(1, nucleotides[i])] = contentCurve;
    }

    // GC content curve
    double* gcContentCurve = new double[mCycles];
    std::memset(gcContentCurve, 0, sizeof(double) * mCycles);
    int gIndex = 'G' & 0x07;
    int cIndex = 'C' & 0x07;
    for(int i = 0; i < mCycles; ++i){
        gcContentCurve[i] = (double)(mCycleBaseContents[gIndex][i] + mCycleBaseContents[cIndex][i])/ (double)mCycleTotalBase[i];
    }
    mContentCurves["GC"] = gcContentCurve;

    if(mKmerLen){
        // K-mer statistics
        mKmerMin = mKmer[0];
        mKmerMax = mKmer[0];
        for(int i = 1; i < mKmerBufLen; ++i){
            mKmerMin = std::min(mKmer[i], mKmerMin);
            mKmerMax = std::max(mKmer[i], mKmerMax);
        }
    }

    mSummarized = true;
}

int Stats::getMeanLength(){
    if(mReads == 0){
        return 0;
    }
    return mLengthSum/mReads;
}

void Stats::statRead(Read* r){
    int len = r->length();
    mLengthSum += len;
    if(mBufLen < len){
        extendBuffer(std::max(len + 100, (int)(len * 1.5)));
    }
    const std::string seqStr = r->seq.seqStr;
    const std::string qualStr = r->quality;
    int kmerVal = -1;
    for(int i = 0; i < len; ++i){
        char base = seqStr[i];
        char qual = qualStr[i];
        int bIndex = base & 0x07;
        const char q20 = '5';
        const char q30 = '?';
        mMinQual = std::min(qual - 33, mMinQual);
        mMaxQual = std::max(qual - 33, mMaxQual);
        if(qual > q30){
            ++mCycleQ20Bases[bIndex][i];
            ++mCycleQ30Bases[bIndex][i];
        }else if(qual > q20){
            ++(mCycleQ20Bases[bIndex][i]);
        }

        ++mCycleBaseContents[bIndex][i];
        mCycleBaseQuality[bIndex][i] += (qual - 33);
        ++mCycleTotalBase[i];
        mCycleTotalQuality[i] += (qual - 33);
        
        if(mKmerLen){
            if(i < mKmerLen - 1){
                continue;
            }
            kmerVal = Evaluator::seq2int(seqStr, i - mKmerLen + 1, mKmerLen, kmerVal);
            if(kmerVal >= 0){
                ++mKmer[kmerVal];
            }
        }
    }
    
    if(mOverRepSampling){
        if(mReads % mOverRepSampling == 0){
            std::set<int> steps = {10, 20, 40, 100, std::min(150, mEvaluatedSeqLen - 2)};
            for(auto& step: steps){
                for(int j = 0; j < len - step; ++j){
                    std::string seq = seqStr.substr(j, step);
                    if(mOverReqSeqCount.count(seq) > 0){
                        ++mOverReqSeqCount[seq];
                        for(size_t p = j; p < seq.length() + j && p < (size_t)mEvaluatedSeqLen; ++p){
                            ++mOverRepSeqDist[seq][p];
                        }
                        j += step;
                    }
                }
            }
        }
    }
    ++mReads;
}

int Stats::getMinReadLength(){
    if(!mSummarized){
        summarize();
    }
    return mMinReadLen;
}

int Stats::getMaxReadLength(){
    return getCycles();
}

int Stats::getMinBaseQual(){
    if(!mSummarized){
        summarize();
    }
    return mMinQual;
}

int Stats::getMaxBaseQual(){
    if(!mSummarized){
        summarize();
    }
    return mMaxQual;
}

int Stats::getCycles(){
    if(!mSummarized){
        summarize();
    }
    return mCycles;
}

size_t Stats::getReads(){
    if(!mSummarized){
        summarize();
    }
    return mReads;
}

size_t Stats::getBases(){
    if(!mSummarized){
        summarize();
    }
    return mBases;
}

size_t Stats::getQ20(){
    if(!mSummarized){
        summarize();
    }
    return mQ20Total;
}

size_t Stats::getQ30(){
    if(!mSummarized){
        summarize();
    }
    return mQ30Total;
}

size_t Stats::getGCNumber(){
    if(!mSummarized){
        summarize();
    }
    return mBaseContents['G' & 0x07] + mBaseContents['C' & 0x07];
}

std::ostream& operator<<(std::ostream& os, const Stats& s){
    os << "total reads: " << s.mReads << "\n";
    os << "total bases: " << s.mBases << "\n";
    os << "Q20   bases: " << s.mQ20Total << "(" << (s.mQ20Total*100.0)/s.mBases << "%)\n";
    os << "Q30   bases: " << s.mQ30Total << "(" << (s.mQ30Total*100.0)/s.mBases << "%)\n";
    return os;
}

bool Stats::overRepPassed(const std::string& seq, size_t count){
    int s = mOverRepSampling;
    switch(seq.length()){
        case 10:
            return s * count > 500;
        case 20:
            return s * count > 200;
        case 40:
            return s * count > 100;
        case 100:
            return s * count > 50;
        default:
            return s * count > 20;
    }
}

bool Stats::isLongRead(){
    return mCycles > 300;
}

jsn::json Stats::reportJson(){
    jsn::json summary;
    summary["TotalReads"] = mReads;
    summary["TotalBases"] = mBases;
    summary["Q20Bases"] = mQ20Total;
    summary["Q30Bases"] = mQ30Total;
    summary["TotalCycles"] = mCycles;
    jsn::json qualContent;
    std::string qualNames[5] = {"A", "T", "C", "G", "Mean"};
    for(int i = 0; i < 5; ++i){
        qualContent[qualNames[i]] = std::vector<double>(mQualityCurves[qualNames[i]], mQualityCurves[qualNames[i]] + mCycles);
    }
    summary["QualityCurves"] = qualContent;
    jsn::json baseContent;
    std::string contentNames[6] = {"A", "T", "C", "G", "N", "GC"};
    for(int i = 0; i < 6; ++i){
        baseContent[contentNames[i]] = std::vector<double>(mContentCurves[contentNames[i]], mContentCurves[contentNames[i]] + mCycles);
    }
    summary["ContentCurves"] = baseContent;
    if(mKmerLen){
        jsn::json kmer;
        for(int i = 0; i < mKmerBufLen; ++i){
            std::string seq = Evaluator::int2seq(i, mKmerLen);
            kmer[seq] = std::to_string(mKmer[i]);
        }
        summary["KmerCount"] = kmer;
    }
    if(mOverRepSampling){
        jsn::json ora;
        for(auto& e : mOverReqSeqCount){
            if(!overRepPassed(e.first, e.second)){
                continue;
            }
            ora[e.first] = e.second;
        }
        summary["OverrepresentedSequences"] = ora;
    }
    return summary;
}

std::vector<CTML::Node> Stats::reportHtml(std::string filteringType, std::string readName){
    std::vector<CTML::Node> r;
    r.push_back(reportHtmlQuality(filteringType, readName));
    r.push_back(reportHtmlContents(filteringType, readName));
    if(mKmerLen){
        r.push_back(reportHtmlKmer(filteringType, readName));
    }
    if(mOverRepSampling){
        r.push_back(reportHtmlORA(filteringType, readName));
    }
    return r;
}

CTML::Node Stats::reportHtmlORA(std::string filteringType, std::string readName){
    double dBases = mBases;
    std::string subSection = filteringType + ": " + readName + ": overrepresented sequences";
    std::string divName = util::replace(subSection, " ", "_");
    divName = util::replace(divName, ":", "_");
    std::string title = "";

    CTML::Node oraSection("div.section_div");
    CTML::Node oraSectionTitle("div.subsection_title");
    CTML::Node oraSectionTitleLink("a", subSection);
    oraSectionTitleLink.SetAttribute("title", "click to hide/show");
    oraSectionTitleLink.SetAttribute("onclick", "showOrHide('" + divName + "')");
    oraSectionTitle.AppendChild(oraSectionTitleLink);
    oraSection.AppendChild(oraSectionTitle);
    CTML::Node oraSectionID("div#" + divName);
    oraSectionID.AppendChild(CTML::Node("div.sub_section_tips", "Sampling rate: 1/" + std::to_string(mOverRepSampling)));
    CTML::Node oraSectionTable("table.summary_table");
    CTML::Node oraSectionTableHeader("tr");
    oraSectionTableHeader.SetAttribute("style", "font-weight:bold;");
    oraSectionTableHeader.AppendChild(CTML::Node("td", "overrepresented sequence"));
    oraSectionTableHeader.AppendChild(CTML::Node("td", "count (% of bases)"));
    oraSectionTableHeader.AppendChild(CTML::Node("td", "distribution: cycle 1 ~ cycle " + std::to_string(mEvaluatedSeqLen)));
    oraSectionTable.AppendChild(oraSectionTableHeader);
    
    int found = 0;
    for(auto& e: mOverReqSeqCount){
        std::string seq = e.first;
        size_t count = e.second;
        if(!overRepPassed(seq, count))
            continue;
        found++;
        double percent = (100.0 * count * seq.length() * mOverRepSampling)/dBases;
        CTML::Node oraSectionTableRow("tr");
        CTML::Node col1("td", seq);
        col1.SetAttribute("width", "400").SetAttribute("style", "word-break:break-all;font-size:8px;");
        CTML::Node col2("td", std::to_string(count) + "(" + std::to_string(percent) + "%)");
        col2.SetAttribute("width", "200");
        CTML::Node col3("td");
        col3.SetAttribute("width", "250");
        CTML::Node cavas("canvas#" + divName + "_" + seq);
        cavas.UseClosingTag(false).SetAttribute("width", "240").SetAttribute("height", "20");
        col3.AppendChild(cavas);
        oraSectionTableRow.AppendChild(col1).AppendChild(col2).AppendChild(col3);
        oraSectionTable.AppendChild(oraSectionTableRow);
    }
    if(found == 0){
        CTML::Node oraSectionTableRowNt("tr");
        CTML::Node col("td", "not found");
        col.SetAttribute("style", "text-align:center").SetAttribute("colspan", "3");
        oraSectionTableRowNt.AppendChild(col);
        oraSectionTable.AppendChild(oraSectionTableRowNt);
    }
    oraSectionID.AppendChild(oraSectionTable);
    CTML::Node oraSectionJS("script");
    oraSectionJS.SetAttribute("language", "javascript");
    // output the JS
    std::stringstream ss;
    ss << "var seqlen = " << mEvaluatedSeqLen << ";" << std::endl;
    ss << "var orp_dist = {" << std::endl;
    bool first = true;
    for(auto& e: mOverReqSeqCount){
        std::string seq = e.first;
        size_t count = e.second;
        if(!overRepPassed(seq, count))
            continue;
        if(!first) {
            ss << "," << std::endl;
        } else
            first = false;
        ss << "\t\"" << divName << "_" << seq << "\":[";
        for(int i=0; i< mEvaluatedSeqLen; ++i){
            if(i !=0 )
                ss << ",";
            ss << mOverRepSeqDist[seq][i];
        }
        ss << "]";
    }
    ss << "\n};" << std::endl;
    ss << "for (seq in orp_dist) {"<< std::endl;
    ss << "    var cvs = document.getElementById(seq);"<< std::endl;
    ss << "    var ctx = cvs.getContext('2d'); "<< std::endl;
    ss << "    var data = orp_dist[seq];"<< std::endl;
    ss << "    var w = 240;"<< std::endl;
    ss << "    var h = 20;"<< std::endl;
    ss << "    ctx.fillStyle='#cccccc';"<< std::endl;
    ss << "    ctx.fillRect(0, 0, w, h);"<< std::endl;
    ss << "    ctx.fillStyle='#0000FF';"<< std::endl;
    ss << "    var maxVal = 0;"<< std::endl;
    ss << "    for(d=0; d<seqlen; d++) {"<< std::endl;
    ss << "        if(data[d]>maxVal) maxVal = data[d];"<< std::endl;
    ss << "    }"<< std::endl;
    ss << "    var step = (seqlen-1) /  (w-1);"<< std::endl;
    ss << "    for(x=0; x<w; x++){"<< std::endl;
    ss << "        var target = step * x;"<< std::endl;
    ss << "        var val = data[Math.floor(target)];"<< std::endl;
    ss << "        var y = Math.floor((val / maxVal) * h);"<< std::endl;
    ss << "        ctx.fillRect(x,h-1, 1, -y);"<< std::endl;
    ss << "    }"<< std::endl;
    ss << "}"<< std::endl;
    oraSectionJS.AppendText(ss.str());
    oraSectionID.AppendChild(oraSectionJS);
    oraSection.AppendChild(oraSectionID);
    return oraSection;
}

CTML::Node Stats::reportHtmlKmer(std::string filteringType, std::string readName) {
    // KMER
    std::string subsection = filteringType + ": " + readName + ": KMER counting";
    std::string divName = util::replace(subsection, " ", "_");
    divName = util::replace(divName, ":", "_");
    
    CTML::Node kmerSection("div.section_div");
    CTML::Node kmerSectionTitle("div.subsection_title");
    CTML::Node kmerSectionTitleLink("a", subsection);
    kmerSectionTitleLink.SetAttribute("title", "click to hide/show' onclick=showOrHide('" + divName + "')");
    kmerSectionTitle.AppendChild(kmerSectionTitleLink);
    kmerSection.AppendChild(kmerSectionTitle);
    CTML::Node kmerSectionID("div#" + divName);
    kmerSectionID.AppendChild(CTML::Node("div.sub_section_tips", "Darker background means larger counts. The count will be shown on mouse over"));
    CTML::Node kmerSectionTable("table.kmer_table");
    kmerSectionTable.SetAttribute("style", "width:680px;");
    CTML::Node kmerSectionTableHeader("tr");
    // the heading row
    kmerSectionTableHeader.AppendChild(CTML::Node("td"));
    for(int h = 0; h < (1 << mKmerLen); ++h){
        CTML::Node row("td", std::to_string(h + 1));
        row.SetAttribute("style", "color:#333333");
        kmerSectionTableHeader.AppendChild(row);
    }
    kmerSectionTable.AppendChild(kmerSectionTableHeader);
    // content
    size_t n = 0;
    for(size_t i = 0; i < (1 << mKmerLen); ++i){
        CTML::Node kmerSectionTableRow("tr");
        CTML::Node firCol("td", std::to_string(i + 1));
        firCol.SetAttribute("style", "color:#333333");
        kmerSectionTableRow.AppendChild(firCol);
        for(int j = 0; j < (1 << mKmerLen); ++j){
            kmerSectionTableRow.AppendChild(makeKmerTD(n++));
        }
        kmerSectionTable.AppendChild(kmerSectionTableRow);
    }
    kmerSectionID.AppendChild(kmerSectionTable);
    kmerSection.AppendChild(kmerSectionID);
    return kmerSection;
}

CTML::Node Stats::makeKmerTD(size_t n){
    std::string seq = Evaluator::int2seq(n, mKmerLen);
    double meanBases = (double)(mBases + 1) / mKmerBufLen;
    double prop = mKmer[n] / meanBases;
    double frac = 0.5;
    if(prop > 2.0){
        frac = (prop-2.0)/20.0 + 0.5;
    }
    else if(prop< 0.5){
        frac = prop;
    }
    frac = std::max(0.01, std::min(1.0, frac));
    int r = (1.0-frac) * 255;
    int g = r;
    int b = r;
    std::stringstream ss;
    if(r<16){
        ss << "0";
    }
    ss << std::hex <<r;
    if(g<16){
        ss << "0";
    }
    ss << std::hex <<g;
    if(b<16){
        ss << "0";
    }
    ss << std::hex << b;
    CTML::Node row("td", seq);
    row.SetAttribute("style", "background:#" + ss.str());
    row.SetAttribute("title", seq + ": " + std::to_string(mKmer[n]) + "\n" + std::to_string(prop) + " times as mean value");
    return row;
}

CTML::Node Stats::reportHtmlQuality(std::string filteringType, std::string readName) {
    // quality
    std::string subsection = filteringType + ": " + readName + ": quality";
    std::string divName = util::replace(subsection, " ", "_");
    divName = util::replace(divName, ":", "_");
    std::string title = "";
    std::string alphabets[5] = {"A", "T", "C", "G", "Mean"};
    std::string colors[5] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(20,20,20,1.0)"};
    std::string json_str = "var data=[";
    int *x = new int[mCycles];
    int total = 0;
    if(!isLongRead()) {
        for(int i = 0; i < mCycles; ++i){
            x[total] = i + 1;
            total++;
        }
    } else {
        const int fullSampling = 40;
        for(int i=0; i < fullSampling && i < mCycles; i++){
            x[total] = i + 1;
            total++;
        }
        // down sampling if it's too long
        if(mCycles > fullSampling) {
            double pos = fullSampling;
            while(true){
                pos *= 1.05;
                if(pos >= mCycles)
                    break;
                x[total] = (int)pos;
                total++;
            }
            // make sure lsat one is contained
            if(x[total-1] != mCycles){
                x[total] = mCycles;
                total++;
            }
        }
    }
    // four bases
    for (int b = 0; b<5; b++) {
        std::string base = alphabets[b];
        json_str += "{";
        json_str += "x:[" + list2string(x, total) + "],";
        json_str += "y:[" + list2string(mQualityCurves[base], total) + "],";
        json_str += "name: '" + base + "',";
        json_str += "mode:'lines',";
        json_str += "line:{color:'" + colors[b] + "', width:1}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'position'";
    json_str += ", tickmode: 'auto', nticks: '";
    json_str += std::to_string(mCycles/5) + "'";
    // use log plot if it's too long
    if(isLongRead()) {
        json_str += ",type:'log'";
    }
    json_str += "},";
    json_str += "yaxis:{title:'quality', tickmode: 'auto', nticks: '20'";
    json_str += "}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";
    delete[] x;
    
    CTML::Node qualSection("div.section_div");
    CTML::Node qualSectionTitle("div.subsection_title");
    CTML::Node qualSectionTitleLink("a", subsection);
    qualSectionTitleLink.SetAttribute("title", "click to hide/show' onclick=showOrHide('" + divName + "')");
    qualSectionTitle.AppendChild(qualSectionTitleLink);
    qualSection.AppendChild(qualSectionTitle);
    CTML::Node qualSectionID("div#" + divName);
    qualSectionID.AppendChild(CTML::Node("div.sub_section_tips", "Value of each position will be shown on mouse over"));
    qualSectionID.AppendChild(CTML::Node("div.figure#plot_" + divName));
    qualSection.AppendChild(qualSectionID);
    CTML::Node qualJSC("script");
    qualJSC.SetAttribute("type", "text/javascript");
    qualJSC.AppendText(json_str);
    qualSection.AppendChild(qualJSC);
    return qualSection;
}

CTML::Node Stats::reportHtmlContents(std::string filteringType, std::string readName) {
    // content
    std::string subsection = filteringType + ": " + readName + ": base contents";
    std::string divName = util::replace(subsection, " ", "_");
    divName = util::replace(divName, ":", "_");
    std::string title = "";

    CTML::Node contentSection("div.section_div");
    CTML::Node contentSectionTitle("div.subsection_title");
    CTML::Node contentSectionTitleClick("a", subsection);
    contentSectionTitleClick.SetAttribute("title", "click to hide/show' onclick=showOrHide('" + divName + "')");
    contentSectionTitle.AppendChild(contentSectionTitleClick);
    contentSection.AppendChild(contentSectionTitle);
    CTML::Node contentSectionID("div#" + divName);
    contentSectionID.AppendChild(CTML::Node("div.sub_section_tips", "Value of each position will be shown on mouse over"));
    contentSectionID.AppendChild(CTML::Node("div.figure#plot_" + divName));
    contentSection.AppendChild(contentSectionID);
    CTML::Node contentSectionJS("script");
    contentSectionJS.SetAttribute("type", "text/javascript");

    std::string alphabets[6] = {"A", "T", "C", "G", "N", "GC"};
    std::string colors[6] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(255, 0, 0, 1.0)", "rgba(20,20,20,1.0)"};
    std::string json_str = "var data=[";
    int *x = new int[mCycles];
    int total = 0;
    if(!isLongRead()) {
        for(int i=0; i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
    } else {
        const int fullSampling = 40;
        for(int i=0; i<fullSampling && i<mCycles; i++){
            x[total] = i+1;
             total++;
        }
        // down sampling if it's too long
        if(mCycles > fullSampling) {
            double pos = fullSampling;
            while(true){
                pos *= 1.05;
                if(pos >= mCycles)
                    break;
                x[total] = (int)pos;
                total++;
            }
            // make sure lsat one is contained
            if(x[total-1] != mCycles){
                x[total] = mCycles;
                total++; 
             }
        } 
    }
    // four bases
    for (int b = 0; b<6; b++) {
        std::string base = alphabets[b];
        long count = 0;
        if(base.size()==1) {
            int b = base[0] & 0x07;
            count = mBaseContents[b];
        } else {
            count = mBaseContents['G' & 0x07] + mBaseContents['C' & 0x07] ;
        }
        std::string percentage = std::to_string((double)count * 100.0 / mBases);
        if(percentage.length()>5)
            percentage = percentage.substr(0,5);
        std::string name = base + "(" + percentage + "%)";

        json_str += "{";
        json_str += "x:[" + list2string(x, total) + "],";
        json_str += "y:[" + list2string(mContentCurves[base], total) + "],";
        json_str += "name: '" + name + "',";
        json_str += "mode:'lines',";
        json_str += "line:{color:'" + colors[b] + "', width:1}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'position'";
    json_str += ", tickmode: 'auto', nticks: '";
    json_str += std::to_string(mCycles/5) + "'";
    // use log plot if it's too long
    if(isLongRead()) {
        json_str += ",type:'log'";
    }
    json_str += "}, yaxis:{title:'base content ratios'";
    json_str += ", tickmode: 'auto', nticks: '20', range: ['0.0', '1.0']";
    json_str += "}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";
    delete[] x;
    contentSectionJS.AppendText(json_str);
    contentSection.AppendChild(contentSectionJS);
    return contentSection;
}

Stats* Stats::merge(std::vector<Stats*>& list){
    if(list.size() == 0){
        return NULL;
    }
    
    int c = 0;
    for(size_t i = 0; i < list.size(); ++i){
        list[i]->summarize();
        c = std::max(c, list[i]->getCycles());
    }

    Stats* s = new Stats(list[0]->mOptions, list[0]->mIsRead2);
    for(size_t i = 0; i < list.size(); ++i){
        int curCycles = list[i]->getCycles();
        s->mReads += list[i]->mReads;
        s->mLengthSum += list[i]->mLengthSum;
        for(int j = 0; j < 8; ++j){
            for(int k = 0; k < c && k < curCycles; ++k){
                s->mCycleQ30Bases[j][k] += list[i]->mCycleQ30Bases[j][k];
                s->mCycleQ20Bases[j][k] += list[i]->mCycleQ20Bases[j][k];
                s->mCycleBaseContents[j][k] += list[i]->mCycleBaseContents[j][k];
                s->mCycleBaseQuality[j][k] += list[i]->mCycleBaseQuality[j][k];
            }
        }

        for(int j = 0; j < c && j < curCycles; ++j){
            s->mCycleTotalBase[j] += list[i]->mCycleTotalBase[j];
            s->mCycleTotalQuality[j] += list[i]->mCycleTotalQuality[j];
        }

        if(s->mKmerLen){
            for(int j = 0; j < s->mKmerBufLen; ++j){
                s->mKmer[j] += list[i]->mKmer[j];
            }
        }

        if(s->mOverRepSampling){
            for(auto& e: s->mOverReqSeqCount){
                s->mOverReqSeqCount[e.first] += list[i]->mOverReqSeqCount[e.first];
                for(int j = 0; j < s->mEvaluatedSeqLen; ++j){
                    s->mOverRepSeqDist[e.first][j] += list[i]->mOverRepSeqDist[e.first][j];
                }
            }
        }
    }

    s->summarize();
    return s;
}

void Stats::initOverRepSeq(){
    std::map<std::string, size_t> *mapORS;
    if(mIsRead2){
        mapORS = &(mOptions->overRepAna.overRepSeqCountR2);
    }else{
        mapORS = &(mOptions->overRepAna.overRepSeqCountR1);
    }
    for(auto& e: *mapORS){
        mOverReqSeqCount[e.first] = 0;
        mOverRepSeqDist[e.first] = new size_t[mEvaluatedSeqLen];
        std::memset(mOverRepSeqDist[e.first], 0, sizeof(size_t) * mEvaluatedSeqLen);
    }
}

void Stats::deleteOverRepSeqDist(){
    for(auto& e: mOverReqSeqCount){
        delete mOverRepSeqDist[e.first];
        mOverRepSeqDist[e.first] = NULL;
    }
}
