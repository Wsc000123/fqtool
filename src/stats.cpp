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
    mQualityCurves["mean"] = meanQualCurve;

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

void Stats::reportJson(std::ofstream& ofs, std::string padding){
    // braces
    ofs << "{\n";
    // basic qc info 
    ofs << padding << "\t\"total_reads\": " << mReads << ",\n";
    ofs << padding << "\t\"total_bases\": " << mBases << ",\n";
    ofs << padding << "\t\"q20_bases\": " << mQ20Total << ",\n";
    ofs << padding << "\t\"q30_bases\": " << mQ30Total << ",\n";
    ofs << padding << "\t\"total_cycles\": " << mCycles << ",\n";
    // quality curves
    std::string qualNames[5] = {"A", "T", "C", "G", "mean"};
    ofs << padding << "\t\"quality_curves\": {\n";
    for(int i = 0; i < 5; ++i){
        std::string name = qualNames[i];
        double* curve = mQualityCurves[name];
        ofs << padding << "\t\t\"" << name << "\":[";
        for(int c = 0; c < mCycles; ++c){
            ofs << curve[c];
            if(c != mCycles - 1){
                ofs << ",";
            }
        }
        ofs << "]";
        if(i != 4){
            ofs << ",";
        }
        ofs << "\n";
    }
    ofs << padding << "\t},\n";
    //content curves
    std::string contentNames[6] = {"A", "T", "C", "G", "N", "GC"};
    ofs << padding << "\t\"content_curves\":{\n";
    for(int i = 0; i < 6; ++i){
        std::string name = contentNames[i];
        double* curve = mContentCurves[name];
        ofs << padding << "\t\t\"" << name << "\":[";
        for(int c = 0; c < mCycles; ++c){
            ofs << curve[c];
            if(c != mCycles - 1){
                ofs << ",";
            }
        }
        ofs << "]";
        if(i != 5){
            ofs << ",";
        }
        ofs << "\n";
    }
    ofs << padding << "\t}";
    //KMER counting(optional)
    std::string maxKmerInt = std::to_string(mKmerMax);
    int maxWidth = maxKmerInt.length();
    if(mKmerLen){
        ofs << ",\n";
        ofs << padding << "\t" << "\"kmer_count\": {\n";
        for(int i = 0; i < mKmerBufLen; ++i){
            std::string seq = Evaluator::int2seq(i, mKmerLen);
            if(i % (1 << mKmerLen)){
                ofs << "\"" << seq << "\":" << std::left << std::setw(maxWidth + 1) << std::to_string(mKmer[i]) +",";
            }else{
                ofs << "\n" << padding << "\t\t\"" << seq << "\":" << std::left << std::setw(maxWidth + 1) << std::to_string(mKmer[i]) + ",";
            }
        }
        ofs << padding << "\t}";
    }
    // over represented seqs(optional)
    if(mOverRepSampling){
        ofs << ",\n";
        ofs << padding << "\t\"overrepresented_sequences\": {\n";
        bool first = true;
        for(auto& e : mOverReqSeqCount){
            if(!overRepPassed(e.first, e.second)){
                continue;
            }
            if(!first){
                ofs << ",\n";
            }else{
                first = false;
            }
            ofs << padding << "\t\t\"" << e.first << "\":" << e.second;
        }
        ofs << padding << "\t}\n";
    }
    ofs << padding << "\n}\n";
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

void Stats::reportHtml(std::ofstream& ofs, std::string filteringType, std::string readName){
    reportHtmlQuality(ofs, filteringType, readName);
    reportHtmlContents(ofs, filteringType, readName);
    if(mKmerLen){
        reportHtmlKmer(ofs, filteringType, readName);
    }
    if(mOverRepSampling){
        reportHtmlORA(ofs, filteringType, readName);
    }
}

void Stats::reportHtmlORA(std::ofstream& ofs, std::string filteringType, std::string readName){
    double dBases = mBases;
    std::string subSection = filteringType + ": " + readName + ": overrepresented sequences";
    std::string divName = util::replace(subSection, " ", "_");
    divName = util::replace(divName, ":", "_");
    std::string title = "";
    
    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('";
    ofs << divName << "')>" << subSection << "</a></div>\n";
    ofs << "<div id='" << divName << "'>\n";
    ofs << "<div class='sub_section_tips'>Sampling rate: 1 / " << mOverRepSampling << "</div>\n";
    ofs << "<table class='summary_table'>\n";
    ofs << "<tr style='font-weight:bold;'><td>overrepresented sequence</td><td>count (% of bases)</td><td>distribution: cycle 1 ~ cycle " << mEvaluatedSeqLen << "</td></tr>"<<std::endl;
    int found = 0;
    for(auto& e: mOverReqSeqCount){
        std::string seq = e.first;
        size_t count = e.second;
        if(!overRepPassed(seq, count))
            continue;
        found++;
        double percent = (100.0 * count * seq.length() * mOverRepSampling)/dBases;
        ofs << "<tr>";
        ofs << "<td width='400' style='word-break:break-all;font-size:8px;'>" << seq << "</td>";
        ofs << "<td width='200'>" << count << " (" << std::to_string(percent) <<"%)</td>";
        ofs << "<td width='250'><canvas id='" << divName << "_" << seq << "' width='240' height='20'></td>";
        ofs << "</tr>" << std::endl;
    }
    if(found == 0)
        ofs << "<tr><td style='text-align:center' colspan='3'>not found</td></tr>" << std::endl;
    ofs << "</table>\n";
    ofs << "</div>\n";

    // output the JS
    ofs << "<script language='javascript'>" << std::endl;
    ofs << "var seqlen = " << mEvaluatedSeqLen << ";" << std::endl;
    ofs << "var orp_dist = {" << std::endl;
    bool first = true;
    for(auto& e: mOverReqSeqCount){
        std::string seq = e.first;
        size_t count = e.second;
        if(!overRepPassed(seq, count))
            continue;
        if(!first) {
            ofs << "," << std::endl;
        } else
            first = false;
        ofs << "\t\"" << divName << "_" << seq << "\":[";
        for(int i=0; i< mEvaluatedSeqLen; ++i){
            if(i !=0 )
                ofs << ",";
            ofs << mOverRepSeqDist[seq][i];
        }
        ofs << "]";
    }
        ofs << "\n};" << std::endl;

    ofs << "for (seq in orp_dist) {"<< std::endl;
    ofs << "    var cvs = document.getElementById(seq);"<< std::endl;
    ofs << "    var ctx = cvs.getContext('2d'); "<< std::endl;
    ofs << "    var data = orp_dist[seq];"<< std::endl;
    ofs << "    var w = 240;"<< std::endl;
    ofs << "    var h = 20;"<< std::endl;
    ofs << "    ctx.fillStyle='#cccccc';"<< std::endl;
    ofs << "    ctx.fillRect(0, 0, w, h);"<< std::endl;
    ofs << "    ctx.fillStyle='#0000FF';"<< std::endl;
    ofs << "    var maxVal = 0;"<< std::endl;
    ofs << "    for(d=0; d<seqlen; d++) {"<< std::endl;
    ofs << "        if(data[d]>maxVal) maxVal = data[d];"<< std::endl;
    ofs << "    }"<< std::endl;
    ofs << "    var step = (seqlen-1) /  (w-1);"<< std::endl;
    ofs << "    for(x=0; x<w; x++){"<< std::endl;
    ofs << "        var target = step * x;"<< std::endl;
    ofs << "        var val = data[Math.floor(target)];"<< std::endl;
    ofs << "        var y = Math.floor((val / maxVal) * h);"<< std::endl;
    ofs << "        ctx.fillRect(x,h-1, 1, -y);"<< std::endl;
    ofs << "    }"<< std::endl;
    ofs << "}"<< std::endl;
    ofs << "</script>"<< std::endl;
}

void Stats::reportHtmlKmer(std::ofstream& ofs, std::string filteringType, std::string readName) {
    // KMER
    std::string subsection = filteringType + ": " + readName + ": KMER counting";
    std::string divName = util::replace(subsection, " ", "_");
    divName = util::replace(divName, ":", "_");
    std::string title = "";

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div  id='" << divName << "'>\n";
    ofs << "<div class='sub_section_tips'>Darker background means larger counts. The count will be shown on mouse over.</div>\n";
    ofs << "<table class='kmer_table' style='width:680px;'>\n";
    ofs << "<tr>";
    ofs << "<td></td>";
    // the heading row
    for(int h = 0; h < (1 << mKmerLen); ++h){
        ofs << "<td style='color:#333333'>" << std::to_string(h + 1) << "</td>";
    }
    ofs << "</tr>\n";
    // content
    size_t n = 0;
    for(size_t i = 0; i < (1 << mKmerLen); ++i){
        ofs << "<tr>";
        ofs << "<td style='color:#333333'>" << std::to_string(i + 1) << "</td>";
        for(int j = 0; j < (1 << mKmerLen); ++j){
            ofs << makeKmerTD(n++);
        }
        ofs << "</tr>\n";
    }
    ofs << "</table>\n";
    ofs << "</div>\n";
}

std::string Stats::makeKmerTD(size_t n){
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
    ss << "<td style='background:#";
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
    ss << std::dec << "' title='"<< seq << ": " << mKmer[n] << "\n" << prop << " times as mean value'>";
    ss << seq << "</td>";
    return ss.str();
}

void Stats::reportHtmlQuality(std::ofstream& ofs, std::string filteringType, std::string readName) {
    // quality
    std::string subsection = filteringType + ": " + readName + ": quality";
    std::string divName = util::replace(subsection, " ", "_");
    divName = util::replace(divName, ":", "_");
    std::string title = "";
    ofs << "<script type=\"text/javascript\" src=\"./plotly.js\"></script>\n";
    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";

    std::string alphabets[5] = {"A", "T", "C", "G", "mean"};
    std::string colors[5] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << std::endl;
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

    ofs << json_str;
    ofs << "</script>" << std::endl;

    delete[] x;
}

void Stats::reportHtmlContents(std::ofstream& ofs, std::string filteringType, std::string readName) {

    // content
    std::string subsection = filteringType + ": " + readName + ": base contents";
    std::string divName = util::replace(subsection, " ", "_");
    divName = util::replace(divName, ":", "_");
    std::string title = "";

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";

    std::string alphabets[6] = {"A", "T", "C", "G", "N", "GC"};
    std::string colors[6] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(255, 0, 0, 1.0)", "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << std::endl;
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

    ofs << json_str;
    ofs << "</script>" << std::endl;

    delete[] x;
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
