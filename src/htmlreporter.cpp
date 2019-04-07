#include "htmlreporter.h"

HtmlReporter::HtmlReporter(Options* opt){
    mOptions = opt;
    mDupHist = NULL;
    mDupRate = 0;
}

HtmlReporter::~HtmlReporter(){
}

void HtmlReporter::setDupHist(size_t* dupHist, double* dupMeanGC, double dupRate){
    mDupHist = dupHist;
    mDupMeanGC = dupMeanGC;
    mDupRate = dupRate;
}

void HtmlReporter::setInsertHist(long* insertHist, int insertSizePeak){
    mInsertHist = insertHist;
    mInsertSizePeak = insertSizePeak;
}

void HtmlReporter::report(FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2) {
    std::ofstream ofs(mOptions->htmlFile);
    htmlutil::printHeader(ofs, "Fastq Preprocess Report");
    printSummary(ofs, result, preStats1, postStats1, preStats2, postStats2);
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('before_filtering')><a name='summary'>Before filtering</a></div>\n";
    ofs << "<div id='before_filtering'>\n";
    if(preStats1){
        preStats1->reportHtml(ofs, "Before filtering", "read1");
    }
    if(preStats2){
        preStats2->reportHtml(ofs, "Before filtering", "read2");
    }
    ofs << "</div>\n";
    ofs << "</div>\n";
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('after_filtering')><a name='summary'>After filtering</a></div>\n";
    ofs << "<div id='after_filtering'>\n";
    if(postStats1){
        postStats1->reportHtml(ofs, "After filtering", "read1");
    }
    if(postStats2){
        postStats2->reportHtml(ofs, "After filtering", "read2");
    }
    ofs << "</div>\n";
    ofs << "</div>\n";
    htmlutil::printFooter(ofs, "Finish Report");
}

void HtmlReporter::printSummary(std::ofstream& ofs, FilterResult* fresult, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2){
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

    std::string sequencingInfo  = mOptions->isPaired()?"paired end":"single end";
    if(mOptions->isPaired()){
        sequencingInfo += " (" + std::to_string(preStats1->getCycles()) + " cycles + " + std::to_string(preStats2->getCycles()) + " cycles)";
    }else{
        sequencingInfo += " (" + std::to_string(preStats1->getCycles()) + " cycles)";
    }
    ofs << std::endl;
    ofs << "<h1 style='text-align:left;'><style='color:#663355;text-decoration:none;'>" + mOptions->reportTitle + "</a>" << std::endl;
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('summary')><a name='summary'>Summary</a></div>\n";
    ofs << "<div id='summary'>\n";
    ofs << "<div class='subsection_title' onclick=showOrHide('general')>General</div>\n";
    ofs << "<div id='general'>\n";
    ofs << "<table class='summary_table'>\n";
    htmlutil::outputTableRow(ofs, "sequencing:", sequencingInfo);

    ofs << "<div class='subsection_title' onclick=showOrHide('before_filtering_summary')>Before filtering</div>\n";
    ofs << "<div id='before_filtering_summary'>\n";
    ofs << "<table class='summary_table'>\n";

    htmlutil::outputTableRow(ofs, "total_reads", preTotalReads);
    htmlutil::outputTableRow(ofs, "total_bases", preTotalBases);
    htmlutil::outputTableRow(ofs, "Q20_bases", preQ20Bases);
    htmlutil::outputTableRow(ofs, "Q30_bases", preQ30Bases);
    htmlutil::outputTableRow(ofs, "Q20_rate", preQ20Rate);
    htmlutil::outputTableRow(ofs, "Q30_rate", preQ30Rate);
    htmlutil::outputTableRow(ofs, "read1_mean_length", preRead1Length);
    if(mOptions->isPaired()){
         htmlutil::outputTableRow(ofs, "read2_mean_length", preRead2Length);
    }
    htmlutil::outputTableRow(ofs, "gc_content", preGCRate);
    ofs << "</table>\n";
    ofs << "</div>\n";
    
    ofs << "<div class='subsection_title' onclick=showOrHide('after_filtering_summary')>After filtering</div>\n";
    ofs << "<div id='after_filtering_summary'>\n";
    ofs << "<table class='summary_table'>\n";
    htmlutil::outputTableRow(ofs, "total_reads", postTotalReads);
    htmlutil::outputTableRow(ofs, "total_bases", postTotalBases);
    htmlutil::outputTableRow(ofs, "Q20_bases", postQ20Bases);
    htmlutil::outputTableRow(ofs, "Q30_bases", postQ30Bases);
    htmlutil::outputTableRow(ofs, "Q20_rate", postQ20Rate);
    htmlutil::outputTableRow(ofs, "Q30_rate", postQ30Rate);
    htmlutil::outputTableRow(ofs, "read1_mean_length", postRead1Length);
    if(mOptions->isPaired()){
         htmlutil::outputTableRow(ofs, "read2_mean_length", postRead2Length);
    }
    htmlutil::outputTableRow(ofs, "gc_content", postGCRate);
    ofs << "</table>\n";
    ofs << "</div>\n";

    if(fresult){
        ofs << "<div class='subsection_title' onclick=showOrHide('filtering_result')>Filtering result</div>\n";
        ofs << "<div id='filtering_result'>\n";
        fresult->reportHtmlBasic(ofs, preTotalReads, preTotalBases);
        ofs << "</div>\n";
    }
    ofs << "</table>\n";
    ofs << "</div>\n";
}
