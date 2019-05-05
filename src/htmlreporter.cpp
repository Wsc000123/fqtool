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
    printHeader(ofs, "Fastq Preprocess Report");
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
    ofs << "<table class='summary_table'>" << std::endl;;

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
    ofs << "</div>" << std::endl;

    if(fresult){
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='subsection_title' onclick=showOrHide('filtering_result')>Filtering result</div>\n";
        ofs << "<div id='filtering_result'>\n";
        fresult->reportHtmlBasic(ofs, preTotalReads, preTotalBases);
        ofs << "</div>\n";
    }

    ofs << "</div>\n";
    ofs << "</div>\n";

    if(fresult && mOptions->adapter.enableTriming){
         ofs << "<div class='section_div'>\n";
         ofs << "<div class='section_title' onclick=showOrHide('adapters')><a name='summary'>Adapters</a></div>\n";
         ofs << "<div id='adapters'>\n";
         fresult->reportAdaptersHtmlSummary(ofs, preTotalBases);
         ofs << "</div>\n";
         ofs << "</div>\n";
    }

    if(mOptions->duplicate.enabled){
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('duplication')><a name='summary'>Duplication</a></div>\n";
        ofs << "<div id='duplication'>\n";
        reportDuplication(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";
    }
}

void HtmlReporter::reportDuplication(std::ofstream& ofs){
    ofs << "<div id='duplication_figure'>\n";
    ofs << "<div class='figure' id='plot_duplication' style='height:400px;'></div>\n";
    ofs << "</div>\n";

    ofs << "\n<script type=\"text/javascript\">" << std::endl;
    std::string json_str = "var data=[";

    int total = mOptions->duplicate.histSize - 2;
    long *x = new long[total];
    double allCount = 0;
    for(int i=0; i<total; i++) {
        x[i] = i+1;
        allCount += mDupHist[i+1];
    }
    double* percents = new double[total];
    memset(percents, 0, sizeof(double)*total);
    if(allCount > 0) {
        for(int i=0; i<total; i++) {
            percents[i] = (double)mDupHist[i+1] * 100.0 / (double)allCount;
        }
    }
    int maxGC = total;
    double* gc = new double[total];
    for(int i=0; i<total; i++) {
        gc[i] = (double)mDupMeanGC[i+1] * 100.0;
        // GC ratio will be not accurate if no enough reads to average
        if(percents[i] <= 0.05 && maxGC == total)
            maxGC = i;
    }
    
    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, total) + "],";
    json_str += "y:[" + Stats::list2string(percents, total) + "],";
    json_str += "name: 'Read percent (%)  ',";
    json_str += "type:'bar',";
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "},";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, maxGC) + "],";
    json_str += "y:[" + Stats::list2string(gc, maxGC) + "],";
    json_str += "name: 'Mean GC ratio (%)  ',";
    json_str += "mode:'lines',";
    json_str += "line:{color:'rgba(255,0,128,1.0)', width:2}\n";
    json_str += "}";
    json_str += "];\n";
    
    json_str += "var layout={title:'duplication rate (" + std::to_string(mDupRate*100.0) + "%)', xaxis:{title:'duplication level'}, yaxis:{title:'Read percent (%) & GC ratio'}};\n";
    json_str += "Plotly.newPlot('plot_duplication', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << std::endl;
    delete[] x;
    delete[] percents;
    delete[] gc;
}

 void HtmlReporter::printCSS(std::ofstream& ofs){
    ofs << "<style type=\"text/css\">" << std::endl;
    ofs << "td {border:1px solid #dddddd;padding:5px;font-size:12px;}" << std::endl;
    ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:800px}" << std::endl;
    ofs << ".col1 {width:240px; font-weight:bold;}" << std::endl;
    ofs << ".adapter_col {width:500px; font-size:10px;}" << std::endl;
    ofs << "img {padding:30px;}" << std::endl;
    ofs << "#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}" << std::endl;
    ofs << "#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helv  etica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}" << std::endl;
    ofs << "a:visited {color: #999999}" << std::endl;
    ofs << ".alignleft {text-align:left;}" << std::endl;
    ofs << ".alignright {text-align:right;}" << std::endl;
    ofs << ".figure {width:800px;height:600px;}" << std::endl;
    ofs << ".header {color:#ffffff;padding:1px;height:20px;background:#000000;}" << std::endl;
    ofs << ".section_title {color:#ffffff;font-size:20px;padding:5px;text-align:left;background:#663355; margin-top:10px;}" << std::endl;
    ofs << ".subsection_title {font-size:16px;padding:5px;margin-top:10px;text-align:left;color:#663355}" << std::endl;
    ofs << "#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << std::endl;
    ofs << ".menu_item {text-align:left;padding-top:5px;font-size:18px;}" << std::endl;
    ofs << ".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}" << std::endl;
    ofs << "#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}" << std::endl;
    ofs << "#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#663355;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << std::endl;
    ofs << ".kmer_table {text-align:center;font-size:8px;padding:2px;}" << std::endl;
    ofs << ".kmer_table td{text-align:center;font-size:8px;padding:0px;color:#ffffff}" << std::endl;
    ofs << ".sub_section_tips {color:#999999;font-size:10px;padding-left:5px;padding-bottom:3px;}" << std::endl;
    ofs << "</style>" << std::endl;
}

void HtmlReporter::printHeader(std::ofstream& ofs, const std::string& title){
    ofs << "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    ofs << "<title>" << title << " </title>";
    printJS(ofs);        
    printCSS(ofs);
    ofs << "</head>";
    ofs << "<body><div id='container'>";
}

void HtmlReporter::printJS(std::ofstream& ofs){
    ofs << "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>" << std::endl;
    ofs << "\n<script type=\"text/javascript\">" << std::endl;
    ofs << "    function showOrHide(divname) {" << std::endl;
    ofs << "        div = document.getElementById(divname);" << std::endl;
    ofs << "        if(div.style.display == 'none')" << std::endl;
    ofs << "            div.style.display = 'block';" << std::endl;
    ofs << "        else" << std::endl;
    ofs << "            div.style.display = 'none';" << std::endl;
    ofs << "    }" << std::endl;
    ofs << "</script>" << std::endl;
}
