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
    CTML::Document hReport;
    printHeader(hReport, "Fastq Preprocess Report");
    printSummary(hReport, result, preStats1, postStats1, preStats2, postStats2);
    CTML::Node preSection("div.section_div");
    CTML::Node preSectionTitle("div.section_title");
    preSectionTitle.SetAttribute("onclick", "showOrHide('before_filtering')");
    CTML::Node preSectionTitleLink("a", "Before filtering");
    preSectionTitleLink.SetAttribute("name", "summary");
    preSectionTitle.AppendChild(preSectionTitleLink);
    preSection.AppendChild(preSectionTitle);
    CTML::Node preSectionID("div#before_filtering");
    if(preStats1){
        std::vector<CTML::Node> v1 = preStats1->reportHtml("Before filtering", "read1");
        for(uint32_t i = 0; i < v1.size(); ++i){
            preSectionID.AppendChild(v1[i]);
        }
    }
    if(preStats2){
        std::vector<CTML::Node> v2 = preStats2->reportHtml("Before filtering", "read2");
        for(uint32_t i = 0; i < v2.size(); ++i){
            preSectionID.AppendChild(v2[i]);
        }
    }
    hReport.AppendNodeToBody(preSection);
    hReport.AppendNodeToBody(preSectionID);
   
    CTML::Node postSection("div.section_div");
    CTML::Node postSectionTitle("div.section_title");
    postSectionTitle.SetAttribute("onclick", "showOrHide('after_filtering')");
    CTML::Node postSectionTitleLink("a", "After filtering");
    postSectionTitleLink.SetAttribute("name", "summary");
    postSectionTitle.AppendChild(postSectionTitleLink);
    postSection.AppendChild(postSectionTitle);
    CTML::Node postSectionID("div#after_filtering");
    if(postStats1){
        std::vector<CTML::Node> v1 = postStats1->reportHtml("After filtering", "read1");
        for(uint32_t i = 0; i < v1.size(); ++i){
            postSectionID.AppendChild(v1[i]);
        }
    }
    if(postStats2){
        std::vector<CTML::Node> v2 = postStats2->reportHtml("After filtering", "read2");
        for(uint32_t i = 0; i < v2.size(); ++i){
            postSectionID.AppendChild(v2[i]);
        }
    }
    postSection.AppendChild(postSectionID);
    hReport.AppendNodeToBody(postSection);
    // software
    CTML::Node softwareSection("div#section_div");
    CTML::Node softwareSectionTitle("div.section_title");
    softwareSection.SetAttribute("onclick", "showOrHide('software')");
    CTML::Node softwareSectionTitleLink("a", "Software Environment");
    softwareSectionTitleLink.SetAttribute("name", "summary");
    softwareSectionTitle.AppendChild(softwareSectionTitleLink);
    softwareSection.AppendChild(softwareSectionTitle);
    CTML::Node softwareSectionID("div#software");
    CTML::Node softwareSectionTable("table.summary_table");
    softwareSectionTable.AppendChild(htmlutil::make2ColRowNode("Version", mOptions->version));
    softwareSectionTable.AppendChild(htmlutil::make2ColRowNode("Command", mOptions->command));
    softwareSectionTable.AppendChild(htmlutil::make2ColRowNode("CWD", mOptions->cwd));
    softwareSectionID.AppendChild(softwareSectionTable);
    hReport.AppendNodeToBody(softwareSection);
    hReport.AppendNodeToBody(softwareSectionID);
    // footer
    CTML::Node footer("div#footer", "Fqtool Report @ " + htmlutil::getCurrentSystemTime());
    hReport.AppendNodeToBody(footer);
    // write to file
    ofs << hReport.ToString(CTML::StringFormatting::SINGLE_LINE);
    ofs.close();
}

void HtmlReporter::printSummary(CTML::Document& d, FilterResult* fresult, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2){
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

    std::string sequencingInfo  = mOptions->isPaired()?"paired end":"single end";
    if(mOptions->isPaired()){
        sequencingInfo += " (" + std::to_string(preStats1->getCycles()) + " cycles + " + std::to_string(preStats2->getCycles()) + " cycles)";
    }else{
        sequencingInfo += " (" + std::to_string(preStats1->getCycles()) + " cycles)";
    }
    // h1 title node
    CTML::Node h1("h1");
    h1.SetAttribute("style", "text-align:left");
    CTML::Node h1Link("a");
    h1Link.SetAttribute("style", "color:#663355;text-decoration:none;").AppendText(mOptions->reportTitle);
    h1.AppendChild(h1Link);
    d.AppendNodeToHead(h1);
    // summary section
    CTML::Node summarySection("div.section_div");
    CTML::Node summaryTitle("div.section_title");
    summaryTitle.SetAttribute("onclick", "showOrHide('summary')");
    CTML::Node summaryTitleLink("a", "Summary");
    summaryTitleLink.SetAttribute("name", "summary");
    summaryTitle.AppendChild(summaryTitleLink);
    summarySection.AppendChild(summaryTitle);
    CTML::Node summaryID("div#summary");
    // summary->general
    CTML::Node generalSection("div.subsection_title", "General");
    generalSection.SetAttribute("onclick", "showOrHide('general')");
    CTML::Node generalID("div#general");
    CTML::Node generalTable("table.summary_table");
    generalTable.AppendChild(htmlutil::make2ColRowNode("Sequencing", sequencingInfo));
    if(mOptions->isPaired()){
        generalTable.AppendChild(htmlutil::make2ColRowNode("Insert Size Peak", mInsertSizePeak));
    }
    if(mOptions->adapter.enableTriming){
        if(!mOptions->adapter.detectedAdapterSeqR1.empty()){
            generalTable.AppendChild(htmlutil::make2ColRowNode("Detected Read1 Adapter", mOptions->adapter.detectedAdapterSeqR1));
        }
        if(!mOptions->adapter.detectedAdapterSeqR2.empty()){
            generalTable.AppendChild(htmlutil::make2ColRowNode("Detected Read2 Adapter", mOptions->adapter.detectedAdapterSeqR2));
        }
    }
    generalID.AppendChild(generalTable);
    summaryID.AppendChild(generalSection);
    summaryID.AppendChild(generalID);
    // summary->preFilter
    CTML::Node preFilterSection("div.subsection_title", "Before Filtering");
    preFilterSection.SetAttribute("onclick", "showOrHide('before_filtering_summary')");
    CTML::Node preFilterID("div#before_filtering_summary");
    CTML::Node preFilterTable("table.summary_table");
    preFilterTable.AppendChild(htmlutil::make2ColRowNode("Total Reads", preTotalReads));
    preFilterTable.AppendChild(htmlutil::make2ColRowNode("Total Bases", preTotalBases));
    preFilterTable.AppendChild(htmlutil::make2ColRowNode("Q20 Bases", std::to_string(preQ20Bases) + "(" + std::to_string(preQ20Rate * 100) + "%)"));
    preFilterTable.AppendChild(htmlutil::make2ColRowNode("Q30 Bases", std::to_string(preQ30Bases) + "(" + std::to_string(preQ30Rate * 100) + "%)"));
    preFilterTable.AppendChild(htmlutil::make2ColRowNode("GC Content", preGCRate));
    preFilterTable.AppendChild(htmlutil::make2ColRowNode("Read1 Mean Length", preRead1Length));
    if(mOptions->isPaired()){
        preFilterTable.AppendChild(htmlutil::make2ColRowNode("Read2 Mean Length", preRead2Length));
    }
    if(mOptions->adapter.enableTriming){
        size_t readWithAdapter = 0;
        for(auto& e: fresult->mAdapter1Count){
            readWithAdapter += e.second;
        }
        double readWithAdapterRate = mOptions->isPaired() ? readWithAdapter * 1.0 / preTotalReads * 2 : readWithAdapter * 1.0 / preTotalReads;
        preFilterTable.AppendChild(htmlutil::make2ColRowNode("Read1 Adapters Left", std::to_string(readWithAdapter) + "(" + std::to_string(readWithAdapterRate * 100) + "%)"));
        if(mOptions->isPaired()){
            readWithAdapter = 0;
            for(auto& e: fresult->mAdapter2Count){
                readWithAdapter += e.second;
            }
            readWithAdapterRate = mOptions->isPaired() ? readWithAdapter * 1.0 / preTotalReads * 2 : readWithAdapter * 1.0 / preTotalReads;
            preFilterTable.AppendChild(htmlutil::make2ColRowNode("Read2 Adapters Left", std::to_string(readWithAdapter) + "(" + std::to_string(readWithAdapterRate * 100) + "%)"));
        }
    }
    preFilterID.AppendChild(preFilterTable);
    summaryID.AppendChild(preFilterSection);
    summaryID.AppendChild(preFilterID);
    // summary->afterFilter
    CTML::Node postFilterSection("div.subsection_title", "After filtering");
    postFilterSection.SetAttribute("onclick", "showOrHide('after_filtering_summary')");
    CTML::Node postFilterID("div#after_filtering_summary");
    CTML::Node postFilterTable("table.summary_table");
    postFilterTable.AppendChild(htmlutil::make2ColRowNode("Total Reads", postTotalReads));
    postFilterTable.AppendChild(htmlutil::make2ColRowNode("Total Bases", postTotalBases));
    postFilterTable.AppendChild(htmlutil::make2ColRowNode("Q20 Bases", std::to_string(postQ20Bases) + "(" + std::to_string(postQ20Rate * 100) + "%)"));
    postFilterTable.AppendChild(htmlutil::make2ColRowNode("Q30 Bases", std::to_string(postQ30Bases) + "(" + std::to_string(postQ30Rate * 100) + "%)"));
    postFilterTable.AppendChild(htmlutil::make2ColRowNode("GC Content", postGCRate));
    postFilterTable.AppendChild(htmlutil::make2ColRowNode("Read1 Mean Length", postRead1Length));
    if(mOptions->isPaired()){
        postFilterTable.AppendChild(htmlutil::make2ColRowNode("Read2 Mean Length", postRead2Length));
    }
    postFilterID.AppendChild(postFilterTable);
    summaryID.AppendChild(postFilterSection);
    summaryID.AppendChild(postFilterID);
    // summary->filterresult
    CTML::Node filterResultSection("div.subsection_title", "Filtering Results");
    filterResultSection.SetAttribute("onclick", "showOrHide('filtering_result')");
    CTML::Node filterResultID("div#filtering_result");
    CTML::Node filterTable = fresult->reportHtmlBasic(preTotalBases, preTotalReads);
    filterResultID.AppendChild(filterTable);
    summaryID.AppendChild(filterResultSection);
    summaryID.AppendChild(filterResultID);
    d.AppendNodeToBody(summarySection);
    d.AppendNodeToBody(summaryID);
    // adapters 
    if(mOptions->adapter.enableTriming){
        d.AppendNodeToBody(fresult->reportAdaptersHtmlSummary(preTotalBases));
    }
    // duplication
    if(mOptions->duplicate.enabled){
        d.AppendNodeToBody(reportDuplication());
    }
}

CTML::Node HtmlReporter::reportDuplication(){
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
    delete[] x;
    delete[] percents;
    delete[] gc;

    CTML::Node dupSection("div.section_div");
    CTML::Node dupSectionTitle("div.section_title");
    dupSectionTitle.SetAttribute("onclick", "showOrHide('duplication')");
    CTML::Node dupSectionLink("a", "Duplication");
    dupSectionLink.SetAttribute("name", "summary");
    dupSectionTitle.AppendChild(dupSectionLink);
    dupSection.AppendChild(dupSectionTitle);

    CTML::Node dupSectionID("div#duplication");
    CTML::Node dupSectionFigID("div#duplication_figure");
    CTML::Node dupSectionFig("div.figure");
    dupSectionFig.SetAttribute("id", "plot_duplication").SetAttribute("style", "height:400px;");
    dupSectionFigID.AppendChild(dupSectionFig);
    dupSectionID.AppendChild(dupSectionFigID);

    CTML::Node dupSectionScript("script");
    dupSectionScript.SetAttribute("type", "text/javascript");
    dupSectionScript.AppendText(json_str);
    dupSection.AppendChild(dupSectionID);
    dupSection.AppendChild(dupSectionScript);
    return dupSection;
}

void HtmlReporter::printHeader(CTML::Document& d, const std::string& title){
    // add meta 
    CTML::Node meta("meta");
    meta.SetAttribute("http-equiv", "content-type");
    meta.SetAttribute("content", "text/html;charset=utf-8");
    meta.UseClosingTag(false);
    d.AppendNodeToHead(meta);
    d.AppendNodeToHead(CTML::Node("title", title));
    // add js
    CTML::Node jsSrc("script");
    jsSrc.SetAttribute("src", "https://cdn.plot.ly/plotly-latest.min.js");
    CTML::Node jsFunc("script");
    jsFunc.SetAttribute("type","text/javascript");
    jsFunc.AppendText("function showOrHide(divname) {\n");
    jsFunc.AppendText("  div = document.getElementById(divname);\n");
    jsFunc.AppendText("  if(div.style.display == 'none')\n");
    jsFunc.AppendText("     div.style.display = 'block';\n");
    jsFunc.AppendText("  else\n");
    jsFunc.AppendText("     div.style.display = 'none';\n");
    jsFunc.AppendText("}\n");
    d.AppendNodeToHead(jsSrc);
    d.AppendNodeToHead(jsFunc);
    // add css
    CTML::Node css("style");
    css.SetAttribute("type", "text/css");
    css.AppendText("td {border:1px solid #dddddd;padding:5px;font-size:12px;}\n");
    css.AppendText("table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:800px}\n");
    css.AppendText(".col1 {width:240px; font-weight:bold;}\n");
    css.AppendText(".adapter_col {width:500px; font-size:10px;}\n");
    css.AppendText("img {padding:30px;}\n");
    css.AppendText("#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}\n");
    css.AppendText("#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;");
    css.AppendText("font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helv  etica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}\n");
    css.AppendText("a:visited {color: #999999}\n");
    css.AppendText(".alignleft {text-align:left;}\n");
    css.AppendText(".alignright {text-align:right;}\n");
    css.AppendText(".figure {width:800px;height:600px;}\n");
    css.AppendText(".header {color:#ffffff;padding:1px;height:20px;background:#000000;}\n");
    css.AppendText(".section_title {color:#ffffff;font-size:20px;padding:5px;text-align:left;background:#663355; margin-top:10px;}\n");
    css.AppendText(".subsection_title {font-size:16px;padding:5px;margin-top:10px;text-align:left;color:#663355}\n");
    css.AppendText("#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}\n");
    css.AppendText(".menu_item {text-align:left;padding-top:5px;font-size:18px;}\n");
    css.AppendText(".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}\n");
    css.AppendText("#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}\n");
    css.AppendText("#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#663355;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}\n");
    css.AppendText(".kmer_table {text-align:center;font-size:8px;padding:2px;}\n");
    css.AppendText(".kmer_table td{text-align:center;font-size:8px;padding:0px;color:#ffffff}\n");
    css.AppendText(".sub_section_tips {color:#999999;font-size:10px;padding-left:5px;padding-bottom:3px;}\n");
    d.AppendNodeToHead(css);
}
