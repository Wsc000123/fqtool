#ifndef HTML_UTIL_H
#define HTML_UTIL_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

/** some useful utilities for html output */
namespace htmlutil{
    /** output a row with two columns as a html table row to ofstream
     * @param ofs std::ofstream to output to
     * @param key first field of a row
     * @param val second field of a row
     */
    template<typename T>
    inline void outputTableRow(std::ofstream& ofs, const std::string& key, T val){
        std::stringstream ss;
        ss << val;
        ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + ss.str() + "</td></tr>\n";
    }

    /** format a big number with KMGTP units
     * @param number a non-negative number to be formatted
     * @return formatted number string
     */
    inline std::string formatNumber(const size_t& number){
        double num = (double)number;
        const std::string unit[6] = {"", "K", "M", "G", "T", "P"};
        int order = 0;
        while(number > 1000.0){
            order += 1;
            num /= 1000.0;
        }
        return std::to_string(num) + " " + unit[order];
    }

    /** get percents number of nuerator/denominator
     * @param numerator numerator
     * @param denominator denominator
     * @return percents number string
     */
    template<typename T>
    std::string getPercentsStr(T numerator, T denominator){
        if(denominator == 0){
            return "0.0";
        }
        return std::to_string((double)numerator * 100.0 / (double) denominator);
    }

    /** print JS of plotly
     * @param ofs reference of std::ofstream
     */ 
    inline void printJS(std::ofstream& ofs){
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

    /** get current system time
     * @return current system time
     */
    inline std::string getCurrentSystemTime(){
        auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        struct tm* ptm = localtime(&tt);
        char date[60] = {0};
        sprintf(date, "%d-%02d-%02d      %02d:%02d:%02d",
                (int)ptm->tm_year + 1900,(int)ptm->tm_mon + 1,(int)ptm->tm_mday,
                (int)ptm->tm_hour,(int)ptm->tm_min,(int)ptm->tm_sec);
        return std::string(date);
    }

    /** print CSS of HTML
     * @param ofs reference of std::ofstream
     */
    inline void printCSS(std::ofstream& ofs){
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

    /** print footer of HTML 
     * @param ofs reference of std::ofstream
     * @param comment comment string
     */
    inline void printFooter(std::ofstream& ofs, const std::string& commentStr){
        ofs << "\n</div>" << std::endl;
        ofs << "<div id='footer'> ";
        ofs << commentStr << std::endl;
        ofs << "date: " << getCurrentSystemTime() << " </div>";
        ofs << "</body></html>";
    }

    /** print header of HTML to ofs
     * @param ofs reference of std::ofstream
     * @param title title of HTML
     */    
    inline void printHeader(std::ofstream& ofs, const std::string& title){
        ofs << "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
        ofs << "<title>" << title << " </title>";
        printJS(ofs);
        printCSS(ofs);
        ofs << "</head>";
        ofs << "<body><div id='container'>";
    }
}

#endif
