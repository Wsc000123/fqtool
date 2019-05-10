#ifndef HTML_UTIL_H
#define HTML_UTIL_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctml.hpp>

/** some useful utilities for html output */
namespace htmlutil{
    /** output a row with two columns as a html table row to ofstream
     * @param ofs std::ofstream to output to
     * @param key first field of a row
     * @param val second field of a row
     */
    template<typename T>
    inline CTML::Node make2ColRowNode(const std::string& key, T val){
        std::stringstream ss;
        ss << val;
        CTML::Node row("tr");
        CTML::Node col1("td.col1", key);
        CTML::Node col2("td.col2", ss.str());
        row.AppendChild(col1).AppendChild(col2);
        return row;
    }

    /** format a big number with KMGTP units
     * @param number a non-negative number to be formatted
     * @return formatted number string
     */
    inline std::string formatNumber(const size_t& number){
        double num = (double)number;
        const std::string unit[6] = {"", "K", "M", "G", "T", "P"};
        int order = 0;
        while(num > 1000.0){
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
}

#endif
