#ifndef JSON_UTIL_H
#define JSON_UTIL_H

#include <fstream>
#include <string>

/** json output utilities */
namespace jsonutil{
    /** output a key-val record to ofstream in json format
     * @param ofs std::ofstream to output to
     * @param padding padding before each record
     * @param key the key of a json record
     * @param val the value of a json record
     */
    template<typename T>
    inline void writeRecord(std::ofstream& ofs, const std::string& padding, const std::string& key, const T& val){
        ofs << padding << "\t\"" << key << "\": " << val << ",\n";
    }
}

#endif
