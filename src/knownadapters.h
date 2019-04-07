#ifndef KNOWN_ADAPTERS_H
#define KNOWN_ADAPTERS_H

#include <map>
#include <string>

/** Get a map of Illumina some known adapters
 * @return map of "adapter seq" -> "adapter name"
 */
inline std::map<std::string, std::string> getKnownAdapter() {
    std::map<std::string, std::string> knownAdapters;
    knownAdapters["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"] = "Illumina TruSeq Adapter Read 1";
    knownAdapters["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"] = "Illumina TruSeq Adapter Read 2";
    knownAdapters["GATCGTCGGACTGTAGAACTCTGAACGTGTAGA"] = "Illumina Small RNA Adapter Read 2";
    return knownAdapters;
}

#endif
