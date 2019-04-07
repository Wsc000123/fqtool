#ifndef FILTER_H
#define FILTER_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "options.h"
#include "common.h"
#include "read.h"

/** Class to do fastq read filter by various standards and methods */
class Filter{
    public:
        Options* mOptions;   ///< Pointer to Options object
        /** Construct a Filter object, negative parameter will turn the corresponding filterr
         * @param opt pointer to Options
         */
        Filter(Options* opt) : mOptions(opt) { }

        /** Destroy a Filter object */
        ~Filter() = default;

        /** Test whether a Read pass filter
         * @param r pointer to a Read object
         * @return 0 if passed
         */
        int passFilter(Read* r);
        /** Test whether a Read pass low complexity filter
         * @param r pointer to a Read object
         * @return true if r is not low complexity
         */
        bool passLowComplexityFliter(Read* r);

        /** Trim and cut a Read by various methods
         * @param r pointer to a Read object
         * @param forceFrontCut force cut length from 5'end, will always happen
         * @param forceTailCut force cut length from 3'end, will always happen
         * @return a Read trim and cut properly
         */
        Read* trimAndCut(Read* r, int forceFrontCut, int forceTailCut);
        
        /** Filter a Read by its first index
         * @param r pointer to a Read object
         * @return true if r's index is succesfully matched in the indexBlackList
         */
        bool filterByIndex(Read* r);
        
        /** Filter a pair of Read by their first index
         * @param r1 pointer to a Read object
         * @param r2 pointer to a Read object
         * @return true if the first index of r1 or r2 succesfully matched in the indexBlackList
         */
        bool filterByIndex(Read* r1, Read* r2);

        /** match a sequnce (target) against a series of sequences (list)
         * @param list a vector of strings to be matched against
         * @param target a string to be tested whether it will match any string in list 
         * @param threshold maximum character difference allowed for a successful match
         * @return true if any successful match occurs between target and string in list
         */
        bool match(const std::vector<std::string>& list, const std::string& target, int threshold);
};

#endif
