#ifndef HTML_REPORTER_H
#define HTML_REPORTER_H

#include <iostream>
#include <fstream>
#include <string>
#include "stats.h"
#include "options.h"
#include "htmlutil.h"
#include "filterresult.h"

/** class to do json report of qc summary information */
class HtmlReporter{
    public:
        Options* mOptions;///< pointer to Options object
        size_t* mDupHist;///< duplication array
        double* mDupMeanGC;///< duplication mean gc array
        double mDupRate;///< duplication rate
        long* mInsertHist;///< insertsize array
        int mInsertSizePeak;///< insert size peak

    public:
        /** construct a HtmlReporter object
         * @param opt pointer to Options object
         */
        HtmlReporter(Options* opt);
        
        /** destroy a HtmlReporter object */
        ~HtmlReporter();

        /** set duplication statistical parameters
         * @param dupHist duplication statistics array
         * @param dupMeanGC duplication mean gc array
         * @param dupRate duplication rate
         */
        void setDupHist(size_t* dupHist, double* dupMeanGC, double dupRate);
        
        /** set insert statistical parameters
         * @param insertHist insert size array
         * @param insertSizePeak insert size peak
         */
        void setInsertHist(long* insertHist, int insertSizePeak);

        /** generate summary info
         * @param ofs reference of std::ofstream
         * @param fresult pointer to FilterResult object
         * @param preStats1 pointer to Stats object
         * @param preStats2 pointer to Stats object
         * @param postStats1 pointer to Stats object
         * @param postStats2 pointer to Stats object
         */
        void printSummary(std::ofstream& ofs, FilterResult* fresult, Stats* preStats1, Stats* postStats1, Stats* preStats2 = NULL, Stats* postStats2 = NULL);
        
        /** generate html report
         * @param fresult pointer to FilterResult object
         * @param preStats1 pointer to Stats object
         * @param preStats2 pointer to Stats object
         * @param postStats1 pointer to Stats object
         * @param postStats2 pointer to Stats object
         */
        void report(FilterResult* fresult, Stats* preStats1, Stats* postStats1, Stats* preStats2 = NULL, Stats* postStats2 = NULL);
};

#endif
