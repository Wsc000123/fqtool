#ifndef COMMON_H
#define COMMON_H

#pragma pack()

namespace COMMONCONST{
    // different filtering results, bigger number means worse
    // if r1 and r2 are both failed, then the bigger one of the two results will be recorded
    static const int PASS_FILTER = 0;           ///< 000000
    static const int FAIL_POLY_X = 4;           ///< 000100
    static const int FAIL_OVERLAP = 8;          ///< 001000
    static const int FAIL_N_BASE = 12;          ///< 001100
    static const int FAIL_LENGTH = 16;          ///< 010000
    static const int FAIL_TOO_LONG = 17;        ///< 010001
    static const int FAIL_QUALITY = 20;         ///< 010100
    static const int FAIL_COMPLEXITY = 24;      ///< 011000
    
    // how many types in total we support
    static const int FILTER_RESULT_TYPES = 32;  ///< 100000

    const static char* FAILED_TYPES[FILTER_RESULT_TYPES] = {
        "passed", "", "", "",
        "failed_polyx_filter", "", "", "",
        "failed_bad_overlap", "", "", "",
        "failed_too_many_n_bases", "", "", "",
        "failed_too_short", "failed_too_long", "", "",
        "failed_quality_filter", "", "", "",
        "failed_low_complexity", "", "", "",
        "", "", "", ""
    };
}

#endif
