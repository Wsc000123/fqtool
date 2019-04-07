#include "overlapanalysis.h"

OverlapResult OverlapAnalysis::analyze(Read* r1, Read* r2, int overlapDiffLimit, int overlapRequire){
    return OverlapAnalysis::analyze(r1->seq, r2->seq, overlapDiffLimit, overlapRequire);
}

OverlapResult OverlapAnalysis::analyze(Seq& s1, Seq& s2, int overlapDiffLimit, int overlapRequire){
    Seq rs2 = ~s2;
    int len1 = s1.length();
    int len2 = rs2.length();
    const char* pstr1 = s1.seqStr.c_str();
    const char* pstr2 = rs2.seqStr.c_str();

    int complete_compare_require = 50; //something werid
    int overlapLen = 0;
    int offset = 0;
    int diff = 0;
    
    // TEMPLATE_LEN >= SEQ_LEN, so the 3' end of s1 and s2 do not have any adaptor sequences
    while(offset < len1 - overlapRequire){
        overlapLen = std::min(len1 - offset, len2);
        diff = 0;
        int i = 0;
        for(i = 0; i < overlapLen; ++i){
            if(pstr1[offset + i] != pstr2[i]){
               ++diff;
               if(diff >= overlapDiffLimit && i < complete_compare_require){
                   break;
               }
            }
        }
        if(diff < overlapDiffLimit || (diff >= overlapDiffLimit && i > complete_compare_require)){
            OverlapResult ovr;
            ovr.overlapped = true;
            ovr.offset = offset;
            ovr.overlapLen = overlapLen;
            ovr.diff = diff;
            return ovr;
        }
        ++offset;
    }

    // TEMPLATE_LEN < SEQ_LEN, so the 3' end of s1 and s2 3' endswith part of adapter sequences
    offset = 0;
    while(offset > overlapRequire - len2){
        overlapLen = std::min(len1, len2 - std::abs(offset));
        diff = 0;
        int i = 0;
        for(i = 0; i < overlapLen; ++i){
            if(pstr1[i] != pstr2[-offset + i]){
                ++diff;
                if(diff >= overlapDiffLimit && i < complete_compare_require){
                    break;
                }
            }
        }

        if(diff < overlapDiffLimit || (diff >= overlapDiffLimit && i > complete_compare_require)){
            OverlapResult ovr;
            ovr.overlapped = true;
            ovr.offset = offset;
            ovr.overlapLen = overlapLen;
            ovr.diff = diff;
            return ovr;
        }
        --offset;
    }
    OverlapResult ovr;
    ovr.overlapped = false;
    ovr.offset = ovr.overlapLen = ovr.diff = 0;
    return ovr;
}
