#include "basecorrector.h"

BaseCorrector::BaseCorrector(){
}

BaseCorrector::~BaseCorrector(){
}

int BaseCorrector::correctByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr){
    OverlapResult ov = OverlapAnalysis::analyze(r1, r2);
    return correctByOverlapAnalysis(r1, r2, fr, ov);
}

int BaseCorrector::correctByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr, OverlapResult& ov){
    if(ov.diff == 0 || ov.diff > 5){
        return 0;
    }
    int ol = ov.overlapLen;
    int start1 = std::max(0, ov.offset);
    int start2 = r2->length() - std::max(0, -ov.offset) - 1;

    const char* seq1 = r1->seq.seqStr.c_str();
    const char* seq2 = r2->seq.seqStr.c_str();
    const char* qual1 = r1->quality.c_str();
    const char* qual2 = r2->quality.c_str();

    const char GOOD_QUAL = util::num2qual(30);
    const char BAD_QUAL = util::num2qual(14);

    int corrected = 0;
    int uncorrected = 0;
    bool r1Corrected = false;
    bool r2Corrected = false;

    for(int i = 0; i < ol; ++i){
        int p1 = start1 + i;
        int p2 = start2 - i;

        if(seq1[p1] != util::complement(seq2[p2])){
            if(qual1[p1] >= GOOD_QUAL && qual2[p2] <= BAD_QUAL){
                r2->seq.seqStr[p2] = util::complement(seq1[p1]);
                r2->quality[p2] = qual1[p1];
                ++corrected;
                r2Corrected = true;
                if(fr){
                    fr->addCorrection(seq2[p2], util::complement(seq1[p1]));
                }
            }else if(qual2[p2] >= GOOD_QUAL && qual1[p1] <= BAD_QUAL){
                r1->seq.seqStr[p1] = util::complement(seq2[p2]);
                r1->quality[p1] = qual2[p2];
                ++corrected;
                r1Corrected = true;
                if(fr){
                    fr->addCorrection(seq1[p1], util::complement(seq2[p2]));
                }
            }else{
                ++uncorrected;
            }
        }
    }

    if(corrected > 0 && fr){
        if(r1Corrected && r2Corrected){
            fr->incCorrectedReads(2);
        }else{
            fr->incCorrectedReads(1);
        }
    }
    return corrected;
}
