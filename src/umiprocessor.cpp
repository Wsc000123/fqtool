#include "umiprocessor.h"

UmiProcessor::UmiProcessor(Options* opt){
    mOptions = opt;
}

UmiProcessor::~UmiProcessor(){
}

void UmiProcessor::process(Read* r1, Read* r2){
    if(!mOptions->umi.enabled){
        return;
    }
    int loc = mOptions->umi.location;
    int len = mOptions->umi.length;
    std::string umi = " OX:Z:";
    std::string qua = " BZ:Z:";
    switch(loc){
        case UMI_LOC_INDEX1:
            umi += r1->firstIndex();
            break;
        case UMI_LOC_INDEX2:
            if(r2){
                umi += r2->firstIndex();
            }
            break;
        case UMI_LOC_READ1:
            umi += r1->seq.seqStr.substr(0, std::min(r1->length(), len));
            qua += r1->quality.substr(0, std::min(r1->length(), len));
            if(!mOptions->umi.notTrimRead){
                r1->trimFront(len + mOptions->umi.skip);
            }
            break;
        case UMI_LOC_READ2:
            if(r2){
                umi += r2->seq.seqStr.substr(0, std::min(r2->length(), len));
                qua += r2->quality.substr(0, std::min(r1->length(), len));
                if(!mOptions->umi.notTrimRead){
                    r2->trimFront(len + mOptions->umi.skip);
                }
            }
            break;
        case UMI_LOC_PER_INDEX:
            umi += r1->firstIndex();
            if(r2){
                umi.append("-" + r2->firstIndex());
            }
            break;
        case UMI_LOC_PER_READ:
            umi += r1->seq.seqStr.substr(0, std::min(r1->length(), len));
            qua += r1->quality.substr(0, std::min(r1->length(), len));
            if(!mOptions->umi.notTrimRead){
                r1->trimFront(len + mOptions->umi.skip);
            }
            if(r2){
                umi.append("-" + r2->seq.seqStr.substr(0, std::min(r2->length(), len)));
                if(!mOptions->umi.notTrimRead){
                    r2->trimFront(len + mOptions->umi.skip);
                }
                qua.append("-" + r2->quality.substr(0, std::min(r1->length(), len)));
            }
            break;
        default:
            break;
    }
    std::string tag = umi;
    if(tag.length() > 6 && qua.length() > 6){
        tag.append(qua);
    }
    if(r1 && tag.length() > 6){
        addTagToName(r1, tag);
    }
    if(r2 && tag.length() > 6){
        addTagToName(r2, tag);
    }
}

void UmiProcessor::addTagToName(Read* r, const std::string& tag){
    std::string::size_type pos = r->name.find_first_of(" ");
    if(pos == std::string::npos){
        r->name = r->name + tag;
    }else{
        if(mOptions->umi.dropOtherComment){
            r->name = r->name.substr(0, pos) + tag;
        }else{
            r->name =r->name.substr(0, pos) + tag + r->name.substr(pos, r->name.length() - pos);
        }
    }
}
