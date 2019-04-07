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
    std::string umi;
    if(loc  == UMI_LOC_INDEX1){
        umi = r1->firstIndex();
    }else if(loc == UMI_LOC_INDEX2 && r2){
        umi = r2->firstIndex();
    }else if(loc == UMI_LOC_READ1){
        umi = r1->seq.seqStr.substr(0, std::min(r1->length(), len));
    }else if(loc == UMI_LOC_READ2 && r2){
        umi = r2->seq.seqStr.substr(0, std::min(r2->length(), len));
    }else if(loc == UMI_LOC_PER_INDEX){
        std::string umiMerged = r1->firstIndex();
        if(r2){
            umiMerged = umiMerged + "_" + r2->lastIndex();
        }
        addUmiToName(r1, umiMerged);
        if(r2){
            addUmiToName(r2, umiMerged);
        }
    }else if(loc == UMI_LOC_PER_READ){
        std::string umi1 = r1->seq.seqStr.substr(0, std::min(r1->length(), len));
        std::string umiMerged = umi1;
        r1->trimFront(umi.length() + mOptions->umi.skip);
        if(r2){
            std::string umi2 = r2->seq.seqStr.substr(0, std::min(r2->length(), len));
            umiMerged = umiMerged + "_" + umi2;
            r2->trimFront(umi2.length() + mOptions->umi.skip);
        }
        addUmiToName(r1, umiMerged);
        if(r2){
            addUmiToName(r2, umiMerged);
        }
    }

    if(loc != UMI_LOC_PER_READ && loc != UMI_LOC_PER_INDEX){
        if(r1 && !umi.empty()){
            addUmiToName(r1, umi);
        }
        if(r2 && !umi.empty()){
            addUmiToName(r2, umi);
        }
    }
}

void UmiProcessor::addUmiToName(Read* r, const std::string& umi){
    std::string tag;
    if(mOptions->umi.prefix.empty()){
        tag = ":" + umi;
    }else{
        tag = ":" + mOptions->umi.prefix + "_" + umi;
    }
    std::string::size_type pos = r->name.find_first_of(" ");
    if(pos == std::string::npos){
        r->name = r->name + tag;
    }else{
        r->name = r->name.substr(0, pos) + tag + r->name.substr(pos, r->name.length() - pos);
    }
}
