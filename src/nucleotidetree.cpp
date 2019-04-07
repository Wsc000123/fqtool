#include "nucleotidetree.h"

NucleotideNode::NucleotideNode(){
    this->count = 0;
    this->base = 'N';
    std::memset(this->children, 0, sizeof(NucleotideNode*) * 8);
}

NucleotideNode::~NucleotideNode(){
    for(size_t i = 0; i < 8; ++i){
        if(this->children[i]){
            delete this->children[i];
            this->children[i] = NULL;
        }
    }
}

void NucleotideNode::dfs(){
    std::printf("%c", this->base);
    std::printf("%d", this->count);
    bool hasChild = false;
    for(size_t i = 0; i < 8; ++i){
        if(this->children[i]){
            this->children[i]->dfs();
            hasChild = true;
        }
    }
    if(!hasChild){
        printf("\n");
    }
}

NucleotideTree::NucleotideTree(){
    this->root = new NucleotideNode();
}

NucleotideTree::~NucleotideTree(){
    delete this->root;
}

void NucleotideTree::addSeq(const std::string& seq){
    NucleotideNode* curNode = this->root;
    for(size_t i = 0; i < seq.length(); ++i){
        if(seq[i] == 'N'){
            break;
        }
        int base = seq[i] & 0x07;
        if(curNode->children[base] == NULL){
            curNode->children[base] = new NucleotideNode();
            curNode->children[base]->base = seq[i];
        }
        curNode->children[base]->count++;
        curNode = curNode->children[base];
    }
}

std::string NucleotideTree::getDominantPath(bool& reachedLeaf){
    std::stringstream ss;
    const double RATIO_THRESHOLD = 0.95;
    const int NUM_THRESHOLD = 50;
    NucleotideNode* curNode = this->root;
    while(true){
        int total = 0;
        for(int i = 0; i < 8; ++i){
            if(curNode->children[i] != NULL){
                total += curNode->children[i]->count;
            }
        }
        if(total < NUM_THRESHOLD){
            break;
        }
        bool hasDominant = false;
        for(int i = 0; i < 8; ++i){
            if(curNode->children[i] == NULL){
                continue;
            }
            if(curNode->children[i]->count / (double)total >= RATIO_THRESHOLD){
                hasDominant = true;
                ss << curNode->children[i]->base;
                curNode = curNode->children[i];
                break;
            }
        }
        if(!hasDominant){
            reachedLeaf = false;
            break;
        }
    }
    return ss.str();
}
