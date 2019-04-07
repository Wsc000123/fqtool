#ifndef NUCLEOTIDETREE_H
#define NUCLEOTIDETREE_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <memory>
#include <sstream>

/** class represent a nucleotide in a sequence */
class NucleotideNode{
    public:
        int count; ///< # of this kind of nucleotides at this node
        char base; ///< base name, can only be aAtTcCgG
        NucleotideNode* children[8]; ///< children/next nucleotide node, 
                                     ///< aA->children[1], cT->children[4], 
                                     ///< gG->children[7], cC->children[3],
                                     ///< by base & 0x7 bitwise and operation 
    
    public:
        /** Construct a NucldotideNode, set count = 0, base = 'N' and memset children[8] all to 0 */
        NucleotideNode();
        /** Destroy a NucleotideNode, free all memory of its children if possible */
        ~NucleotideNode();
        /** Depth-First-Show NucleotideNodes, print base and count of each node by DFS */
        void dfs();
};

/** class represent a series of nucleotide sequences started at the same 5' position 
 * the root is just a virtual nucleotide N which is not part of these nucleotide sequences
 */
class NucleotideTree{
    public:
        NucleotideNode* root; ///< the root node with N as its base and only acts as root

    public:
        /** Construct a NucleotideTree with only one node, the root */
        NucleotideTree();

        /** Destroy a NucleotideTree, delete the root */
        ~NucleotideTree();
        
        /** Add a nucleotide sequence to this NucleotideTree
         * from the first Nucleotide in the sequence to the first 'N' base or the end
         * add Nucldotide as NucleotideNode sequentially to the children based on 0x07 bitwise and operation
         * @param seq nucleotide sequence to be added
         */
        void addSeq(const std::string& seq);
        
        /** Get the dominant path(most common sequence) from the NucleotideTree
         * a dominant path will return only if each node in the path has at least total 50 nucleotides
         * and a single base dominate more than 0.95 nucleotides at this node, this 
         * base is the common base at this position then
         * @param reachedLeaf no dominant path found, and haven't reach the leaf
         * @return the longest dominant path if exists
         */
        std::string getDominantPath(bool& reachedLeaf);
};

#endif
