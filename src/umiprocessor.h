#ifndef UMI_PROCESSOR_H
#define UMI_PROCESSOR_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include "options.h"
#include "read.h"

/** Class to process umi */
class UmiProcessor{
    public:
        /** construct a UmiProcessor */
        UmiProcessor(Options* opt);

        /** destruct a UmiProcessor */
        ~UmiProcessor();

        /** find the umi according position arguments provided\n
         * and add umi to read name
         * @param r1 pointer to Read (read1)
         * @param r2 pointer to Read (read2)
         */
        void process(Read* r1, Read* r2 = NULL);
        /** add umi string to read 
         * @param r pointer to Read
         * @param umi umi string
         */ 
        void addUmiToName(Read* r, const std::string& umi);
    private:
        Options* mOptions;
        int UMI_LOC_INDEX1 = 1;///<
        int UMI_LOC_INDEX2 = 2;///<
        int UMI_LOC_READ1 = 3;///< 
        int UMI_LOC_READ2 = 4;///<
        int UMI_LOC_PER_INDEX = 5;///<
        int UMI_LOC_PER_READ = 6;///<
};

#endif
