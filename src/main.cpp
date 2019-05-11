#include <string>
#include "CLI.hpp"
#include "options.h"
#include "evaluator.h"
#include "processor.h"

int main(int argc, char** argv){
    std::string sysCMD = std::string(argv[0]) + " -h";
    if(argc == 1){
        std::system(sysCMD.c_str());
        return 0;
    }
    
    Options* opt = new Options();
    // I/O
    CLI::App app("program: " + std::string(argv[0]) + "\nversion: " + opt->version + "\nupdated: " + opt->compile);
    app.get_formatter()->column_width(50);
    CLI::Option* pin1 = app.add_option("-i", opt->in1, "read1 input file name")->required(true)->check(CLI::ExistingFile)->group("IO");
    app.add_option("-o", opt->out1, "read1 output file name")->required(true)->group("IO");
    CLI::Option* pin2 = app.add_option("-I", opt->in2, "read2 input file name")->needs(pin1)->check(CLI::ExistingFile)->group("IO");
    app.add_option("-O", opt->out2, "read2 output file name")->needs(pin2)->group("IO");
    app.add_option("--r1", opt->unpaired1, "output read1 whose mate failed QC")->group("IO");
    app.add_option("--r2", opt->unpaired2, "output read2 whose mate failed QC")->group("IO");
    app.add_option("--fo", opt->failedOut, "output failed QC reads")->group("IO");
    CLI::Option* pmerge = app.add_flag("-m", opt->mergePE.enabled, "merge overlapped readpair")->needs(pin2)->group("Merge");
    app.add_flag("--dm", opt->mergePE.discardUnmerged, "discard unmerged reads")->needs(pmerge)->group("Merge");
    app.add_option("--mo", opt->mergePE.out, "merged output")->needs(pmerge)->group("Merge");
    app.add_flag("--p64", opt->phred64, "input fastq is phred64")->group("IO");
    app.add_option("-z", opt->compression, "gzip output compress level")->check(CLI::Range(1, 9))->group("IO");
    app.add_flag("--iin", opt->interleavedInput, "input fastq interleaved")->excludes(pin2)->group("IO");
    // duplication
    CLI::Option* pdupana = app.add_flag("-d", opt->duplicate.enabled, "enable duplication analysis")->group("Duplication");
    app.add_option("--dl", opt->duplicate.keylen, "duplication analysis key length")->check(CLI::Range(12, 31))->needs(pdupana)->group("Duplication");
    app.add_option("--dh", opt->duplicate.histSize, "duplicate analysis hist size")->check(CLI::Range(1, 10000))->needs(pdupana)->group("Duplication");
    // adapter
    CLI::Option* pcutadapter = app.add_flag("-a", opt->adapter.enableTriming, "enable adapter trimming")->group("Adapter");
    app.add_option("--a1", opt->adapter.inputAdapterSeqR1, "adapter for read1")->needs(pcutadapter)->group("Adapter");
    app.add_option("--a2", opt->adapter.inputAdapterSeqR2, "adapter for read2")->needs(pcutadapter)->group("Adapter");
    app.add_flag("--dpea", opt->adapter.enableDetectForPE, "detect PE adapters")->needs(pin2)->group("Adapter");
    // trimming
    app.add_option("-f", opt->trim.front1, "bases trimmed in read1 front", true)->check(CLI::Range(0, 1000))->group("Trim");
    app.add_option("-t", opt->trim.tail1, "bases trimmed in read1 tail", true)->check(CLI::Range(0, 1000))->group("Trim");
    app.add_option("-b", opt->trim.maxLen1, "read1 max length allowed", true)->check(CLI::Range(0, 1000))->group("Trim");
    app.add_option("-F", opt->trim.front2, "bases trimmed in read2 front", true)->check(CLI::Range(0, 1000))->group("Trim");
    app.add_option("-T", opt->trim.tail2, "#bases trimmed in read2 tail", true)->check(CLI::Range(0, 1000))->group("Trim");
    app.add_option("-B", opt->trim.maxLen2, "read2 max length allowed", true)->check(CLI::Range(0, 1000))->group("Trim");
    // polyG tail trimming
    CLI::Option* ptrimG = app.add_flag("-g", opt->polyGTrim.enabled, "enable polyG trim")->group("PolyX");
    app.add_option("--pgm", opt->polyGTrim.minLen, "minimum length to detect polyG", true)->needs(ptrimG)->group("PolyX");
    // polyX tail trimming
    CLI::Option* ptrimX = app.add_flag("-x", opt->polyXTrim.enabled, "enable polyX trim")->group("PolyX");
    app.add_option("--pxm", opt->polyXTrim.minLen, "minimum length to detect polyG", true)->needs(ptrimX)->group("PolyX");
    // cutting by quality
    CLI::Option* pcutfront = app.add_flag("--cf", opt->qualitycut.enableFront, "window max mean quality to cut from 5'->3'")->group("Cut");
    CLI::Option* pcuttail = app.add_flag("--ct", opt->qualitycut.enableTail, "window max mean quality fo cut from 3'->5'")->group("Cut");
    CLI::Option* pcutright = app.add_flag("--cr", opt->qualitycut.enableRright, "drop window and right part when cut from 5'->3'")->group("Cut");
    app.add_option("-W", opt->qualitycut.windowSizeShared, "window size for cut slidind", true)->check(CLI::Range(0, 1000))->group("Cut");
    app.add_option("-M", opt->qualitycut.qualityShared, "min mean quality to drop window/bases", true)->check(CLI::Range(1, 36))->group("Cut");
    app.add_option("--cwf", opt->qualitycut.windowSizeFront, "window size for --cf'", true)->check(CLI::Range(0, 1000))->needs(pcutfront)->group("Cut");
    app.add_option("--cwt", opt->qualitycut.windowSizeTail, "window size for --ct'")->check(CLI::Range(0, 1000))->needs(pcuttail)->group("Cut");
    app.add_option("--cur", opt->qualitycut.windowSizeRight, "window size for --cr", true)->check(CLI::Range(0, 1000))->needs(pcutright)->group("Cut");
    app.add_option("--cqf", opt->qualitycut.qualityFront, "mean quality for --cf ", true)->check(CLI::Range(1, 36))->needs(pcutfront)->group("Cut");
    app.add_option("--cqt", opt->qualitycut.qualityTail, "mean quality for --ct", true)->check(CLI::Range(1, 36))->needs(pcuttail)->group("Cut");
    app.add_option("--cqr", opt->qualitycut.qualityRight, "mean quality for --cr", true)->check(CLI::Range(1, 36))->needs(pcuttail)->group("Cut");
    // quality filtering
    CLI::Option* pqfilter = app.add_flag("-q", opt->qualFilter.enabled, "enable quality filter")->group("Qual");
    app.add_option("-Q", opt->qualFilter.lowQualityLimit, "minimum ASCII Quality for qualified bases,", true)->needs(pqfilter)->check(CLI::Range(33, 75))->group("Qual");
    app.add_option("-U", opt->qualFilter.lowQualityRatio, "maximum low quality ratio allowed in one read", true)->check(CLI::Range(0, 1))->needs(pqfilter)->group("Qual");
    app.add_option("-N", opt->qualFilter.nBaseLimit, "maximum N bases allowed in one read", true)->needs(pqfilter)->group("Qual");
    app.add_option("-e", opt->qualFilter.averageQualityLimit, "average quality needed for one read")->needs(pqfilter)->group("Qual");
    // length filtering
    CLI::Option* plenfilter = app.add_flag("-l", opt->lengthFilter.enabled, "enable length filter")->group("Length");
    app.add_option("--lm", opt->lengthFilter.minReadLength, "min length required for a read", true)->check(CLI::Range(0, 1000))->needs(plenfilter)->group("Length");
    app.add_option("--lM", opt->lengthFilter.maxReadLength, "max length allowed for a read", true)->check(CLI::Range(0, 1000))->needs(plenfilter)->group("Length");
    // low complexity filtering
    CLI::Option* plowcfilter = app.add_flag("-y", opt->complexityFilter.enabled, "enable low complexity filter")->group("Complexity");
    app.add_option("-Y", opt->complexityFilter.threshold, "min complexity required for a read", true)->check(CLI::Range(0, 1))->needs(plowcfilter)->group("Complexity");
    // index filtering
    CLI::Option* pindexfilter = app.add_flag("--fi", opt->indexFilter.enabled, "enable index filtering")->group("Index");
    app.add_option("--fi1", opt->indexFilter.index1File, "index1 file to filter")->check(CLI::ExistingFile)->needs(pindexfilter)->group("Index");
    app.add_option("--fi2", opt->indexFilter.index2File, "index2 file to filetr")->check(CLI::ExistingFile)->needs(pindexfilter)->group("Index");
    app.add_option("--fit", opt->indexFilter.threshold, "max ed to validate index matcha", true)->check(CLI::Range(0, 10))->needs(pindexfilter)->group("Index");
    // base correction
    app.add_flag("-c", opt->correction.enabled, "enable base correction in PE reads")->group("Correction");
    app.add_option("--olr", opt->overlapRequire, "min overlap length needed for overlap analysis", true)->check(CLI::Range(0, 1000))->group("Correction");
    app.add_option("--odl", opt->overlapDiffLimit, "max ed to validate overlap", true)->check(CLI::Range(0, 10))->group("Correction");
    // umi processing
    CLI::Option* pumi = app.add_flag("-u", opt->umi.enabled, "enable UMI preprocess")->group("UMI");
    app.add_option("--uloc", opt->umi.location, "umi location,0->none,1->index1,2->index2,3->read1,4->read2,5->perindex,6->perread", true)->check(CLI::Range(1, 6))->needs(pumi)->group("UMI");
    app.add_option("--ulen", opt->umi.length, "umi length", true)->check(CLI::Range(0, 1000))->needs(pumi)->group("UMI");
    app.add_option("--uskip", opt->umi.skip, "bases to skip after umi", true)->check(CLI::Range(0, 1000))->needs(pumi)->group("UMI");
    app.add_flag("--udrop", opt->umi.dropOtherComment, "drop other comment information")->needs(pumi)->group("UMI");
    app.add_flag("--unt", opt->umi.notTrimRead, "do not trim reads")->needs(pumi)->group("UMI");
    // overrepresentation sequence analysis
    CLI::Option* pora = app.add_flag("--ora", opt->overRepAna.enabled, "enable ORA")->group("ORA");
    app.add_option("--osa", opt->overRepAna.sampling, "ORA sampling steps", true)->check(CLI::Range(1, 10000))->needs(pora)->group("ORA");
    // reporting 
    app.add_option("-J", opt->jsonFile, "json format report file", true)->group("Report");
    app.add_option("-H", opt->htmlFile, "html format report file", true)->group("Report");;
    // threading
    app.add_option("-w", opt->thread, "worker thread number", true)->check(CLI::Range(1, 16))->group("System");
    // output split
    CLI::Option* split_by_fn = app.add_flag("-s", opt->split.byFileNumber, "split output by file number")->excludes(pmerge)->group("Split");
    app.add_option("--sfn", opt->split.number, "total split output file number")->needs(split_by_fn)->group("Split");
    CLI::Option* split_by_ln = app.add_flag("-S", opt->split.byFileLines, "max line of each output file")->excludes(split_by_fn)->excludes(pmerge)->group("Split");
    app.add_option("--sfl", opt->split.size, "split output file line limit")->needs(split_by_ln)->group("Split");
    app.add_option("--sfd", opt->digits, "digits for sequential output filename", true)->check(CLI::Range(1, 10))->group("Split");
    // buffer size options
    app.add_option("--MR", opt->bufSize.maxPacksInReadPackRepo, "max packs in repo", true)->check(CLI::Range(1, 1000000))->group("System");
    app.add_option("--MP", opt->bufSize.maxReadsInPack, "max read/pairs in pack", true)->check(CLI::Range(1, 1000000))->group("System");
    app.add_option("--MM", opt->bufSize.maxPacksInMemory, "max packs in memory", true)->check(CLI::Range(1, 1000000))->group("System");
    // parse args
    CLI_PARSE(app, argc, argv);
    // update options
    opt->update();
    // validate options
    opt->validate();
    // evaluate read length
    Evaluator eva(opt);
    eva.evaluateReadLen();
    // evaluate read number
    eva.evaluateReadNum();
    if(opt->split.byFileNumber){
        opt->split.size = std::max(opt->est.readsNum / opt->split.number, 1);
        util::loginfo("total reds: " + std::to_string(opt->est.readsNum) + " split size: " + std::to_string(opt->split.size), opt->logmtx);
    }
    if(opt->overRepAna.enabled){
        eva.evaluateOverRepSeqs();
    }
    // evaluate adapter sequence
    if(opt->adapter.enableDetectForPE){
        eva.evaluateAdapterSeq(false);
        eva.evaluateAdapterSeq(true);
    }
    // setup processor
    Processor p(opt);
    p.process();
}
