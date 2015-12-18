// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// vector_join
std::string vector_join(const std::vector<std::string>& v, const std::string& token);
RcppExport SEXP mango_vector_join(SEXP vSEXP, SEXP tokenSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type token(tokenSEXP);
    __result = Rcpp::wrap(vector_join(v, token));
    return __result;
END_RCPP
}
// string_split
std::vector<std::string> string_split(const std::string& s, const std::string& delimiter);
RcppExport SEXP mango_string_split(SEXP sSEXP, SEXP delimiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const std::string& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type delimiter(delimiterSEXP);
    __result = Rcpp::wrap(string_split(s, delimiter));
    return __result;
END_RCPP
}
// parseFastq_gzip
std::vector< int > parseFastq_gzip(std::string fastq1, std::string fastq2, std::string basename, int minlength, int maxlength, bool keepempty, bool verbose, std::string linker1, std::string linker2);
RcppExport SEXP mango_parseFastq_gzip(SEXP fastq1SEXP, SEXP fastq2SEXP, SEXP basenameSEXP, SEXP minlengthSEXP, SEXP maxlengthSEXP, SEXP keepemptySEXP, SEXP verboseSEXP, SEXP linker1SEXP, SEXP linker2SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type fastq1(fastq1SEXP);
    Rcpp::traits::input_parameter< std::string >::type fastq2(fastq2SEXP);
    Rcpp::traits::input_parameter< std::string >::type basename(basenameSEXP);
    Rcpp::traits::input_parameter< int >::type minlength(minlengthSEXP);
    Rcpp::traits::input_parameter< int >::type maxlength(maxlengthSEXP);
    Rcpp::traits::input_parameter< bool >::type keepempty(keepemptySEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< std::string >::type linker1(linker1SEXP);
    Rcpp::traits::input_parameter< std::string >::type linker2(linker2SEXP);
    __result = Rcpp::wrap(parseFastq_gzip(fastq1, fastq2, basename, minlength, maxlength, keepempty, verbose, linker1, linker2));
    return __result;
END_RCPP
}
// parseFastq
std::vector< int > parseFastq(std::string fastq1, std::string fastq2, std::string basename, int minlength, int maxlength, bool keepempty, bool verbose, std::string linker1, std::string linker2);
RcppExport SEXP mango_parseFastq(SEXP fastq1SEXP, SEXP fastq2SEXP, SEXP basenameSEXP, SEXP minlengthSEXP, SEXP maxlengthSEXP, SEXP keepemptySEXP, SEXP verboseSEXP, SEXP linker1SEXP, SEXP linker2SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type fastq1(fastq1SEXP);
    Rcpp::traits::input_parameter< std::string >::type fastq2(fastq2SEXP);
    Rcpp::traits::input_parameter< std::string >::type basename(basenameSEXP);
    Rcpp::traits::input_parameter< int >::type minlength(minlengthSEXP);
    Rcpp::traits::input_parameter< int >::type maxlength(maxlengthSEXP);
    Rcpp::traits::input_parameter< bool >::type keepempty(keepemptySEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< std::string >::type linker1(linker1SEXP);
    Rcpp::traits::input_parameter< std::string >::type linker2(linker2SEXP);
    __result = Rcpp::wrap(parseFastq(fastq1, fastq2, basename, minlength, maxlength, keepempty, verbose, linker1, linker2));
    return __result;
END_RCPP
}
// mergeTwoBam
void mergeTwoBam(std::string inputBam1, std::string inputBam2, std::string outputBam);
RcppExport SEXP mango_mergeTwoBam(SEXP inputBam1SEXP, SEXP inputBam2SEXP, SEXP outputBamSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type inputBam1(inputBam1SEXP);
    Rcpp::traits::input_parameter< std::string >::type inputBam2(inputBam2SEXP);
    Rcpp::traits::input_parameter< std::string >::type outputBam(outputBamSEXP);
    mergeTwoBam(inputBam1, inputBam2, outputBam);
    return R_NilValue;
END_RCPP
}
// buildBedpe
void buildBedpe(std::string sam1, std::string sam2, std::string bedpefile);
RcppExport SEXP mango_buildBedpe(SEXP sam1SEXP, SEXP sam2SEXP, SEXP bedpefileSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type sam1(sam1SEXP);
    Rcpp::traits::input_parameter< std::string >::type sam2(sam2SEXP);
    Rcpp::traits::input_parameter< std::string >::type bedpefile(bedpefileSEXP);
    buildBedpe(sam1, sam2, bedpefile);
    return R_NilValue;
END_RCPP
}
// buildBedpefromBam
void buildBedpefromBam(std::string bam1, std::string bam2, std::string bedpefile);
RcppExport SEXP mango_buildBedpefromBam(SEXP bam1SEXP, SEXP bam2SEXP, SEXP bedpefileSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type bam1(bam1SEXP);
    Rcpp::traits::input_parameter< std::string >::type bam2(bam2SEXP);
    Rcpp::traits::input_parameter< std::string >::type bedpefile(bedpefileSEXP);
    buildBedpefromBam(bam1, bam2, bedpefile);
    return R_NilValue;
END_RCPP
}
// findPairs
void findPairs(std::string overlapfile, std::string petpairsfile, std::string interactionfile, std::string peakscount, std::string peaksfileslopdepth);
RcppExport SEXP mango_findPairs(SEXP overlapfileSEXP, SEXP petpairsfileSEXP, SEXP interactionfileSEXP, SEXP peakscountSEXP, SEXP peaksfileslopdepthSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type overlapfile(overlapfileSEXP);
    Rcpp::traits::input_parameter< std::string >::type petpairsfile(petpairsfileSEXP);
    Rcpp::traits::input_parameter< std::string >::type interactionfile(interactionfileSEXP);
    Rcpp::traits::input_parameter< std::string >::type peakscount(peakscountSEXP);
    Rcpp::traits::input_parameter< std::string >::type peaksfileslopdepth(peaksfileslopdepthSEXP);
    findPairs(overlapfile, petpairsfile, interactionfile, peakscount, peaksfileslopdepth);
    return R_NilValue;
END_RCPP
}
// splitBedbyChrom
std::vector<std::string> splitBedbyChrom(std::string bedfile, std::string outnamebase);
RcppExport SEXP mango_splitBedbyChrom(SEXP bedfileSEXP, SEXP outnamebaseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type bedfile(bedfileSEXP);
    Rcpp::traits::input_parameter< std::string >::type outnamebase(outnamebaseSEXP);
    __result = Rcpp::wrap(splitBedbyChrom(bedfile, outnamebase));
    return __result;
END_RCPP
}
// makeDistanceFile
void makeDistanceFile(std::string bedpefilesortrmdup, std::string distancefile, int mindist, int maxdist);
RcppExport SEXP mango_makeDistanceFile(SEXP bedpefilesortrmdupSEXP, SEXP distancefileSEXP, SEXP mindistSEXP, SEXP maxdistSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type bedpefilesortrmdup(bedpefilesortrmdupSEXP);
    Rcpp::traits::input_parameter< std::string >::type distancefile(distancefileSEXP);
    Rcpp::traits::input_parameter< int >::type mindist(mindistSEXP);
    Rcpp::traits::input_parameter< int >::type maxdist(maxdistSEXP);
    makeDistanceFile(bedpefilesortrmdup, distancefile, mindist, maxdist);
    return R_NilValue;
END_RCPP
}
// joinchromfiles
void joinchromfiles(std::vector<std::string> sortedchromfiles, std::string bedpefilesort);
RcppExport SEXP mango_joinchromfiles(SEXP sortedchromfilesSEXP, SEXP bedpefilesortSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type sortedchromfiles(sortedchromfilesSEXP);
    Rcpp::traits::input_parameter< std::string >::type bedpefilesort(bedpefilesortSEXP);
    joinchromfiles(sortedchromfiles, bedpefilesort);
    return R_NilValue;
END_RCPP
}
// DeterminePeakDepthsC
void DeterminePeakDepthsC(std::string temppeakoverlap, std::string peaksfileslopdepth);
RcppExport SEXP mango_DeterminePeakDepthsC(SEXP temppeakoverlapSEXP, SEXP peaksfileslopdepthSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type temppeakoverlap(temppeakoverlapSEXP);
    Rcpp::traits::input_parameter< std::string >::type peaksfileslopdepth(peaksfileslopdepthSEXP);
    DeterminePeakDepthsC(temppeakoverlap, peaksfileslopdepth);
    return R_NilValue;
END_RCPP
}
// removeDups
std::vector< std::string > removeDups(std::string bedpein, std::string outnamebase, double distancesplit);
RcppExport SEXP mango_removeDups(SEXP bedpeinSEXP, SEXP outnamebaseSEXP, SEXP distancesplitSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type bedpein(bedpeinSEXP);
    Rcpp::traits::input_parameter< std::string >::type outnamebase(outnamebaseSEXP);
    Rcpp::traits::input_parameter< double >::type distancesplit(distancesplitSEXP);
    __result = Rcpp::wrap(removeDups(bedpein, outnamebase, distancesplit));
    return __result;
END_RCPP
}
// splitBedpe
std::vector<std::string> splitBedpe(std::string bedpein, std::string outnamebase, bool printreads, bool printpets, bool skipstars, bool skipinter);
RcppExport SEXP mango_splitBedpe(SEXP bedpeinSEXP, SEXP outnamebaseSEXP, SEXP printreadsSEXP, SEXP printpetsSEXP, SEXP skipstarsSEXP, SEXP skipinterSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type bedpein(bedpeinSEXP);
    Rcpp::traits::input_parameter< std::string >::type outnamebase(outnamebaseSEXP);
    Rcpp::traits::input_parameter< bool >::type printreads(printreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type printpets(printpetsSEXP);
    Rcpp::traits::input_parameter< bool >::type skipstars(skipstarsSEXP);
    Rcpp::traits::input_parameter< bool >::type skipinter(skipinterSEXP);
    __result = Rcpp::wrap(splitBedpe(bedpein, outnamebase, printreads, printpets, skipstars, skipinter));
    return __result;
END_RCPP
}
// buildTagAlign
void buildTagAlign(std::string bedpefile, std::string TagAlignfile);
RcppExport SEXP mango_buildTagAlign(SEXP bedpefileSEXP, SEXP TagAlignfileSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type bedpefile(bedpefileSEXP);
    Rcpp::traits::input_parameter< std::string >::type TagAlignfile(TagAlignfileSEXP);
    buildTagAlign(bedpefile, TagAlignfile);
    return R_NilValue;
END_RCPP
}
// external_sort
void external_sort(std::string inputfile, std::string outputfile);
RcppExport SEXP mango_external_sort(SEXP inputfileSEXP, SEXP outputfileSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type inputfile(inputfileSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputfile(outputfileSEXP);
    external_sort(inputfile, outputfile);
    return R_NilValue;
END_RCPP
}
// everyotherline
void everyotherline(std::string overlapin, std::string overlapout);
RcppExport SEXP mango_everyotherline(SEXP overlapinSEXP, SEXP overlapoutSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type overlapin(overlapinSEXP);
    Rcpp::traits::input_parameter< std::string >::type overlapout(overlapoutSEXP);
    everyotherline(overlapin, overlapout);
    return R_NilValue;
END_RCPP
}
// AddQvals
void AddQvals(std::string interactionfile, std::string interactionfilefinal, std::vector<double> Q, double maxPval);
RcppExport SEXP mango_AddQvals(SEXP interactionfileSEXP, SEXP interactionfilefinalSEXP, SEXP QSEXP, SEXP maxPvalSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type interactionfile(interactionfileSEXP);
    Rcpp::traits::input_parameter< std::string >::type interactionfilefinal(interactionfilefinalSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type Q(QSEXP);
    Rcpp::traits::input_parameter< double >::type maxPval(maxPvalSEXP);
    AddQvals(interactionfile, interactionfilefinal, Q, maxPval);
    return R_NilValue;
END_RCPP
}
