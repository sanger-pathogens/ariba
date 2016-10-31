#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include "Python.h"
#include "fml.h"


class Assembly
{
public:
    Assembly(int numberOfUnitigs, fml_utg_t *unitigs, unsigned short minCountIn, unsigned short overlapIn, std::string namePrefixIn);
    void printStats(std::ostream& outStream) const;
    void toFile(std::ostream& outStream) const;

    uint32_t numberOfContigs;
    unsigned short minCount;
    unsigned short overlap;
    std::string namePrefix;
    uint32_t longestContig;
    float meanLength;
    std::vector<std::string> sequences;
};


Assembly::Assembly(int numberOfUnitigs, fml_utg_t *unitigs, unsigned short minCountIn, unsigned short overlapIn, std::string namePrefixIn)
{
    numberOfContigs = numberOfUnitigs;
    minCount = minCountIn;
    overlap = overlapIn;
    namePrefix = namePrefixIn;
    longestContig = 0;
    meanLength = 0;
    uint32_t lengthSum = 0;

    if (numberOfContigs == 0)
    {
        return;
    }

    for (uint32_t i = 0; i < numberOfContigs; ++i)
    {
        const fml_utg_t *unitig = &unitigs[i];
        sequences.push_back(unitig->seq);
        longestContig = ((uint32_t) unitig->len) > longestContig ? unitig->len : longestContig;
        lengthSum += unitig->len;
    }

    meanLength = 1.0 * lengthSum / numberOfContigs;
}


void Assembly::printStats(std::ostream& outStream) const
{
    outStream << overlap << '\t' << minCount << '\t' << numberOfContigs << '\t' << meanLength << '\t' << longestContig << std::endl;
}


void Assembly::toFile(std::ostream& ofs) const
{
    for (unsigned int i = 0; i < sequences.size(); i++)
    {
        ofs << ">" << namePrefix << ".l" << overlap << ".c" << minCount << ".ctg." << i + 1 << '\n'
            << sequences[i] << '\n';
    }
}


bool assemblyCompare(const Assembly&  lhs, const Assembly& rhs) {
    if (lhs.longestContig != rhs.longestContig)
    {
        return lhs.longestContig > rhs.longestContig;
    }
    if (lhs.meanLength != rhs.meanLength)
    {
        return lhs.meanLength > rhs.meanLength;
    }
    if (lhs.numberOfContigs != rhs.numberOfContigs)
    {
        return lhs.numberOfContigs < rhs.numberOfContigs;
    }

    return lhs.minCount < rhs.minCount;
}


int assemble(char *readsFile, char *fastaOut, char* logfileOut, char *contigNamePrefix);


static PyObject * main_wrapper(PyObject * self, PyObject * args)
{
  char *readsFile;
  char *fastaOut;
  char *logOut;
  char *contigNamePrefix;
  int gotFromMain = 1;

  // parse arguments
  if (!PyArg_ParseTuple(args, "ssss", &readsFile, &fastaOut, &logOut, &contigNamePrefix)) {
      return NULL;
  }

  gotFromMain = assemble(readsFile, fastaOut, logOut, contigNamePrefix);
  return PyLong_FromLong((long) gotFromMain);
}


static PyMethodDef fermiliteMethods[] = {
   { "fermilite_ariba", main_wrapper, METH_VARARGS, "fermilite ariba" },
   { NULL, NULL, 0, NULL }
};


static struct PyModuleDef fermiliteModule = {
   PyModuleDef_HEAD_INIT,
   "fermilite_ariba",   /* name of module */
   NULL, /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   fermiliteMethods
};


PyMODINIT_FUNC
PyInit_fermilite_ariba(void)
{
    return PyModule_Create(&fermiliteModule);
}


int assemble(char *readsFile, char *fastaOut, char* logfileOut, char* contigNamePrefix)
{
    fml_opt_t opt;
    int n_seqs, n_utg;
    bseq1_t *seqs;
    fml_opt_init(&opt);
    opt.max_cnt = 10000;
    opt.mag_opt.flag |= MAG_F_AGGRESSIVE;
    std::vector<unsigned short> minCounts;
    minCounts.push_back(4);
    minCounts.push_back(17);
    minCounts.push_back(30);
    std::vector<unsigned short> overlaps;
    overlaps.push_back(6);
    overlaps.push_back(15);
    unsigned short assemblyCount = 0;
    std::ofstream ofs_stats(logfileOut);

    if (!ofs_stats.good())
    {
        std::cerr << "[ariba_fermilite] Error opening log output file '" << logfileOut << "'. Cannot continue" << std::endl;
        return 1;
    }

    ofs_stats << "Fermilite assembly stats:\n"
        << "Overlap\tMin_count\tContig_number\tMean_length\tLongest" << std::endl;

    std::ofstream ofs_fa(fastaOut);
    if (!ofs_fa.good())
    {
        ofs_stats.close();
        std::cerr << "[ariba_fermilite] Error opening fasta output file '" << fastaOut << "'. Cannot continue" << std::endl;
        return 1;
    }

    for (std::vector<unsigned short>::iterator overlapsIter = overlaps.begin(); overlapsIter != overlaps.end(); overlapsIter++)
    {
        for (std::vector<unsigned short>::iterator minCountIter = minCounts.begin(); minCountIter != minCounts.end(); minCountIter++)
        {
            opt.min_cnt = *minCountIter;
            opt.min_asm_ovlp = *overlapsIter;

            // need to get the reads from the file every time, instead of before
            // the loop because fml_assemble() destroys them :(
            seqs = bseq_read(readsFile, &n_seqs);
            if (seqs && n_seqs > 0) {
                fml_utg_t *utg;
                utg = fml_assemble(&opt, n_seqs, seqs, &n_utg);
                Assembly a(n_utg, utg, *minCountIter, *overlapsIter, contigNamePrefix);
                a.printStats(ofs_stats);
                a.toFile(ofs_fa);
                fml_utg_destroy(n_utg, utg);
                if (a.numberOfContigs > 0)
                {
                    assemblyCount++;
                }
            }
        }
    }

    ofs_stats.close();
    ofs_fa.close();

    if (assemblyCount == 0)
    {
        ofs_stats << "Didn't get any assemblies!\n";
        std::remove(fastaOut);
        return 1;
    }

    return 0;
}
