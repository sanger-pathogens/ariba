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
    Assembly(int numberOfUnitigs, fml_utg_t *unitigs, unsigned short minCountIn);
    void printStats(std::ostream& outStream) const;
    void toFile(std::string filename) const;

    uint32_t numberOfContigs;
    unsigned short minCount;
    uint32_t longestContig;
    float meanLength;
    std::vector<std::string> sequences;
};


Assembly::Assembly(int numberOfUnitigs, fml_utg_t *unitigs, unsigned short minCountIn)
{
    numberOfContigs = numberOfUnitigs;
    minCount = minCountIn;
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
    outStream << minCount << '\t' << numberOfContigs << '\t' << meanLength << '\t' << longestContig << std::endl;
}


void Assembly::toFile(std::string filename) const
{
    std::ofstream ofs(filename.c_str());
    if (!ofs.good())
    {
        std::cerr << "[ariba_fermilite] Error opening sequence output file '" << filename << "'. Cannot continue" << std::endl;
        exit(1);
    }

    for (unsigned int i = 0; i < sequences.size(); i++)
    {
        ofs << ">contig." << i << '\n'
            << sequences[i] << '\n';
    }

    ofs.close();
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


int assemble(char *readsFile, char *fastaOut, char* logfileOut);


static PyObject * main_wrapper(PyObject * self, PyObject * args)
{
  char *readsFile;
  char *fastaOut;
  char *logOut;
  int gotFromMain = 1;

  // parse arguments
  if (!PyArg_ParseTuple(args, "sss", &readsFile, &fastaOut, &logOut)) {
      return NULL;
  }

  gotFromMain = assemble(readsFile, fastaOut, logOut);
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


int assemble(char *readsFile, char *fastaOut, char* logfileOut)
{
    fml_opt_t opt;
    int n_seqs, n_utg;
    bseq1_t *seqs;
    fml_opt_init(&opt);
    opt.max_cnt = 10000;
    opt.min_asm_ovlp = 15;
    opt.mag_opt.flag |= MAG_F_AGGRESSIVE;
    std::vector<unsigned short> minCounts;
    minCounts.push_back(4);
    minCounts.push_back(8);
    minCounts.push_back(12);
    minCounts.push_back(16);
    minCounts.push_back(20);
    minCounts.push_back(25);
    minCounts.push_back(30);
    std::vector<Assembly> assemblies;
    std::ofstream ofs(logfileOut);

    if (!ofs.good())
    {
        std::cerr << "[ariba_fermilite] Error opening log output file '" << logfileOut << "'. Cannot continue" << std::endl;
        return 1;
    }

    ofs << "Fermilite assembly stats:\n"
        << "Min_count\tContig_number\tMean_length\tLongest" << std::endl;

    for (std::vector<unsigned short>::iterator minCountIter = minCounts.begin(); minCountIter != minCounts.end(); minCountIter++)
    {
        opt.min_cnt = *minCountIter;

        // need to get the reads from the file every time, instead of before
        // the loop because fml_assemble() destroys them :(
        seqs = bseq_read(readsFile, &n_seqs);
        if (seqs && n_seqs > 0) {
            fml_utg_t *utg;
            utg = fml_assemble(&opt, n_seqs, seqs, &n_utg);
            Assembly a(n_utg, utg, *minCountIter);
            assemblies.push_back(a);
            a.printStats(ofs);
            fml_utg_destroy(n_utg, utg);
        }
    }

    if (assemblies.size() == 0 || assemblies[0].numberOfContigs == 0)
    {
        ofs << "Didn't get any assemblies!\n";
        return 1;
    }

    std::sort(assemblies.begin(), assemblies.end(), &assemblyCompare);
    ofs << "Best assembly is from min_count " << assemblies[0].minCount << '\n';
    assemblies[0].toFile(fastaOut);
    ofs.close();
    return 0;
}
