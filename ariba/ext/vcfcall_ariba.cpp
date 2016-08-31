#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <algorithm>
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "Python.h"


void openFileWrite(std::string filename, std::ofstream& ofs);
void split(const std::string& s, char delim, std::vector<std::string>& elems);
void split(const std::string& s, char delim, std::vector<uint32_t> &elems);
std::string getKey(const std::string& key, const std::string& s);
bool adStringPassesFilter(std::string& adString, uint32_t minTotalDepth, uint32_t minSecondDepth, float maxAlleleFreq);

int run(char* infileIn, char* outprefixIn, uint32_t minTotalDepth, uint32_t minSecondDepth, float maxAlleleFreq);


static PyObject * main_wrapper(PyObject * self, PyObject * args)
{
  char *infile;
  char *outprefix;
  uint32_t minSecondDepth;
  uint32_t minTotalDepth;
  float maxAlleleFreq;
  int gotFromMain = 1;

  // parse arguments
  if (!PyArg_ParseTuple(args, "ssiif", &infile, &outprefix, &minTotalDepth, &minSecondDepth, &maxAlleleFreq)) {
      return NULL;
  }

  gotFromMain = run(infile, outprefix, minTotalDepth, minSecondDepth, maxAlleleFreq);
  return PyLong_FromLong((long) gotFromMain);
}


static PyMethodDef vcfcallMethods[] = {
   { "vcfcall_ariba", main_wrapper, METH_VARARGS, "vcfcall ariba" },
   { NULL, NULL, 0, NULL }
};


static struct PyModuleDef vcfcallModule = {
   PyModuleDef_HEAD_INIT,
   "vcfcall_ariba",   /* name of module */
   NULL, /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   vcfcallMethods
};

PyMODINIT_FUNC
PyInit_vcfcall_ariba(void)
{
    return PyModule_Create(&vcfcallModule);
}

void openFileWrite(std::string filename, std::ofstream& ofs)
{
    ofs.open(filename.c_str());
    if (!ofs.good())
    {
        std::cerr << "[ariba vcfcall] Error opening output file '" << filename << "'. Cannot continue" << std::endl;
        exit(1);
    }
}


void split(const std::string& s, char delim, std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    elems.clear();

    while(getline(ss, item, delim))
    {
        elems.push_back(item);
    }
}


void split(const std::string& s, char delim, std::vector<uint32_t> &elems)
{
    std::stringstream ss(s);
    std::string item;
    elems.clear();

    while(getline(ss, item, delim))
    {
        elems.push_back(atof(item.c_str()));
    }
}

std::string getKey(const std::string& key, const std::string& s)
{
    std::size_t found = s.find(key);
    if (found == std::string::npos)
    {
        return "";
    }

    std::size_t foundSemicolon = s.find(";", found);

    if (foundSemicolon == std::string::npos || found + key.size() >= foundSemicolon)
    {
        return "";
    }

    return s.substr(found + key.size(), foundSemicolon - (found + key.size()));
}


bool adStringPassesFilter(std::string& adString, uint32_t minTotalDepth, uint32_t minSecondDepth, float maxAlleleFreq)
{
    std::vector<uint32_t> depths;
    split(adString, ',', depths);
    uint32_t firstDepth = depths[0];
    std::sort(depths.begin(), depths.end());
    uint32_t totalDepth = 0;

    for (std::vector<uint32_t>::const_iterator p = depths.begin(); p != depths.end(); p++)
    {
        totalDepth += *p;
    }

    bool secondDepthOk = (depths.size() == 1 || (depths.size() > 1 && depths[depths.size() - 2] >= minSecondDepth));
    bool maxDepthOk = (totalDepth >= minTotalDepth && 1.0 * depths.back() / totalDepth <= maxAlleleFreq);
    return ( firstDepth < depths.back() || (secondDepthOk && maxDepthOk) );
}



int run(char* infileIn, char* outprefixIn, uint32_t minTotalDepth, uint32_t minSecondDepth, float maxAlleleFreq)
{
    std::string infile(infileIn);
    std::string outprefix(outprefixIn);
    std::string readDepthOut = outprefix + ".read_depths";
    std::string variantOut = outprefix + ".vcf";
    std::string contigDepthOut = outprefix + ".contig_depths";
    std::ofstream readDepthOutFh;
    std::ofstream variantOutFh;
    std::ofstream contigDepthOutFh;
    openFileWrite(readDepthOut, readDepthOutFh);
    openFileWrite(variantOut, variantOutFh);
    std::ifstream ifs;
    ifs.open(infile.c_str());
    std::string line;
    std::vector<std::string> fields;
    std::vector<std::string> varNucs;
    std::string refString;
    std::string varString;
    std::string adString;
    std::string dpString;
    std::map<std::string, uint32_t> totalDepths;

    if (!ifs.good())
    {
        std::cerr << "[ariba vcfcall] Error opening vcf file '" << infile << "'. Cannot continue" << std::endl;
        exit(1);
    }

    while(getline(ifs, line))
    {
        if (line[0] == '#')
        {
            continue;
        }

        split(line, '\t', fields);
        adString = getKey("AD=", fields[7]);
        dpString = getKey("DP=", fields[7]);

        if (adString.size() == 0 || dpString.size() == 0)
        {
            std::cerr << "Error getting AD=... or DP=... from vcf line. Cannot continue. Line was:\n" << line << std::endl;
            return 1;
        }

        if (fields[7].find("INDEL") != std::string::npos)
        {
            refString = fields[3];
            varString = fields[4];
        }
        else
        {
            if (adString.substr(adString.size() - 2).compare(",0") == 0)
            {
                adString.resize(adString.size() - 2);
            }

            if (fields[4].compare("<*>"))
            {
                split(fields[4], ',', varNucs);
                varString = "";
                for (std::vector<std::string>::const_iterator p = varNucs.begin(); p != varNucs.end(); p++)
                {
                    if (p->compare("<*>"))
                    {
                        varString += *p + ",";
                    }
                }
                varString.resize(varString.size() - 1);
            }
            else
            {
                varString = ".";
            }

            uint32_t depth = atof(dpString.c_str());
            totalDepths[fields[0]] += depth;
        }

        readDepthOutFh << fields[0]
                       << '\t' << fields[1]
                       << '\t' << fields[3]
                       << '\t' << varString
                       << '\t' << dpString
                       << '\t' << adString << '\n';

        if (varString.compare(".") && adStringPassesFilter(adString, minTotalDepth, minSecondDepth, maxAlleleFreq))
        {
            variantOutFh << line << '\n';
        }
    }

    ifs.close();
    readDepthOutFh.close();
    variantOutFh.close();

    openFileWrite(contigDepthOut, contigDepthOutFh);
    for (std::map<std::string, uint32_t>::const_iterator p = totalDepths.begin(); p != totalDepths.end(); p++)
    {
        contigDepthOutFh << p->first << '\t' << p->second << '\n';
    }
    contigDepthOutFh.close();

    return 0;
}

