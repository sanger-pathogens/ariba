#include <iostream>
#include <set>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include <map>
#include <vector>
#include <algorithm>
#include "minimap.h"
#include "kseq.h"
#include "Python.h"
KSEQ_INIT(gzFile, gzread)


typedef std::vector<std::pair<uint32_t, bool> > MapPositionVector;

void loadClusters(std::string& filename, std::map<std::string, std::string>& refnameToCluster);
void chooseCluster(std::string outfile, std::map<std::string, uint64_t>& refnameToScore, std::map<std::string, std::string>& refnameToCluster);
void writeClusterCountsFile(std::string outfile, const std::map<std::string, uint64_t>& readCounters, const std::map<std::string, uint64_t>& baseCounters);
void writeInsertHistogramFile(std::string outfile, const std::map<uint32_t, uint32_t>& insertHist);
void writeProperPairsFile(std::string outfile, uint32_t properPairs);
bool readMappingOk(const mm_reg1_t* r, const mm_idx_t* mi, const kseq_t *ks1, uint32_t endTolerance);

int run_minimap(char *clustersFileIn, char *refFileIn, char *readsFile1In, char *readsFile2In, char *outprefixIn);

static PyObject * main_wrapper(PyObject * self, PyObject * args)
{
  char *clustersFile;
  char *refFile;
  char *readsFile1;
  char *readsFile2;
  char *outprefix;
  int gotFromMain = 1;

  // parse arguments
  if (!PyArg_ParseTuple(args, "sssss", &clustersFile, &refFile, &readsFile1, &readsFile2, &outprefix)) {
      return NULL;
  }

  gotFromMain = run_minimap(clustersFile, refFile, readsFile1, readsFile2, outprefix);
  return PyLong_FromLong((long) gotFromMain);
}


static PyMethodDef minimapMethods[] = {
   { "minimap_ariba", main_wrapper, METH_VARARGS, "minimap ariba" },
   { NULL, NULL, 0, NULL }
};


static struct PyModuleDef minimapModule = {
   PyModuleDef_HEAD_INIT,
   "minimap_ariba",   /* name of module */
   NULL, /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   minimapMethods
};

PyMODINIT_FUNC
PyInit_minimap_ariba(void)
{
    return PyModule_Create(&minimapModule);
}



int run_minimap(char *clustersFileIn, char *refFileIn, char *readsFile1In, char *readsFile2In, char *outprefixIn)
{
    mm_verbose = 0;
    std::map<std::string, uint64_t> refnameToScore;
    std::map<std::string, std::string> refnameToCluster;
    std::string clustersFile(clustersFileIn);
    loadClusters(clustersFile, refnameToCluster);
    std::map<std::string, uint64_t> readCounters;
    std::map<std::string, uint64_t> baseCounters;
    std::map<uint32_t, uint32_t> insertHist;
    uint32_t properPairs = 0;
    std::string outprefix(outprefixIn);
    std::string readsOutfile = outprefix + ".reads";
    std::ofstream ofs;
    ofs.open(readsOutfile.c_str());
    if (!ofs.good())
    {
        std::cerr << "[ariba_minimap] Error opening reads output file '" << readsOutfile << "'. Cannot continue" << std::endl;
        return 1;
    }

    // open query file for reading; you may use your favorite FASTA/Q parser
    gzFile infile1 = gzopen(readsFile1In, "r");
    if (!infile1)
    {
        std::cerr << "[ariba_minimap] Error opening file " << readsFile1In << std::endl;
        return 1;
    }
    gzFile infile2 = gzopen(readsFile2In, "r");
    if (!infile2)
    {
        std::cerr << "[ariba_minimap] Error opening file " << readsFile2In << std::endl;
        return 1;
    }
    kseq_t *ks1 = kseq_init(infile1);
    kseq_t *ks2 = kseq_init(infile2);

    // create index for target; we are creating one index for all target sequence
    int n_threads = 1;
    int w = 10, k = 15;
    mm_idx_t *mi = mm_idx_build(refFileIn, w, k, n_threads);
    if (!mi)
    {
        std::cerr << "[ariba_minimap] Error indexing" << std::endl;
        return 1;
    }

    // mapping
    mm_mapopt_t opt;
    mm_mapopt_init(&opt); // initialize mapping parameters
    mm_tbuf_t *tbuf1 = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
    mm_tbuf_t *tbuf2 = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
    while (kseq_read(ks1) >= 0) { // each kseq_read() call reads one query sequence
        const mm_reg1_t *reg1, *reg2;
        int j, n_reg1, n_reg2;
        int ks2ok = kseq_read(ks2);
        if (ks2ok == 0)
        {
            std::cerr << "Error getting mate of read. Cannot continue" << std::endl;
            return 1;
        }

        // get all hits for the forward and reverse reads
        reg1 = mm_map(mi, ks1->seq.l, ks1->seq.s, &n_reg1, tbuf1, &opt, 0);
        reg2 = mm_map(mi, ks2->seq.l, ks2->seq.s, &n_reg2, tbuf2, &opt, 0);

        if (n_reg1 > 0 || n_reg2 > 0)
        {
            std::map<std::string, MapPositionVector> positions1;
            std::map<std::string, MapPositionVector> positions2;
            std::set<std::string> refnames;

            for (j  =0; j < n_reg1; ++j)
            {
                const mm_reg1_t *r = &reg1[j];
                if (readMappingOk(r, mi, ks1, (int) 1.1 * k))
                {
                    refnames.insert(mi->name[r->rid]);
                    refnameToScore[mi->name[r->rid]] += r->cnt;
                    uint32_t coord = r->rev ? std::max(r->rs, r->re) : std::min(r->rs, r->re);
                    positions1[mi->name[r->rid]].push_back(std::make_pair(coord, r->rev));
                }
            }
            for (j  =0; j < n_reg2; ++j)
            {
                const mm_reg1_t *r = &reg2[j];
                if (readMappingOk(r, mi, ks2, (int) 1.1 * k))
                {
                    refnames.insert(mi->name[r->rid]);
                    refnameToScore[mi->name[r->rid]] += r->cnt;
                    uint32_t coord = r->rev ? std::max(r->rs, r->re) : std::min(r->rs, r->re);
                    positions2[mi->name[r->rid]].push_back(std::make_pair(coord, r->rev));
                }
            }

            bool foundProperPair = false;
            std::set<std::string> usedClusters;
            for (std::set<std::string>::const_iterator iter = refnames.begin(); iter != refnames.end(); iter++)
            {
                std::string cluster = refnameToCluster[*iter];

                // do not write a read pair to the same cluster more than once
                if (usedClusters.find(cluster) != usedClusters.end())
                {
                    continue;
                }

                usedClusters.insert(cluster);
                readCounters[cluster]++;
                ofs << cluster << '\t' << readCounters[cluster] << '\t' << ks1->seq.s << '\t' << ks1->qual.s << '\n';
                readCounters[cluster]++;
                ofs << cluster << '\t' << readCounters[cluster] << '\t' << ks2->seq.s << '\t' << ks2->qual.s << '\n';
                baseCounters[cluster] += ks1->seq.l + ks2->seq.l;

                // get insert size info, if reads mapped as proper pair
                if (positions1.find(*iter) != positions1.end() && positions2.find(*iter) != positions2.end())
                {
                    if (positions1[*iter].size() != 1 || positions2[*iter].size() != 1 || positions1[*iter][0].second == positions2[*iter][0].second)
                    {
                        continue;
                    }

                    uint32_t insertSize;

                    if (positions1[*iter][0].second && positions1[*iter][0].first > positions2[*iter][0].first)
                    {
                        insertSize = positions1[*iter][0].first - positions2[*iter][0].first + 1;
                    }
                    else if (positions2[*iter][0].second && positions2[*iter][0].first > positions1[*iter][0].first)
                    {
                        insertSize = positions2[*iter][0].first - positions1[*iter][0].first + 1;
                    }
                    else
                    {
                        continue;
                    }

                    insertHist[insertSize + 2*k]++;
                    foundProperPair = true;
                }
            }

            if (foundProperPair)
            {
                properPairs++;
            }
        }
    }
    mm_tbuf_destroy(tbuf1);
    mm_tbuf_destroy(tbuf2);

    // deallocate index and close the query file
    mm_idx_destroy(mi);
    kseq_destroy(ks1);
    kseq_destroy(ks2);
    gzclose(infile1);
    gzclose(infile2);
    ofs.close();
    chooseCluster(outprefix + ".cluster2representative", refnameToScore, refnameToCluster);
    writeClusterCountsFile(outprefix + ".clusterCounts", readCounters, baseCounters);
    writeInsertHistogramFile(outprefix + ".insertHistogram", insertHist);
    writeProperPairsFile(outprefix + ".properPairs", properPairs);
    return 0;
}


void loadClusters(std::string& filename, std::map<std::string, std::string>& refnameToCluster)
{
    std::ifstream ifs;
    std::string line;
    ifs.open(filename.c_str());
    if (!ifs.good())
    {
        std::cerr << "Error opening clusters file '" << filename << "'. Cannot continue" << std::endl;
        exit(1);
    }

    while(getline(ifs, line))
    {
        std::stringstream ss(line);
        std::string cluster;
        std::string seqname;
        getline(ss, cluster, '\t');
        if (cluster.size() == 0)
        {
            std::cerr << "Error reading clusters file at the following line\n" << line << std::endl;
            exit(1);
        }

        while (getline(ss, seqname, '\t'))
        {
            refnameToCluster[seqname] = cluster;
        }
    }

    ifs.close();
}


void chooseCluster(std::string outfile, std::map<std::string, uint64_t>& refnameToScore, std::map<std::string, std::string>& refnameToCluster)
{
    std::map<std::string, uint64_t> bestClusterScore;
    std::map<std::string, std::string> bestCluster;
    std::map<std::string, uint64_t>::iterator iter;
    for (iter = refnameToScore.begin(); iter != refnameToScore.end(); iter++)
    {
        std::string cluster = refnameToCluster[iter->first];
        if (bestClusterScore.find(cluster) == bestClusterScore.end() || bestClusterScore[cluster] < iter->second)
        {
            bestClusterScore[cluster] = iter->second;
            bestCluster[cluster] = iter->first;
        }
    }

    std::ofstream ofs;
    ofs.open(outfile.c_str());
    if (!ofs.good())
    {
        std::cerr << "Error opening output best cluster file '" << outfile << "'. Cannot continue" << std::endl;
        exit(1);
    }

    for (iter = bestClusterScore.begin(); iter != bestClusterScore.end(); iter++)
    {
        ofs << iter->first << '\t' << bestCluster[iter->first] << '\n';
    }

    ofs.close();
}


void writeClusterCountsFile(std::string outfile, const std::map<std::string, uint64_t>& readCounters, const std::map<std::string, uint64_t>& baseCounters)
{
    std::ofstream ofs;
    ofs.open(outfile.c_str());
    if (!ofs.good())
    {
        std::cerr << "Error opening output cluster reads/bases counts file '" << outfile << "'. Cannot continue" << std::endl;
        exit(1);
    }

    for (std::map<std::string, uint64_t>::const_iterator iter = readCounters.begin(); iter != readCounters.end(); iter++)
    {
        if (baseCounters.find(iter->first) == baseCounters.end())
        {
            std::cerr << "Error writing cluster counts file" << std::endl;
            exit(1);
        }
        ofs << iter->first << '\t' << iter->second << '\t' << baseCounters.find(iter->first)->second << '\n';
    }

    ofs.close();
}


void writeInsertHistogramFile(std::string outfile, const std::map<uint32_t, uint32_t>& insertHist)
{
    std::ofstream ofs;
    ofs.open(outfile.c_str());
    if (!ofs.good())
    {
        std::cerr << "Error opening output insert histogram file '" << outfile << "'. Cannot continue" << std::endl;
        exit(1);
    }

    for (std::map<uint32_t, uint32_t>::const_iterator iter = insertHist.begin(); iter != insertHist.end(); iter++)
    {
        ofs << iter->first << '\t' << iter->second << '\n';
    }

    ofs.close();
}


void writeProperPairsFile(std::string outfile, uint32_t properPairs)
{
    std::ofstream ofs;
    ofs.open(outfile.c_str());
    if (!ofs.good())
    {
        std::cerr << "Error opening output proper pairs count file '" << outfile << "'. Cannot continue" << std::endl;
        exit(1);
    }

    ofs << properPairs << '\n';
    ofs.close();
}


bool readMappingOk(const mm_reg1_t* r, const mm_idx_t* mi, const kseq_t *ks, uint32_t endTolerance)
{
    // coords are same style as python (0-based, end is one past the end)
    assert (r->qs < r->qe && r->rs <  r->re);

    if (r->qe - r->qs < std::min((unsigned) 50, (int) 0.5 * ks->seq.l))
    {
        return false;
    }

    uint32_t refLength = mi->len[r->rid];
    bool startOk;
    bool endOk;
    if (r->rev)
    {
        startOk = (r->qs < endTolerance || refLength - r->re < endTolerance);
        endOk = (ks->seq.l - r->qe < endTolerance || r->rs < endTolerance);
    }
    else
    {
        startOk = (r->qs < endTolerance || r->rs < endTolerance);
        endOk = (ks->seq.l - r->qe < endTolerance || refLength - r->re < endTolerance);
    }

    return (startOk && endOk);
}
