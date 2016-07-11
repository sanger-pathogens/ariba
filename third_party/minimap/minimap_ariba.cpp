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
KSEQ_INIT(gzFile, gzread)



typedef std::vector<std::pair<uint32_t, bool> > MapPositionVector;

void loadClusters(std::string& filename, std::map<std::string, std::string>& refnameToCluster);
void chooseCluster(std::string outfile, std::map<std::string, uint64_t>& refnameToScore, std::map<std::string, std::string>& refnameToCluster);
void writeClusterCountsFile(std::string outfile, const std::map<std::string, uint64_t>& readCounters, const std::map<std::string, uint64_t>& baseCounters);
void writeInsertHistogramFile(std::string outfile, const std::map<uint32_t, uint32_t>& insertHist);
void writeProperPairsFile(std::string outfile, uint32_t properPairs);


int main(int argc, char *argv[])
{
    if (argc < 6) {
        std::cerr << "Usage: minimap_ariba <clusters> <target.fa> <query_fwd> <query_rev> <outprefix>\n"
                  << "This is NOT intended to be run directly! Should be called during ARIBA pipeline\n";
        return 1;
    }

    mm_verbose = 0;
    std::map<std::string, uint64_t> refnameToScore;
    std::map<std::string, std::string> refnameToCluster;
    std::string clustersFile(argv[1]);
    loadClusters(clustersFile, refnameToCluster);
    std::map<std::string, uint64_t> readCounters;
    std::map<std::string, uint64_t> baseCounters;
    std::map<uint32_t, uint32_t> insertHist;
    uint32_t properPairs = 0;
    std::string outprefix(argv[5]);
    std::string readsOutfile = outprefix + ".reads";
    std::ofstream ofs;
    ofs.open(readsOutfile.c_str());
    if (!ofs.good())
    {
        std::cerr << "Error opening reads output file '" << readsOutfile << "'. Cannot continue" << std::endl;
        return 1;
    }

    // open query file for reading; you may use your favorite FASTA/Q parser
    gzFile infile1 = gzopen(argv[3], "r");
    assert(infile1);
    gzFile infile2 = gzopen(argv[4], "r");
    assert(infile2);
    kseq_t *ks1 = kseq_init(infile1);
    kseq_t *ks2 = kseq_init(infile2);

    // create index for target; we are creating one index for all target sequence
    int n_threads = 1;
    int w = 10, k = 15;
    mm_idx_t *mi = mm_idx_build(argv[2], w, k, n_threads);
    assert(mi);

    // mapping
    mm_mapopt_t opt;
    mm_mapopt_init(&opt); // initialize mapping parameters
    mm_tbuf_t *tbuf1 = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
    mm_tbuf_t *tbuf2 = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
    while (kseq_read(ks1) >= 0) { // each kseq_read() call reads one query sequence
        assert(kseq_read(ks2) >= 0);
        const mm_reg1_t *reg1, *reg2;
        int j, n_reg1, n_reg2;
        std::map<std::string, MapPositionVector> positions1;
        std::map<std::string, MapPositionVector> positions2;

        // get all hits for the forward and reverse reads
        reg1 = mm_map(mi, ks1->seq.l, ks1->seq.s, &n_reg1, tbuf1, &opt, 0);
        reg2 = mm_map(mi, ks2->seq.l, ks2->seq.s, &n_reg2, tbuf2, &opt, 0);
        if (n_reg1 > 0 || n_reg2 > 0)
        {
            std::set<std::string> refnames;
            for (j  =0; j < n_reg1; ++j)
            {
                const mm_reg1_t *r = &reg1[j];
                refnames.insert(mi->name[r->rid]);
                refnameToScore[mi->name[r->rid]] += r->cnt;
                uint32_t coord = r->rev ? std::max(r->rs, r->re) : std::min(r->rs, r->re);
                positions1[mi->name[r->rid]].push_back(std::make_pair(coord, r->rev));
            }
            for (j  =0; j < n_reg2; ++j)
            {
                const mm_reg1_t *r = &reg2[j];
                refnames.insert(mi->name[r->rid]);
                refnameToScore[mi->name[r->rid]] += r->cnt;
                uint32_t coord = r->rev ? std::max(r->rs, r->re) : std::min(r->rs, r->re);
                positions2[mi->name[r->rid]].push_back(std::make_pair(coord, r->rev));
            }

            bool foundProperPair = false;

            for (std::set<std::string>::const_iterator iter = refnames.begin(); iter != refnames.end(); iter++)
            {
                std::string cluster = refnameToCluster[*iter];
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

                    insertHist[insertSize]++;
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
        assert(getline(ss, cluster, '\t'));
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
        assert ( baseCounters.find(iter->first) != baseCounters.end());
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
