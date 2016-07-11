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
#include "minimap.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)


void loadClusters(std::string& filename, std::map<std::string, std::string>& refnameToCluster);

void chooseCluster(std::string outfile, std::map<std::string, uint64_t>& refnameToScore, std::map<std::string, std::string>& refnameToCluster);


int main(int argc, char *argv[])
{
    if (argc < 6) {
        std::cerr << "Usage: minimap_ariba <clusters> <target.fa> <query_fwd> <query_rev> <outprefix>\n"
                  << "This is NOT intended to be run directly! Should be called during ARIBA pipeline\n";
        return 1;
    }

    std::map<std::string, uint64_t> refnameToScore;
    std::map<std::string, std::string> refnameToCluster;
    std::string clustersFile(argv[1]);
    loadClusters(clustersFile, refnameToCluster);
    std::map<std::string, uint64_t> readCounters;
    std::string outprefix(argv[5]);
    std::string readsOutfile = outprefix + ".reads";
    std::ofstream ofs;
    ofs.open(readsOutfile.c_str());
    if (!ofs.good())
    {
        std::cerr << "Error opening reads output file '" << readsOutfile << "'. Cannot continue" << std::endl;
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
    mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
    while (kseq_read(ks1) >= 0) { // each kseq_read() call reads one query sequence
        assert(kseq_read(ks2) >= 0);
        const mm_reg1_t *reg1, *reg2;
        int j, n_reg1, n_reg2;
        // get all hits for the forward and reverse reads
        reg1 = mm_map(mi, ks1->seq.l, ks1->seq.s, &n_reg1, tbuf, &opt, 0);
        reg2 = mm_map(mi, ks2->seq.l, ks2->seq.s, &n_reg2, tbuf, &opt, 0);
        if (n_reg1 > 0 || n_reg2 > 0)
        {
            std::set<std::string> refnames;
            for (j  =0; j < n_reg1; ++j)
            {
                const mm_reg1_t *r = &reg1[j];
                refnames.insert(mi->name[r->rid]);
                refnameToScore[mi->name[r->rid]] += r->cnt;
            }
            for (j  =0; j < n_reg2; ++j)
            {
                const mm_reg1_t *r = &reg2[j];
                refnames.insert(mi->name[r->rid]);
                refnameToScore[mi->name[r->rid]] += r->cnt;
            }

            for (std::set<std::string>::const_iterator iter = refnames.begin(); iter != refnames.end(); iter++)
            {
                std::string cluster = refnameToCluster[*iter];
                readCounters[cluster]++;
                ofs << cluster << '\t' << readCounters[cluster] << '\t' << ks1->seq.s << '\t' << ks1->qual.s << '\n';
                readCounters[*iter]++;
                ofs << cluster << '\t' << readCounters[cluster] << '\t' << ks2->seq.s << '\t' << ks2->qual.s << '\n';
            }

        }
    }
    mm_tbuf_destroy(tbuf);

    // deallocate index and close the query file
    mm_idx_destroy(mi);
    kseq_destroy(ks1);
    kseq_destroy(ks2);
    gzclose(infile1);
    gzclose(infile2);
    ofs.close();
    chooseCluster(outprefix + ".clusters", refnameToScore, refnameToCluster);
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
    }

    for (iter = bestClusterScore.begin(); iter != bestClusterScore.end(); iter++)
    {
        ofs << iter->first << '\t' << bestCluster[iter->first] << '\n';
    }

    ofs.close();
}
