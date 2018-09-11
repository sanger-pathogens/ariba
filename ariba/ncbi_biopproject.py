BIOPROJECT = "PRJNA313047"
#download the database as genbank
cmd =f"esearch -db nucleotide -query '{BIOPROJECT}' | efetch -db nuccore -format gbwithparts > {BIOPROJECT}.gbk"
from pathlib import Path
import os
if not Path(f"{BIOPROJECT}.gbk").exists():
    os.system(cmd)
#pull out the records as fasta
from Bio.Alphabet import generic_dna, generic_protein
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

with open(f"{BIOPROJECT}.gbk", "r") as input_handle:
    for gb_record in SeqIO.parse(input_handle, "genbank"):
        n=0
        record_new=[]
        for index, feature in enumerate(gb_record.features):
            if feature.type == 'CDS':
                n+=1
                gb_feature = gb_record.features[index]
                id = None
                try:
                    id = gb_feature.qualifiers["allele"]
                except:
                    try:
                        try:
                            id = gb_feature.qualifiers["gene"]
                        except:
                            id = gb_feature.qualifiers["locus_tag"]
                    except KeyError:
                        import sys
                        sys.stderr(f"gb_feature.qualifer not found")
                        
                accession = gb_record.id
                product = gb_feature.qualifiers['product']
                locus_tag = gb_feature.qualifiers['locus_tag']
                desc = '| '
                seq_out = Seq(str(gb_feature.extract(gb_record.seq)), generic_dna)
                record_new.append(SeqRecord(seq_out,
                             id=f"{str(id[0])}.{str(accession)}",
                             description=""))
        if len(record_new) == 1:
                print(record_new[0].format("fasta").rstrip())
