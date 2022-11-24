import gffutils
import pandas as pd
import numpy as np
import tqdm
import json
import os

from multiprocessing import Process
from transform import gene2trans, trans2gene

PROCESS_NUMBER = 55
PATH_TO_RNA22_DB = "/home/dude/huge/dude/rbp-miRNA/rna22/rna22_energy_seed_pvalue.tsv"
PATH_TO_GTF_DB = "/home/dude/huge/dude/rbp-miRNA/genome/Homo_sapiens.GRCh38.96.db"
PATH_TO_OUTPUT_DIR = "/home/dude/huge/dude/rbp-miRNA/rna22/annotate/"


rna22_db = pd.read_csv(PATH_TO_RNA22_DB, header=None, sep="\t")
# gtf_db = gffutils.FeatureDB(PATH_TO_GTF_DB)

GENE_STRUCTURES = [
    "five_prime_utr",
    "three_prime_utr",
    "transcript",
    "exon",
    # "start_codon",
    # "stop_codon",
    "CDS",
    # "Selenocysteine",
    # "gene",
]

rna22_db = rna22_db[[0, 1, 2, 5, 6, 8, 9, 10, 15]]
rna22_db["miRNA"] = rna22_db[0].apply(lambda x: "-".join(x.split("_")))
rna22_db["gene_id"] = rna22_db[1].str.split("_").str[0]
rna22_db["transcript_id"] = rna22_db[1].str.split("_").str[1]
rna22_db["chr"] = rna22_db[1].str.split("_").str[2]

rna22_db["strand"] = rna22_db[1].str.split("_").str[3]
rna22_db.loc[rna22_db["strand"] == "-1", "strand"] = "-"
rna22_db.loc[rna22_db["strand"] == "1", "strand"] = "+"

rna22_db["transcriptomic_start"] = rna22_db[2].str.split("_").str[2]
rna22_db["transcriptomic_end"] = rna22_db[2].str.split("_").str[3]

rna22_db = rna22_db.drop(columns=[0, 1, 2])
rna22_db = rna22_db.rename(columns={
    5: "binding_score",
    6: "binding_seq", 
    8: "mRNA_binding_pattern", 
    9: "miRNA_binding_pattern",
    10: "heteroduplex_complementary_number",
    15: "pvalue"
})

def annotate(rna22_db, output_path):
    gtf_db = gffutils.FeatureDB(PATH_TO_GTF_DB)
    
    # output = open("/home/dude/huge/dude/rbp-miRNA/rna22/rna22_energy_seed_pvalue_annotated.tsv", "w")
    output = open(output_path, "w")
    annotation = [
        "miRNA",
        "chr",
        "strand",
        "genomic_start",
        "genomic_end",
        "binding_score",
        "binding_seq", 
        "mRNA_binding_pattern", 
        "miRNA_binding_pattern",
        "heteroduplex_complementary_number",
        "pvalue",
        "binding_gene_name",
        "binding_gene_id",
        "binding_transcript_id",
        "binding_structure_type",
        "binding_structure_start",
        "binding_structure_end",
        "transcriptomic_start",
        "transcriptomic_end"
    ]

    for column in annotation[:-1]:
        output.write(column + "\t")
    output.write(annotation[-1] + "\n")

    tqdrator = tqdm.tqdm(rna22_db.iterrows(), total=rna22_db.shape[0])
    for _, peak in tqdrator:
        chrom = peak['chr']
        strand = peak['strand']
        trans_start = peak['transcriptomic_start']
        trans_end = peak['transcriptomic_end']
        transcript_id = peak['transcript_id']

        start = trans2gene(str(chrom), str(strand),
                int(trans_start), gtf_db, transcript_id)
        end = trans2gene(str(chrom), str(strand),
                int(trans_end), gtf_db, transcript_id)

        region = f"{chrom}:{start}-{end}"
        for ftype in GENE_STRUCTURES:
            for gene_structure in gtf_db.region(region=region, featuretype=ftype, strand=strand):
                append = [
                    peak["miRNA"],
                    peak["chr"],
                    peak["strand"],
                    start,
                    end,
                    peak["binding_score"],
                    peak["binding_seq"], 
                    peak["mRNA_binding_pattern"], 
                    peak["miRNA_binding_pattern"],
                    peak["heteroduplex_complementary_number"],
                    peak["pvalue"],
                    gene_structure["gene_name"][0],
                    peak["gene_id"],
                    peak["transcript_id"],
                    ftype,
                    gene_structure.start,
                    gene_structure.end,
                    peak["transcriptomic_start"],
                    peak["transcriptomic_end"]
                ]

                for value in append[:-1]:
                    output.write(f"{value}\t")
                output.write(f"{append[-1]}\n")

    output.close()

def merge(output_path, dir_path):
    output = open(output_path, "w")
    header_flag = True
    for path in os.listdir(dir_path):
        inp = open(dir_path + path, "r")
        
        header = next(inp, None)
        if header_flag:
            output.write(header.rstrip("\n") + "\n")
            header_flag = False
        
        line = next(inp, None)
        while line:
            if line != "\n":
                output.write(line.rstrip("\n") + "\n")
            line = next(inp, None)

        inp.close()
        
    output.close()

if __name__ == "__main__":
    processes = []
    for i in range(PROCESS_NUMBER):
        start = len(rna22_db) // PROCESS_NUMBER * i
        end = len(rna22_db) // PROCESS_NUMBER * (i + 1)

        if i == PROCESS_NUMBER - 1:
            end = len(rna22_db)
        
        processes.append(Process(
            target=annotate,
            args=(
                rna22_db.iloc[start:end],
                f"{PATH_TO_OUTPUT_DIR}{i}.tsv"
            )
        ))

    for i in range(PROCESS_NUMBER):
        processes[i].start()

    for i in range(PROCESS_NUMBER):
        processes[i].join()
        print(f"Process {i} finished")

    merge("/home/dude/huge/dude/rbp-miRNA/rna22/rna22_energy_seed_pvalue_annotated.tsv",
        PATH_TO_OUTPUT_DIR)
