import gffutils
import pandas as pd
import numpy as np
import tqdm
import json
import os
import copy

from multiprocessing import Process
from transform import gene2trans, trans2gene

PROCESS_NUMBER = 55
PATH_TO_CLIP_DB = "/home/dude/huge/dude/rbp-miRNA/clip/clip_db.tsv"
PATH_TO_GTF_DB = "/home/dude/huge/dude/rbp-miRNA/genome/Homo_sapiens.GRCh38.96.db"
PATH_TO_OUTPUT_DIR = "/home/dude/huge/dude/rbp-miRNA/clip/annotate/"


clip_db = pd.read_csv(PATH_TO_CLIP_DB, sep="\t")
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

def annotate(clip_db, output_path):
    gtf_db = gffutils.FeatureDB(PATH_TO_GTF_DB)
    
    # output = open("/home/dude/huge/dude/rbp-miRNA/clip/clip_db_annotated.tsv", "w")
    output = open(output_path, "w")
    for column in clip_db.columns:
        if column == "start":
            column = "genomic_start"
        if column == "end":
            column = "genomic_end"
        output.write(column + "\t")

    annotation = [
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

    counter = 0
    tqdrator = tqdm.tqdm(clip_db.iterrows(), total=clip_db.shape[0])
    for _, peak in tqdrator:
        chrom = peak['chr']
        strand = peak['strand']
        start, end = peak['start'], peak['end']
        region = f"{chrom}:{start}-{end}"

        for ftype in GENE_STRUCTURES:
            for gene_structure in gtf_db.region(region=region, featuretype=ftype, strand=strand):
                # binding region inside gene_structure
                # if (start < gene_structure.start) or (stop > gene_structure.stop):
                #     continue
                append = list(peak.values)
                append.extend([
                    gene_structure["gene_name"][0],
                    gene_structure["gene_id"][0],
                    gene_structure["transcript_id"][0],
                    ftype,
                    gene_structure.start,
                    gene_structure.stop
                ])

                transcript_id = gene_structure["transcript_id"][0]
                trans_start = gene2trans(str(chrom), str(strand),
                        int(start), gtf_db, transcript_id)
                trans_end = gene2trans(str(chrom), str(strand),
                        int(end), gtf_db, transcript_id)

                append.extend([trans_start, trans_end])

                for value in append[:-1]:
                    output.write(f"{value}\t")
                output.write(f"{append[-1]}\n")

                counter += 1

        tqdrator.set_description(f"{counter}")

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
        start = len(clip_db) // PROCESS_NUMBER * i
        end = len(clip_db) // PROCESS_NUMBER * (i + 1)

        if i == PROCESS_NUMBER - 1:
            end = len(clip_db)
        
        processes.append(Process(
            target=annotate,
            args=(
                clip_db.iloc[start:end],
                f"{PATH_TO_OUTPUT_DIR}{i}.tsv"
            )
        ))

    for i in range(PROCESS_NUMBER):
        processes[i].start()

    for i in range(PROCESS_NUMBER):
        processes[i].join()
        print(f"Process {i} finished")

    merge("/home/dude/huge/dude/rbp-miRNA/clip/clip_db_annotated.tsv",
        PATH_TO_OUTPUT_DIR)


