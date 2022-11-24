import pandas as pd
import numpy as np
import json
import gffutils
import tqdm
import os

from multiprocessing import Process


PROCESS_NUMBER = 55
PATH_TO_CLIP_DB = "/home/dude/huge/dude/rbp-miRNA/clip/clip_db_annotated.tsv"
PATH_TO_RNA22_DB = "/home/dude/huge/dude/rbp-miRNA/rna22/rna22_energy_seed_pvalue_annotated.tsv"
PATH_TO_OUTPUT_DIR = "/home/dude/huge/dude/rbp-miRNA/results/merge/"

GENE_STRUCTURES = [
    "five_prime_utr",
    "three_prime_utr",
    "CDS"
]

CONDITION = [
    "clip_transcriptomic_start",
    "clip_transcriptomic_end",
    "rna22_transcriptomic_start",
    "rna22_transcriptomic_end",
    "clip_genomic_start",
    "clip_genomic_end",
    "rna22_genomic_start",
    "rna22_genomic_end"
]

MUTUAL_COLUMNS = [
    "chr",
    "strand",
    "binding_gene_name",
    "binding_gene_id",
    "binding_transcript_id",
    # "binding_structure_type",
    # "binding_structure_start",
    # "binding_structure_end",
]


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

def cross(rna22_groups, clip_db, output_path):
    clip_groups = clip_db.groupby(MUTUAL_COLUMNS)
    clip_keys = set(clip_groups.groups.keys())

    output = open(output_path, "w")
    flag = True

    success_number = 0
    tqdrator = tqdm.tqdm(rna22_groups, total=len(rna22_groups))
    for group_name, rna22_group in tqdrator:
        if not group_name in clip_keys:
            continue
        
        success_number += 1
        clip_group = clip_groups.get_group(group_name)
        clip_group = clip_group.drop(columns=MUTUAL_COLUMNS)

        append = rna22_group.merge(clip_group, how="cross")
        conditions = [
            (append[cond] >= 0) for cond in CONDITION
        ]
        
        condition = conditions[0]
        for cond in conditions[1:]:
            condition = condition & cond

        append = append.loc[condition]

        tqdrator.set_description(f"{success_number}")

        if flag:
            for column in append.columns[:-1]:
                output.write(f"{column}\t")
            output.write(f"{append.columns[-1]}\n")
            flag = False

        for _, row in append.iterrows():
            for value in row.values[:-1]:
                output.write(f"{value}\t")
            output.write(f"{row.values[-1]}\n")


if __name__ == "__main__":
    rna22_db = pd.read_csv(PATH_TO_RNA22_DB, sep="\t")
    rna22_db = rna22_db.loc[rna22_db["binding_structure_type"].isin(GENE_STRUCTURES)]
    
    clip_db = pd.read_csv(PATH_TO_CLIP_DB, sep="\t")
    clip_db = clip_db.loc[clip_db["binding_structure_type"].isin(GENE_STRUCTURES)]

    rna22_columns = MUTUAL_COLUMNS + [column for column in rna22_db.columns if not column in MUTUAL_COLUMNS]
    clip_columns = MUTUAL_COLUMNS + [column for column in clip_db.columns if not column in MUTUAL_COLUMNS]

    rna22_db = rna22_db[rna22_columns]
    rna22_columns = dict(zip(rna22_columns,
        MUTUAL_COLUMNS + \
        ["rna22_" + column for column in rna22_columns if not column in MUTUAL_COLUMNS]
    ))
    rna22_db.rename(columns=rna22_columns, inplace=True)
    rna22_groups = rna22_db.groupby(MUTUAL_COLUMNS)

    clip_db = clip_db[clip_columns]
    clip_columns = dict(zip(clip_columns,
        MUTUAL_COLUMNS + \
        ["clip_" + column for column in clip_columns if not column in MUTUAL_COLUMNS]
    ))
    clip_db.rename(columns=clip_columns, inplace=True)
    clip_groups = clip_db.groupby(MUTUAL_COLUMNS)

    rna22_groups = list(rna22_groups)

    processes = []
    for i in range(PROCESS_NUMBER):
        start = len(rna22_groups) // PROCESS_NUMBER * i
        end = len(rna22_groups) // PROCESS_NUMBER * (i + 1)

        if i == PROCESS_NUMBER - 1:
            end = len(rna22_groups)

        processes.append(Process(
            target=cross,
            args=(
                rna22_groups[start:end],
                clip_db,
                f"{PATH_TO_OUTPUT_DIR}{i}.tsv"
            )
        ))

    for i in range(PROCESS_NUMBER):
        processes[i].start()

    for i in range(PROCESS_NUMBER):
        processes[i].join()
        print(f"Process {i} finished")

    merge("/home/dude/huge/dude/rbp-miRNA/results/rna22_energy_seed_pvalue_clip.tsv",
        PATH_TO_OUTPUT_DIR)
