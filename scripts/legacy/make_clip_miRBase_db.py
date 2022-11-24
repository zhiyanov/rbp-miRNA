import gffutils
import pandas as pd
import tqdm
import json
import os

PATH_TO_DATA = "/home/dude/huge/bulk/ENCONDE_eCLIP-seq/data/HepG2"
# PATH_TO_DATA = "/home/dude/huge/bulk/ENCONDE_eCLIP-seq/data/K562"

mirbase_db = pd.read_csv("../miRBase/miRBase_22.1.tsv", sep="\t")
tqdrator = tqdm.tqdm(os.listdir(f"{PATH_TO_DATA}"))

result_table = []
for rbp in tqdrator:
    peaks = pd.read_csv(
        f"{PATH_TO_DATA}/{rbp}/narrow_peak.bed",
        header=None,
        sep="\t",
        names=[
            "chr", "start", "stop",
            "dataset_label", "1000", "strand",
            "log2(fold-enrichment)",
            "-log10(p-value)",
            "1-1", "2-1"
        ]
    )

    # stringent cutoff
    peaks = peaks.loc[
        (peaks["-log10(p-value)"] >= 5) & \
        (peaks["log2(fold-enrichment)"] >= 3)
    ]
    
    # peaks["chr"] = peaks["chr"].str.lstrip("chr")
    peaks["protein"] = [rbp] * len(peaks)
    
    peaks = peaks[[
        "protein",
        "chr", "strand", "start", "stop", 
        "log2(fold-enrichment)",
        "-log10(p-value)",
    ]]
    
    count = len(peaks)
    for i, row in peaks.iterrows():
        count -= 1
        
        chrom = row['chr']
        strand = row['strand']
        start, stop = row['start'], row['stop']

        condition = (mirbase_db["Chr"] == chrom) & (
            ((mirbase_db["Start"] >= start) & (mirbase_db["Start"] <= stop)) | \
            ((mirbase_db["End"] >= start) & (mirbase_db["End"] <= stop)) | \
            ((mirbase_db["Start"] <= start) & (mirbase_db["End"] >= start)) | \
            ((mirbase_db["Start"] <= stop) & (mirbase_db["End"] >= stop))
        ) & (mirbase_db["Strand"] == strand)
            
        for j, mirow in mirbase_db.loc[condition].iterrows():
            merged = list(row.values)
            merged.extend(list(mirow.values))
            result_table.append(merged)

        tqdrator.set_description(f"{len(result_table)}")

columns = [
    "protein",
    "bind_chr", "bind_strand", "bind_start", "bind_stop", 
    "log2(fold-enrichment)",
    "-log10(p-value)"
]
columns.extend(mirbase_db.columns.to_list())
result_table = pd.DataFrame(
    result_table,
    columns=columns
)

result_table.to_csv("../results/HepG2_clip_miRBase.csv", sep=",", index=None)
# result_table.to_csv("../results/K562_clip_miRBase.csv", sep=",", index=None)
