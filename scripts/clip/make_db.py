import pandas as pd
import tqdm
import os

PATH_TO_DATA = "/home/dude/huge/bulk/ENCONDE_eCLIP-seq/data"


samples = list(filter(
    lambda dir: os.path.isdir(f"{PATH_TO_DATA}/{dir}"),
    os.listdir(PATH_TO_DATA)))

result = pd.DataFrame(columns=[
    "sample", "rbp",
    "chr","strand",
    "start", "end",
    "log2(fold-enrichment)",
    "-log10(p-value)",
])

for sample in samples:
    tqdrator = tqdm.tqdm(os.listdir(f"{PATH_TO_DATA}/{sample}"))
    for rbp in tqdrator:
        peaks = pd.read_csv(
            f"{PATH_TO_DATA}/{sample}/{rbp}/narrow_peak.bed",
            header=None,
            sep="\t",
            names=[
                "chr", "start", "end",
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
        
        peaks["chr"] = peaks["chr"].str.lstrip("chr")
        peaks["sample"] = [sample] * len(peaks)
        peaks["rbp"] = [rbp] * len(peaks)
        
        result = pd.concat([result, peaks[
            result.columns.to_list()
        ]])

result.to_csv("../clip/clip_db.tsv", sep="\t", index=None)
