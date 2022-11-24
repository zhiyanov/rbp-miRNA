import gffutils
import pandas as pd
import tqdm
import json
import os

# PATH_TO_DATA = "/home/dude/huge/bulk/ENCONDE_eCLIP-seq/data/HepG2"
PATH_TO_DATA = "/home/dude/huge/bulk/ENCONDE_eCLIP-seq/data/K562"

FEATURE_TYPES = [
    "five_prime_utr",
    "three_prime_utr",
    "transcript",
    # "exon",
    # "start_codon",
    # "stop_codon",
    "CDS",
    # "Selenocysteine",
    # "gene",
]


gene_db = gffutils.FeatureDB('../genome/Homo_sapiens.GRCh38.96.db')
tqdrator = tqdm.tqdm(os.listdir(f"{PATH_TO_DATA}"))

result_dict = {}
for rbp in tqdrator:
    result_dict[rbp] = {}

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
    
    peaks["chr"] = peaks["chr"].str.lstrip("chr")
    peaks["protein"] = [rbp] * len(peaks)
    
    peaks = peaks[[
        "protein",
        "chr", "start", "stop", "strand",
        "log2(fold-enrichment)",
        "-log10(p-value)",
    ]]
    
    count = len(peaks)
    for i, row in peaks.iterrows():
        count -= 1

        region = f"{row['chr']}:{row['start']}-{row['stop']}"
        strand = row['strand']
        start, stop = row['start'], row['stop']
        
        for ftype in FEATURE_TYPES:
            result_dict[rbp][ftype] = []
            
            tqdrator.set_description(rbp + " " + str(count))

            for gene_structure in gene_db.region(region=region, featuretype=ftype, strand=strand):
                if (start < gene_structure.start) or (stop > gene_structure.stop):
                    continue

                result_dict[rbp][ftype].append({
                    "rbp_start": start,
                    "rbp_stop": stop,
                    "chr": row["chr"],
                    "strand": row["strand"],
                    "log2(fold-enrichment)": row["log2(fold-enrichment)"],
                    "-log10(p-value)": row["-log10(p-value)"],
                    # "exone_id": gene_structure["exone_id"],
                    "transcript_id": gene_structure["transcript_id"],
                    "gene_id": gene_structure["gene_id"],
                    "gene_name": gene_structure["gene_name"],
                    "start": gene_structure.start,
                    "stop": gene_structure.stop
                })

                print(ftype, result_dict[rbp][ftype][-1]["gene_name"])

    # with open('../clip/HepG2.json', 'w') as fp:
    with open('../clip/K562.json', 'w') as fp:
        json.dump(result_dict, fp)
    fp.close()
