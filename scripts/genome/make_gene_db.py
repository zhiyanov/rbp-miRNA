import gffutils
import pandas as pd
import tqdm

db = gffutils.create_db(
    '/home/dude/huge/bulk/Ensembl_96/Homo_sapiens.GRCh38.96.gtf',
    dbfn='/home/dude/huge/dude/rbp-miRNA/genome/Homo_sapiens.GRCh38.96.db',
    verbose=True,
    merge_strategy='error',
    disable_infer_transcripts=True,
    disable_infer_genes=True,
)
