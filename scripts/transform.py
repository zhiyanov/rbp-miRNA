import pandas as pd
import json
import gffutils

# gtf_db = gffutils.FeatureDB('../genome/Homo_sapiens.GRCh38.96.db')

def trans2gene(chrom, strand, coord, gtf_db, transcript_id):
    transcript = gtf_db[transcript_id]
    
    if not (transcript.strand == strand):
        return -1
    if not (transcript.seqid == chrom):
        return -2
    
    position = 0
    for exon in gtf_db.children(
            transcript,
            featuretype="exon",
            order_by="start"):
        exon_length = exon.stop - exon.start + 1
        if (coord > position) and (coord <= position + exon_length):
           return exon.start + (coord - position) - 1
       
        position += exon.stop - exon.start + 1

    return -3

def gene2trans(chrom, strand, coord, gtf_db, transcript_id=None):
    if transcript_id == None:
        region = f"{chrom}:{coord}-{coord}"
        transcripts = gtf_db.region(
                region=region,
                strand=strand,
                featuretype="transcript")
        if not len(transcripts):
            return -1
        transcript = transcripts[0]
    else:
        transcript = gtf_db[transcript_id]
    
    if not (transcript.strand == strand):
        return -2
    if not (transcript.seqid == chrom):
        return -3
    if coord < transcript.start:
        return -4
    if coord > transcript.stop:
        return -5

    position = 0
    for exon in gtf_db.children(
            transcript,
            featuretype="exon",
            order_by="start"):
        exon_length = exon.stop - exon.start + 1
        if (coord >= exon.start) and (coord <= exon.stop):
           return position + coord - exon.start + 1
       
        position += exon_length

    return -6
