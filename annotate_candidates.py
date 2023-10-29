#!/usr/bin/env python

from gtfparse import read_gtf
import pandas as pd
import numpy as np

candidates = pd.read_csv("candidates-2020-12-18.from_spreadsheet.primary.tsv", sep="\t")

gtfFilename = f"/home/bullaugh/resources/mm10/gencode.vM22.annotation.gtf.gz"
gtf = read_gtf(gtfFilename);
exons = gtf[gtf["feature"] == "exon"]
genes = gtf[gtf["feature"] == "gene"]

# Identify exons with an endpoint shared with refStart
#sel = (sameChr["start"] == candidateLeft) | sameChr["end"] == candidateLeft
# Identify exons with an endpoint shared with refEnd
#sel = (sameChr["start"] == candidateRight) | sameChr["end"] == candidateRight

def zipStrTuples(a, b, sep=""):
    return([f"{a[i]}{sep}{b[i]}" for i in range(len(a))])

class Annotator:
    def __init__(self, candidate):
        self.chr = candidate["chr"]
        self.fields = {}
        self.left = candidate["refStart"]
        self.right = candidate["refEnd"]
        self.fieldNames = []
    def annotate(self):
        global exons
        global genes
        self.chrExons  = exons[exons["seqname"] == self.chr]
        self.chrGenes  = genes[genes["seqname"] == self.chr]
        # Return a list of genes that contain either the refStart or refEnd 
        sel = (((self.chrGenes["start"] <= self.left) & (self.left <= self.chrGenes["end"])) | \
            ((self.chrGenes["start"] <= self.right) & (self.right <= self.chrGenes["end"])))
        zipped = zipStrTuples(list(self.chrGenes["gene_name"][sel]), list(self.chrGenes["strand"][sel]))
        genesContainingCirc = ",".join(zipped)
        self.append(f"ContainingGenes", genesContainingCirc)
        # Find the closest exon defined by the exon start and the refStart
        self.nearestFeature("start", self.left, "exonLeftCircLeft")
        # Find the closest exon defined by the exon end and the refStart
        self.nearestFeature("end", self.left, "exonRightCircLeft")
        # Find the closest exon defined by the exon start and the refEnd
        self.nearestFeature("start", self.right, "exonLeftCircRight")
        # Find the closest exon defined by the exon end and the refEnd
        self.nearestFeature("end", self.right, "exonRightCircRight")
        # Identify the exon that contains refStart (if there is one) and for which
        # the exon start (or end) is closest to refStart
        sel = (self.chrExons["start"] <= self.left) & (self.left <= self.chrExons["end"])
        self.nearestFeature("start", self.left, "exonContainingCircLeftNearStart", sel)
        self.nearestFeature("end", self.left, "exonContainingCircLeftNearEnd", sel)
        # Same for containing refEnd
        sel = (self.chrExons["start"] <= self.right) & (self.right <= self.chrExons["end"])
        self.nearestFeature("start", self.right, "exonContainingCircRightNearStart", sel)
        self.nearestFeature("end", self.right, "exonContainingCircRightNearEnd", sel)
    def nearestFeature(self, gtfColumn, coord, label, sel=None):
        df = self.chrExons
        if not sel is None:
            if sum(sel) == 0:
                self.append(f"{label}GeneName", "NA")
                self.append(f"{label}ExonNumber", "NA")
                self.append(f"{label}LeftDist", "NA")
                return
            df = df[sel]
        dist = (df[gtfColumn] - coord).to_numpy()
        idx = np.absolute(dist).argmin()
        self.append(f"{label}GeneName", df.iloc[idx]["gene_name"])
        self.append(f"{label}ExonNumber", df.iloc[idx]["exon_number"])
        self.append(f"{label}LeftDist", dist[idx])
    def append(self, key, value):
        self.fieldNames.append(key)
        self.fields[key] = value
    def printHeader(self):
        line = ""
        sep = ""
        for key in self.fieldNames:
            line = f"{line}{sep}{key}"
            sep="\t"
        print(line)
    def print(self):
        line = ""
        sep = ""
        for key in self.fieldNames:
            line = f"{line}{sep}{self.fields[key]}"
            sep = "\t"
        print(line)

for i in range(len(candidates)):
    annotation = Annotator(candidates.iloc[i])
    annotation.annotate()
    if i == 0:
        annotation.printHeader()
    annotation.print()



