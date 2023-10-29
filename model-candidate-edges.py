#!/usr/bin/env python

import argparse
import sys
import os
import pandas as pd
import tempfile
import subprocess
import stringdist
import random

mouseFasta = "/home/bullaugh/resources/mm10/GRCm38-M22/GRCm38.primary_assembly.genome.fa"

parser = argparse.ArgumentParser(description='look at the edges of candidates to look for matches')
parser.add_argument('-bed', type=str, help='Bed file of candidates', required=True)
parser.add_argument('-null', type=int, help='sample N samples from null distribution')
parser.add_argument('-k', type=int, help='kmer size for computing levenstein distance', default=10)
parser.add_argument('-rc', action='store_true', help='reverse complement one end')
parser.add_argument('-inside', action='store_true', help='consider both inside kmers instead')
parser.add_argument('-outside', action='store_true', help='consider both outside kmers instead')
parser.add_argument('-window', type=int, help='window within which to sample null start/stop', default=100)
parser.add_argument('-log', type=str, help='Write table of kmers to log file')
args = parser.parse_args()

def log(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def fail(msg):
    log(msg)
    sys.exit(1)

def getSequence(bed, fastaFn):
    # Write a temporary bed file for the sequence we want
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        try:
            content = "\n".join(bed["chrom"] + "\t" + bed["start"].astype(str) + "\t" + bed["stop"].astype(str))
            tmp.write((content + "\n").encode('utf8'))
            tmp.close()
            # Use the bedfile to get the sequence in tabular format
            cmd = ["bedtools", "getfasta", "-bed", tmp.name, "-tab", "-fi", fastaFn]
            with subprocess.Popen(cmd, stdout=subprocess.PIPE) as proc:
                return pd.read_csv(proc.stdout, sep="\t", names=("label", "seq"))["seq"]
        finally:
            os.remove(tmp.name)

def writeLogFile(df, args):
    with open(args.log, "w") as fp:
        print(args, file=fp)
        df.to_csv(fp, sep="\t")

def levenshteinDistance(seqPair, rc=False):
    dists = []
    for i in range(seqPair.shape[0]):
        left = seqPair["left"][i]
        right = seqPair["right"][i]
        if rc:
            right = revcomp(right)
        dist = stringdist.levenshtein(left, right)
        dists.append(dist)
    return(pd.Series(dists))

tr_table = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '|': '|'}

def revcomp(seq):
    complement = [tr_table[ch] for ch in seq]
    complement.reverse()
    return "".join(complement)

# Check arguments
if args.inside and args.outside:
    fail("cannot specify both -inside and -outside")

b = pd.read_csv(args.bed, sep="\t", names=("chrom", "start", "stop"))
k = args.k

if not args.null is None:
    N = b.shape[0]
    nullChrom = []
    nullStart = []
    nullStop = []
    for _ in range(args.null):
        row = random.randrange(N)
        start = b["start"][row]
        stop = b["stop"][row]
        startOffset = random.randrange(start - args.window, start + args.window)
        stopOffset = random.randrange(start - args.window, start + args.window)
        nullChrom.append(b["chrom"][row])
        nullStart.append(startOffset)
        nullStop.append(stopOffset)
    nullChrom = pd.Series(nullChrom).rename("chrom")
    nullStart = pd.Series(nullStart).rename("start")
    nullStop = pd.Series(nullStop).rename("stop")
    b = pd.concat([nullChrom, nullStart, nullStop], axis=1)

leftFlank = pd.concat([b["chrom"], b["start"]-k, b["start"].rename("stop")], axis=1)
leftInterior = pd.concat([b["chrom"], b["start"], (b["start"]+k).rename("stop")], axis=1)
rightInterior = pd.concat([b["chrom"], (b["stop"]-k).rename("start"), b["stop"]], axis=1)
rightFlank = pd.concat([b["chrom"], b["stop"].rename("start"), b["stop"]+k], axis=1)

distsList = []

needLeftFlank = not args.inside
needRightFlank = not args.inside
needLeftInterior = not args.outside
needRightInterior = not args.outside

if needLeftFlank:
    leftFlank = leftFlank.join(getSequence(leftFlank, mouseFasta))
if needRightInterior:
    rightInterior = rightInterior.join(getSequence(rightInterior, mouseFasta))

if needLeftInterior:
    leftInterior = leftInterior.join(getSequence(leftInterior, mouseFasta))
if needRightFlank:
    rightFlank = rightFlank.join(getSequence(rightFlank, mouseFasta))

if args.inside:
    seqPair = pd.concat([leftInterior["seq"].rename("left"), rightInterior["seq"].rename("right")], axis=1)
    dists = levenshteinDistance(seqPair, args.rc)
    distsList.append(dists)
    if args.log:
        writeLogFile(pd.concat([seqPair, dists], axis=1), args)
elif args.outside:
    seqPair = pd.concat([leftFlank["seq"].rename("left"), rightFlank["seq"].rename("right")], axis=1)
    distsList.append(levenshteinDistance(seqPair, args.rc))
    if args.log:
        writeLogFile(pd.concat([seqPair, dists], axis=1), args)
else:
    seqPair1 = pd.concat([leftFlank["seq"].rename("left"), rightInterior["seq"].rename("right")], axis=1)
    dists1 = levenshteinDistance(seqPair1, args.rc)
    distsList.append(dists1)
    seqPair2 = pd.concat([leftInterior["seq"].rename("left"), rightFlank["seq"].rename("right")], axis=1)
    dists2 = levenshteinDistance(seqPair2, args.rc)
    distsList.append(dists2)
    if args.log:
        writeLogFile(pd.concat([seqPair1, dists1, seqPair2, dists2], axis=1), args)

mins = pd.concat(distsList, axis=1).min(axis=1)
counts = [0] * (k+1)
dist = pd.Series([i for i in range(k+1)])
for _, value in mins.items():
    counts[value] += 1

counts = pd.Series(counts)
print("\n".join(dist.astype(str) + "\t" + counts.astype(str)))



