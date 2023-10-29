#!/usr/bin/env python

# This script takes a list of chimeric junction files on STDIN

import sys
import re
import uuid
import subprocess
import pickle
import argparse
import os

parser = argparse.ArgumentParser(description='List the junctions that are compatible with cirularization.')
parser.add_argument('-run', type=str, help='which special run')
args = parser.parse_args()

def log(msg):
    print(msg, file=sys.stderr, flush=True)

def splitCigar(cigar):
    pieces = re.findall("-?[0-9]+[pA-Z]", cigar)
    if "".join(pieces) != cigar:
        raise RuntimeError(f"Malformed {cigar}")
    return [(p[-1], int(p[0:-1])) for p in pieces]

candidates = 0
segsOutOfOrder = 0
missingGap = 0
entries = 0
junctionsBetweenReads = 0
seg1StartTooLow = 0
backwards = 0
comments = 0
otherChr = 0
opposingStrands = 0
notMainChr = 0
tooFar = 0
pFirstCount = 0
pSecondCount = 0
docs = {}
mitochondria = 0
hotspot = 0
readToDocIds = {}
backSplices = {}

validRealignmentPaths = ("mouse-dendrites",)

inputFiles = 0
for line in sys.stdin:
    line = line.rstrip()
    log(f"processing: {line}")
    pathParts = line.split("/")
    barcode = None
    if args.run in ("mouse-dendrites",):
        if len(pathParts) != 4:
            raise RuntimeError(f"malformed path: {line}")
    else:
        raise RuntimeError(f"unsupported run {args.run}")
    if pathParts[0] != "realignment" or (not (pathParts[1] in validRealignmentPaths)):
        raise RuntimeError(f"incorrect path prefix: {pathParts[0]}: {line}")
    sample = pathParts[2]
    chimericFn = line
    inputFiles += 1
    readNames = []
    with open(chimericFn, "r") as fp:
        lineNum = 0
        for line in fp:
            line = line.rstrip()
            lineNum += 1
            if lineNum == 1 and line[0:10] == "chr_donorA":
                # header line
                continue
            if line[0] == "#":
                comments += 1
                continue
            fields = line.split("\t")
            # the last column, readgrp, might not be present
            if len(fields) != 21 and len(fields) != 20:
                raise RuntimeError(f"malformed line ({len(fields)}) {chimericFn}:{lineNum}: {line}")
            entries += 1
            (chrom, dIntronStart, dStrand, aChr, aIntronStart, aStrand, junctionType, rep1, rep2, read,
                seg1Start, cigar1, seg2Start, cigar2, numChimAln, maxPossAlnScore,
                nonChimAlnScore, thisChimAlnScore, bestallChimAlnScore, PEmergedBool) = fields[0:20]
            dIntronStart = int(dIntronStart)
            aIntronStart = int(aIntronStart)
            seg1Start = int(seg1Start)
            seg2Start = int(seg2Start)
            junctionType = int(junctionType)
            # Ignore lines we don't expect to be back-spliced junctions.
            if dStrand != aStrand:
                opposingStrands += 1
                continue
            if chrom != aChr:
                otherChr += 1
                continue
            if re.search("_", chrom) or re.search("ERCC", chrom):
                notMainChr += 1
                continue
            if chrom == "chrM":
                mitochondria += 1
                continue
            if junctionType < 0:
                junctionsBetweenReads += 1
                continue
            if abs(dIntronStart - aIntronStart) > 5000:
                tooFar += 1
                continue
            if dStrand == "+":
                if dIntronStart < aIntronStart:
                    backwards += 1
                    continue
                if seg1Start < aIntronStart:
                    seg1StartTooLow += 1
                    continue
                if seg1Start < seg2Start:
                    segsOutOfOrder += 1
                    continue
            else:
                if dIntronStart > aIntronStart:
                    backwards += 1
                    continue
                if seg1Start > aIntronStart:
                    seg1StartTooLow += 1
                    continue
                if seg1Start > seg2Start:
                    segsOutOfOrder += 1
                    continue
            cigar1Split = splitCigar(cigar1)
            cigar2Split = splitCigar(cigar2)
            pFirst = None
            if any([t[0] == "p" for t in cigar1Split]):
                # first segment contains parts of both reads.
                pFirst = True
                pFirstCount += 1
            elif any([t[0] == "p" for t in cigar2Split]):
                # second segment contains parts of both reads.
                pFirst = False
                pSecondCount += 1
            else:
                missingGap += 1
                continue
            # The two reads are split across three chunks. The chunks are ordered
            # according to the order across the chimeric molecule, not the order in
            # the reference genome.
            # TODO: Do I need to do something different for the minus strand when assembling chunks?
            chunks = ({}, {}, {})
            if pFirst:
                chunks[0]['start'] = seg1Start
                chunks[1]['start'] = seg1Start
                for (op, n) in cigar1Split:
                    if op == "p": break
                    chunks[1]['start'] += n
                chunks[2]['start'] = seg2Start
            else:
                chunks[0]['start'] = seg1Start
                chunks[1]['start'] = seg2Start
                chunks[2]['start'] = seg2Start
                for (op, n) in cigar1Split:
                    if op == "p": break
                    chunks[2]['start'] += n
            bedKey = None
            if dStrand == "+":
                output = f"{chrom}\t{aIntronStart}\t{dIntronStart-1}"
            else:
                output = f"{chrom}\t{dIntronStart}\t{aIntronStart-1}"
            output = f"{output}\t{sample}\t{pFirst}\t{dStrand}\t{junctionType}\t{read}"
            output = f"{output}\t{seg1Start}\t{cigar1}\t{seg2Start}\t{cigar2}\t{nonChimAlnScore}\t{thisChimAlnScore}\t{bestallChimAlnScore}\t{PEmergedBool}"
            print(output)

# END    
