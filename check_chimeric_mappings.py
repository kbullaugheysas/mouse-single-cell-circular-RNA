#!/usr/bin/env python

# Filter chimeric mappings stdin to include only those compatible with circularization

import argparse
import sys
import os
import re
import uuid
from subprocess import Popen, PIPE
import math

parser = argparse.ArgumentParser(description='Filter chimeric mappings for those compatible with circularization')
parser.add_argument('-limit', type=int, help='limit number of lines of input')
parser.add_argument('-reads', type=str, help='file name containing tab-separated reads', required=True)
parser.add_argument('-refs', type=str, help='file name containing tab-separated reference sequences')
params = parser.parse_args()

totalMatches = 0
totalMismatches = 0
trTable = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

def log(msg):
    print(msg, file=sys.stderr)

def revcomp(seq):
    complement = [trTable[ch] for ch in seq]
    complement.reverse()
    return "".join(complement)

def splitCigar(cigarStr):
    cigarPieces = re.sub(r"([A-Zp])", ":\\1 ", cigarStr).rstrip().split(" ")
    cigarTuples = [x.split(":") for x in cigarPieces]
    return [(int(n), op) for (n,op) in cigarTuples]

def countMatches(x, y):
    if len(x) != len(y):
        log(x)
        log(y)
        raise RuntimeError("can only count matches in equal-length sequences")
    return len([i for i in range(len(x)) if x[i] == y[i]])

def preliminaryCheck(r):
    (dChr, dPos, dSt, aChr, aPos, aSt, junc) = (
            r['dChr'], r['dPos'], r['dSt'], r['aChr'], r['aPos'], r['aSt'], r['junc'])
    if aChr != dChr:
        return "mismatch-chromosome"
    if dSt != aSt:
        return "mismatch-strand"
    if junc == -1:
        return "encompassing-junction"
    if not re.match("^chr[0-9XY]*$", dChr):
        return "non-primary-chromosome"
    return "okay"

def processRecord(r):
    global totalMatches, totalMismatches
    log(f"processing record:")
    log(r['line'])
    (dChr, dPos, dSt, aChr, aPos, aSt, junc, seg1Pos, cig1, seg2Pos, cig2) = (
            r['dChr'], r['dPos'], r['dSt'], r['aChr'], r['aPos'], r['aSt'], r['junc'], r['seg1Pos'],
            r['cig1'], r['seg2Pos'], r['cig2'])
    ref = r['ref']
    log(f"ref: {ref}")
    mate1 = r['mate1']
    mate2 = r['mate2']
    if r['dSt'] == "-":
        # Mate 2 has already been reverse complemented. So the two mates are the same strand, but to make
        # coordinate math easier, we work with reference stand bases by reverse complementing both mates
        # when on the minus strand.
        mate1 = revcomp(mate1)
        mate2 = revcomp(mate2)
    lowBound = None
    highBound = None
    # Find bounds of the aligned bases
    if dPos < aPos:
        lowBound = dPos
        highBound = aPos
    else:
        lowBound = aPos
        highBound = dPos
    # Check that segments don't go out of bounds
    if seg1Pos < lowBound:
        return "seg1Pos-low"
    if highBound < seg1Pos:
        return "seg1Pos-high"
    if seg2Pos < lowBound:
        return "seg2Pos-low"
    if highBound < seg2Pos:
        return "seg2Pos-high"
    r['pFirst'] = "p" in cig1
    cig1 = splitCigar(cig1)
    cig2 = splitCigar(cig2)
    c = 0
    currentRead = None
    nextRead = None
    if r['pFirst']:
        log("starting with mate 2")
        currentRead = mate2
        nextRead = mate1
    else:
        log("starting with mate 2")
        currentRead = mate1
        nextRead = mate2
    refStart = r['refStart']
    log(f"refStart: {refStart}")
    for (refCoord, cig) in ((seg1Pos, cig1), (seg2Pos, cig2)):
        c += 1
        offset = 0
        log(f"cig{c}={cig}, refCoord={refCoord}")
        log(f"current read: {currentRead}")
        log(f"next read: {nextRead}")
        for i, (n,op) in enumerate(cig):
            if op == "S":
                offset += n
            elif op == "M":
                if len(currentRead) < offset+n:
                    raise RuntimeError(f"current read too short")
                if len(ref) < refCoord - refStart + n:
                    raise RuntimeError(f"ref sequence is too short")
                segMatch = currentRead[offset:(offset+n)]
                refMatch = ref[(refCoord-refStart):(refCoord-refStart+n)]
                log(f"segMatch: {segMatch}")
                log(f"refMatch: {refMatch}")
                matches = countMatches(refMatch, segMatch)
                totalMatches += matches
                mismatches = n - matches
                totalMismatches += mismatches
                if mismatches > 5:
                    log(f"warning: got {mismatches} mismatches, which is too many")
                    return "too-many-mismatches"
                log(f"matched: {matches} out of {n}")
                offset += n
                refCoord += n
            elif op == "I":
                offset += n
            elif op == "D":
                refCoord += n
            elif op == "N":
                refCoord += n
            elif op == "p":
                if offset != len(currentRead):
                    #log(f"warning: only consumed {offset} out of {len(currentRead)} bases in mate at insert")
                    return "incomplete"
                refCoord += n
                # Swap mates when we get to the paired-end fragment insertion
                (currentRead, nextRead) = (nextRead, currentRead)
                offset = 0
            else:
                raise RuntimeError(f"invalid opcode {op} at line {lineNum}")
            if refCoord < lowBound:
                return f"cig{c}-low"
            if highBound < refCoord:
                return "cig{c}-high"
        # Swap mates again. This is only relevant after the first segment and if we're on the minus strand, when there's still another segment.
        if c == 1 and r['dSt'] == "-":
            if offset != len(currentRead):
                #log(f"warning: only consumed {offset} out of {len(currentRead)} bases in mate at end of read")
                return "incomplete"
            (currentRead, nextRead) = (nextRead, currentRead)
            log(f"switching to read {currentRead}")
    if r['dSt'] == "-":
        if refCoord != refStart + r['length']:
            log(f"warning: segment doesn't end at end of ref sequence")
    return "okay"

lineNum = 0
skipped = 0
okay = 0
records = []
for line in sys.stdin:
    if (not params.limit is None) and lineNum >= params.limit:
        log(f"termiating early after {lineNum} lines due to -limit")
        break
    lineNum += 1
    line = line.rstrip()
    fields = line.split('\t')
    if len(fields) < 14:
        raise RuntimeError(f"expecting at least 22 fields at line {lineNum} got {len(fields)}")
    (dChr, dPos, dSt, aChr, aPos, aSt, junc, rep1, rep2, read, seg1Pos, cig1, seg2Pos, cig2) = fields[0:14]
    junc = int(junc)
    record = {}
    record['dChr'] = dChr
    record['dSt'] = dSt
    record['aChr'] = aChr
    record['aSt'] = aSt
    record['junc'] = junc
    record['read'] = read
    record['cig1'] = cig1
    record['cig2'] = cig2
    # Convert coordinates to zero-indexed
    record['dPos'] = int(dPos)-1
    record['aPos'] = int(aPos)-1
    record['seg1Pos'] = int(seg1Pos)-1
    record['seg2Pos'] = int(seg2Pos)-1
    record['line'] = line
    result = preliminaryCheck(record)
    if result == "okay":
        records.append(record)
    else:
        log(f"skipping due to {result}:")
        log(line)
        skipped += 1
        continue
    record['id'] = str(uuid.uuid4())
    left = record['dPos']+1
    right = record['aPos']
    if record['aPos'] < record['dPos']:
        left = record['aPos']+1
        right = record['dPos']
    bedEntryName = f"{record['dChr']}:{left}-{right}"
    record['bedName'] = bedEntryName
    record['refStart'] = left
    record['length'] = right - left

log(f"initially accepting {len(records)} entries")

def readBedtoolsOutput(fp, saveFp=None):
    refSeqs = {}
    lineNum = 0
    for line in fp:
        lineNum += 1
        if not type(line) is str:
            line = line.decode("utf-8")
        line = line.rstrip()
        fields = line.split('\t')
        if len(fields) != 2:
            raise RuntimeError(f"expecting 2 bed fields at line {lineNum} got {len(fields)}")
        (bedName, seq) = fields
        refSeqs[bedName] = seq.upper()
        if not saveFp is None:
            print(line, file=saveFp)
    return refSeqs

# If we're given a -refs argument but the file doesn't exist, we use bedtools to generate it. If this
# argument is given and the file exits, the contents of the file are used. If the argument is not given
# then we read directly from bedtools and don't save the result.
refsFp = None
refsMode = "tmp"
refSeqs = None
if not params.refs is None:
    if os.path.isfile(params.refs):
        refsFp = open(params.refs, "r")
        refsMode = "read"
        log(f"reading bedtools output from {params.refs}")
    else:
        refsMode = "write"
        log(f"writing bedtools output to {params.refs}")
        refsFp = open(params.refs, "w")

if refsMode == "read":
    refSeqs = readBedtoolsOutput(refsFp)
else:
    # Write a temporary bed file
    bedFn = f"/home/bullaugh/tmp/{uuid.uuid4()}.bed"
    log(f"bed: {bedFn}")
    with open(bedFn, "w") as bedFp:
        for r in records:
            left = r['refStart']
            right = left + r['length']
            entry = f"{r['dChr']}\t{left}\t{right}"
            print(entry, file=bedFp)
    # Fetch the reads
    genomeFa = "/home/bullaugh/resources/mm10/chromosomes/mm10.fa"
    with Popen(["bedtools", "getfasta", "-fi", genomeFa, "-bed", bedFn, "-tab"], stdout=PIPE) as proc:
        if refsMode == "write":
            refSeqs = readBedtoolsOutput(proc.stdout, refsFp)
        elif refsMode == "tmp":
            refSeqs = readBedtoolsOutput(proc.stdout)
        else:
            raise RuntimeError(f"Invalid refsMode {refsMode}")
if not refsFp is None:
    refsFp.close()

# Read in the reads data. This was generated with fqseq -names, but the third column has
# already been reverse complemented.
lineNum = 0
reads = {}
with open(params.reads, "r") as fp:
    for line in fp:
        lineNum += 1
        line = line.rstrip()
        fields = line.split("\t")
        if len(fields) != 3:
            raise RuntimeError(f"expecting 3 reads fields at line {lineNum} got {len(fields)}")
        (read, left, right) = fields
        reads[read] = (left, right)

log(f"read {len(reads)} reads")

for record in records:
    record['ref'] = refSeqs[record['bedName']]
    (mate1, mate2) = reads[record['read']]
    record['mate1'] = mate1
    record['mate2'] = mate2
    result = processRecord(record)
    if result != "okay":
        log(f"skipping due to {result}:")
        log(line)
        skipped += 1
        continue
    log("okay")
    refLen = record['length']
    if refLen != len(record['ref']):
        raise RuntimeError("reflen mismatch")
    halfLen = 100
    if refLen < halfLen*2:
        halfLen = math.floor(refLen/2)
    fragment = record['ref'][(refLen-halfLen):refLen] + record['ref'][0:halfLen]
    strand = record['dSt']
    if strand == "-":
        fragment = revcomp(fragment)
    print(f"{record['dChr']}\t{record['refStart']}\t{record['refStart']+refLen}\t{strand}\t{fragment}")
    okay += 1

log(f"processed: {lineNum}")
log(f"skipped: {skipped}")
log(f"okay: {okay}")
log(f"mismatches: {totalMismatches}")
log(f"matches: {totalMatches}")
if totalMatches+totalMismatches > 0:
    log(f"mismatch rate: {totalMismatches/(totalMismatches+totalMatches)}")

# END
