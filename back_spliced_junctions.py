#!/usr/bin/env python

# This script takes a list of chimeric junction files on STDIN

import sys
import re
import uuid
import subprocess
import pickle
import argparse
import os

genomeConfig = {}
genomeConfig['rat'] = {}
genomeConfig['rat']['genome'] = "/home/bullaugh/resources/rn6.fa"
genomeConfig['rat']['hotspots'] = []

genomeConfig['human'] = {}
genomeConfig['human']['genome'] = "/home/bullaugh/resources/hg38.fa"
genomeConfig['human']['hotspots'] = []

genomeConfig['mouse'] = {}
genomeConfig['mouse']['genome'] = "/home/bullaugh/resources/mm10/chromosomes/mm10.fa"
genomeConfig['mouse']['hotspots'] = [
        {'chr': 'chr2', 'left': 98660000, 'right': 98668000},
        {'chr': 'chr14', 'left': 19415000, 'right': 19419000}]

parser = argparse.ArgumentParser(description='Assemble evidence of cirular RNAs.')
parser.add_argument('-pickle', type=str, help='file to pickle results to', required=True)
parser.add_argument('-hotspots-only', action='store_true', dest='hotspots', help='only include notspots instead of excluding them')
parser.add_argument('-no-filter-hotspots', action='store_true', dest='noFilterHotspots', help="don't filter out hostspot candidates")
parser.add_argument('-species', type=str, help='species can be rat|human|mouse', required=True)
parser.add_argument('-run', type=str, help='which special run')
args = parser.parse_args()
print(args)

def log(msg):
    print(msg, file=sys.stderr, flush=True)

class DiscreteDistribution:
    def __init__(self, min=None, max=None):
        self.distr = {}
        self.min = min
        self.max = max
    def add(self, value):
        if len(self.distr) == 0:
            self.min = value
            self.max = value
        if not value in self.distr:
            self.distr[value] = 0
        self.distr[value] += 1
        if value < self.min:
            self.min = value
        if self.max < value:
            self.max = value
    def print(self, min=None, max=None):
        if min is None:
            min = self.min
        if max is None:
            max = self.max
        for i in range(min, max+1):
            count = self.distr.get(i, 0)
            print(f"{i}: {count}")

def splitCigar(cigar):
    pieces = re.findall("-?[0-9]+[pA-Z]", cigar)
    if "".join(pieces) != cigar:
        raise RuntimeError(f"Malformed {cigar}")
    return [(p[-1], int(p[0:-1])) for p in pieces]
def joinCigar(cigar):
    return "".join([b + str(n) for b,n in cigar])

def inHotspot(d, a, low, high):
    return (low <= a and a <= high) or (low <= d and d <= high)

if not args.species in genomeConfig:
    log(f"unsupported species: '{args.species}'")
    sys.exit(1)

genome = genomeConfig[args.species]['genome']
if not os.path.exists(genome):
    log(f"can't find genome file '{genome}'")
    sys.exit(1)

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

validRealignmentPaths = ("mm10-filtered-fastq", "mm10-filtered-fastq-STAR-2.7.2b", "TiPeR_and_SCAP", "mouse-dendrites", "mouse-astrocytes-and-neurons")

inputFiles = 0
tmpFilename = f"tmp-{uuid.uuid4()}.bed"
with open(tmpFilename, "w") as bed:
    log(f"writing to tmp file {tmpFilename}")
    for line in sys.stdin:
        line = line.rstrip()
        log(f"processing: {line}")
        pathParts = line.split("/")
        barcode = None
        if args.run in ("mouse-dendrites", "mouse-astrocytes-and-neurons"):
            if len(pathParts) != 4:
                raise RuntimeError(f"malformed path: {line}")
        elif args.run == "unfiltered":
            if len(pathParts) != 4:
                raise RuntimeError(f"malformed path: {line}")
        else:
            if len(pathParts) != 5:
                raise RuntimeError(f"malformed path: {line}")
            barcode = pathParts[3]
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
                isHotspot = False
                for hs in genomeConfig[args.species]['hotspots']:
                    if chrom == hs['chr'] and (not args.noFilterHotspots):
                        isHotspot = isHotspot or inHotspot(dIntronStart, aIntronStart, hs['left'], hs['right'])
                if args.hotspots:
                    if isHotspot:
                        hotspot += 1
                    else:
                        # We're focusing exclusively on hotspots, and this isn't one, so skip it.
                        continue
                else:
                    if isHotspot:
                        # We're excluding hotspots and this is a hotspot.
                        hotspot += 1
                        continue
                cigar1 = splitCigar(cigar1)
                cigar2 = splitCigar(cigar2)
                pFirst = None
                if any([t[0] == "p" for t in cigar1]):
                    # first segment contains parts of both reads.
                    pFirst = True
                    pFirstCount += 1
                elif any([t[0] == "p" for t in cigar2]):
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
                    for (op, n) in cigar1:
                        if op == "p": break
                        chunks[1]['start'] += n
                    chunks[2]['start'] = seg2Start
                else:
                    chunks[0]['start'] = seg1Start
                    chunks[1]['start'] = seg2Start
                    chunks[2]['start'] = seg2Start
                    for (op, n) in cigar1:
                        if op == "p": break
                        chunks[2]['start'] += n
                docId = str(uuid.uuid4())
                junctionKey = None
                bedKey = None
                if dStrand == "+":
                    junctionKey = (chrom, aIntronStart, dIntronStart)
                    bedKey = f"{chrom}\t{aIntronStart}\t{dIntronStart-1}"
                else:
                    junctionKey = (chrom, dIntronStart, aIntronStart)
                    bedKey = f"{chrom}\t{dIntronStart}\t{aIntronStart-1}"
                splice = None
                if junctionKey in backSplices:
                    splice = backSplices[junctionKey]
                else:
                    seqId = str(uuid.uuid4())
                    print(f"{bedKey}\t{seqId}", file=bed)
                    splice = {
                        'chr': chrom,
                        'dIntronStart':  dIntronStart,
                        'aIntronStart':  aIntronStart,
                        'docIds': [],
                        'seqId': seqId,
                    }
                    backSplices[junctionKey] = splice
                splice['docIds'].append(docId)
                doc = {
                    'id': docId,
                    'sample': sample,
                    'barcode': barcode,
                    'chunks': chunks,
                    'line': line,
                    'read': read,
                    'pFirst': pFirst,
                    'strand':  dStrand,
                    'junctionType':  junctionType,
                    'rep1':  rep1,
                    'rep2':  rep2,
                    'read':  read,
                    'seg1Start':  seg1Start,
                    'cigar1':  cigar1,
                    'seg2Start':  seg2Start,
                    'cigar2':  cigar2,
                    'nonChimAlnScore': nonChimAlnScore,
                    'thisChimAlnScore': thisChimAlnScore,
                    'bestallChimAlnScore': bestallChimAlnScore,
                    'PEmergedBool': PEmergedBool,
                }
                docs[docId] = doc
                if not read in readToDocIds:
                    readToDocIds[read] = []
                readToDocIds[read].append(docId)
                readNames.append(read)
                candidates += 1

        # Make a list of all the read names we want.
        readNamesData = "\n".join(readNames)
        # For mouse-astrocytes-and-neurons, we pre-filtered the fastq files for the chimeric reads.
        # Use fqfilter to extract just the reads we need
        if args.run == "mouse-astrocytes-and-neurons":
            subdir = f"{args.run}/chimeric-reads"
            cmd = ["zcat", f"{subdir}/{sample}.tsv.gz"]
            procInput = None
        elif args.run == "mouse-dendrites":
            subdir = f"{args.run}/filtered-fastq"
            read1Fn = f"{subdir}/{sample}_1.fq.gz"
            read2Fn = f"{subdir}/{sample}_2.fq.gz"
            cmd = ["fqfilter", "-short-name", "-reads", "stdin", "-tab", read1Fn, read2Fn]
            procInput = subprocess.PIPE
        elif args.run == "unfiltered":
            read1Fn = f"unfiltered-fastq/{sample}.unaligned_1.fq.gz"
            read2Fn = f"unfiltered-fastq/{sample}.unaligned_2.fq.gz"
            cmd = ["fqfilter", "-short-name", "-reads", "stdin", "-tab", read1Fn, read2Fn]
            procInput = subprocess.PIPE
        else:
            subdir = f"filtered-fastq/{sample}/{barcode}"
            read1Fn = f"{subdir}/filtered_1.fq"
            read2Fn = f"{subdir}/filtered_2.fq"
            cmd = ["fqfilter", "-short-name", "-reads", "stdin", "-tab", read1Fn, read2Fn]
            procInput = subprocess.PIPE

        with subprocess.Popen(cmd, stdin=procInput, stdout=subprocess.PIPE) as proc:
            if not procInput is None:
                proc.stdin.write(readNamesData.encode('utf-8'))
                proc.stdin.close()
            readCount = 0
            extraReads = 0
            for line in proc.stdout:
                line = line.decode('utf-8').rstrip()
                fields = line.split("\t")
                if len(fields) != 3:
                    raise RuntimeError(f"malformed fqfilter output: {line}")
                (name, read1, read2) = fields
                if name in readToDocIds:
                    for docId in readToDocIds[name]:
                        doc = docs[docId]
                        doc['read1'] = read1
                        doc['read2'] = read2
                    readCount += 1
                else:
                    extraReads += 1
            log(f"found {extraReads} extra reads in fqfilter output")
            log(f"found {readCount} reads from {sample}/{barcode}")


log(f"found {inputFiles} input files containing {comments} comments and {entries} non-comment entries")
log(f"missing gap between segments: {missingGap}")
log(f"segs out of order {segsOutOfOrder} times")
log(f"seg1 split {pFirstCount} times and seg2 split {pSecondCount} times")
log(f"excluding {otherChr} mappings on different chromosomes")
log(f"excluding {opposingStrands} mappings on opposing strands")
log(f"excluding {notMainChr} mappings from non-main chromosomes")
log(f"excluding {junctionsBetweenReads} mappings where the junction was not hit by a read")
log(f"excluding {seg1StartTooLow} mappings because seg1 start too early")
log(f"excluding {backwards} mappings because donar and acceptor are backwards")
log(f"excluding {tooFar} mappings because donor and acceptor are too far apart")
log(f"excluding {mitochondria} chrM mappings")
log(f"excluding {hotspot} mappings from hotspots")
log(f"identified {len(backSplices)} distinct backsplice junctions")
log(f"found {candidates} candidates")

hasSeq = len([d for d in docs.values() if "read1" in d])
if hasSeq == len(docs):
    log("all docs have read sequences")
else:
    missing = len(docs) - hasSeq
    log(f"Warning: {missing} docs are missing read sequences")

log("extracting fasta sequence")
cmd = ["bedtools", "getfasta", "-fi", genome,
        "-bed", tmpFilename, "-tab", "-name"]
dna = {}
with subprocess.Popen(cmd, stdout=subprocess.PIPE) as proc:
    for line in proc.stdout:
        line = line.decode('utf-8').rstrip()
        fields = line.split("\t")
        if len(fields) != 2:
            raise RuntimeError("Malformed bedtools getfasta tabular output")
        (seqId, chunk) = fields
        if seqId in dna:
            raise RuntimeError(f"Already saw seqId {seqId}")
        dna[seqId] = chunk

log("making a backup pickle")
with open(args.pickle, "wb") as fp:
    pickle.dump({'docs': docs, 'backSplices': backSplices}, fp, pickle.HIGHEST_PROTOCOL)

# Put the sequence into the backSplice dicts
for splice in backSplices.values():
    seqId = splice['seqId']
    if seqId in dna:
        chunk = dna[seqId]
        splice['dna'] = chunk
    else:
        log(f"Warning: failed to find dna for {seqId}")

log("making a pickle")
with open(args.pickle, "wb") as fp:
    pickle.dump({'docs': docs, 'backSplices': backSplices}, fp, pickle.HIGHEST_PROTOCOL)

backSpliceCountDistr = DiscreteDistribution()
for key,splice in backSplices.items():
    backSpliceCountDistr.add(len(splice['docIds']))

log("distribution of back-splice hits:")
backSpliceCountDistr.print(max=20)

# END    

