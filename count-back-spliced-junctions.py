#!/usr/bin/env python

import sys
import re
import uuid
import subprocess
import pickle
import argparse
import zlib
import os
import math
import random
import cProfile

random.seed(4)

def log(msg):
    print(msg, file=sys.stderr)
    sys.stderr.flush()

def compressionRatio(s):
    return len(s) / len(zlib.compress(s.encode('utf-8')))

class SoftClippingNotAllowed(Exception):
    def __init__( self, msg):
        Exception.__init__(self, msg)
class ShouldBeAtEndOfSequence(Exception):
    def __init__( self, msg):
        Exception.__init__(self, msg)

dnaPairings = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
iupac = {
    ("A","G"): "R",
	("G","A"): "R",
	("C","T"): "Y",
	("T","C"): "Y",
	("G","C"): "S",
	("C","G"): "S",
	("A","T"): "W",
	("T","A"): "W",
	("G","T"): "K",
	("T","G"): "K",
	("A","C"): "M",
	("C","A"): "M",
}

def revcomp(dna):
    return "".join(list(map(lambda x: dnaPairings[x], dna))[::-1])

class Column:
    def __init__(self, n):
        self.length = n
        self.col = [" "] * n
        self.warning = " "
        self.offset = 0
        self.read = None
        self.conflicts = 0
        self.matches = None
        self.misses = None
        self.encompassing = False
        self.overlapping = 0
        self.aligned = 0
        self.sample = None
        self.barcode = None
        self.otherStats = {}
        self.starts = []
        self.ends = []
        self.junctionType = None
    def addOtherStats(self, doc):
        self.otherStats['nonChimAlnScore'] = doc['nonChimAlnScore']
        self.otherStats['thisChimAlnScore'] = doc['thisChimAlnScore']
        self.otherStats['bestallChimAlnScore'] = doc['bestallChimAlnScore']
        self.otherStats['PEmergedBool'] = doc['PEmergedBool']
        #self.otherStats['readgrp'] = doc['readgrp']
        self.rep1 = doc['rep1']
        self.rep2 = doc['rep2']
        self.junctionType = doc['junctionType']
    def warn(self, message):
        self.log(message)
        self.warning = "W"
    def placeLargeGap(self, n):
        # This is a large gap. Could be a spliced intron. Mark the beginning and end with '*'
        if n <= 0:
            raise RuntimeError(f"N cigar term has count {n}")
        # Place the initial star
        self.place("*")
        if n >= 2:
            # Skip the intervening slots
            self.advance(n-2)
        if n > 1:
            # place the final star
            self.place("*")
    def place(self, thing, segment=None):
        if not segment is None:
            # If we haven't recorded a start or end for this segment yet, note them now.
            if len(self.starts) == segment:
                self.starts.append(self.offset)
                self.ends.append(self.offset)
            # Update the end point for the end.
            self.ends[segment] = self.offset
        existing = self.col[self.offset]
        if existing == " ":
            # Nothing here, just place new thing
            self.col[self.offset] = thing
            if thing in dnaPairings and thing != "N":
                self.aligned += 1
        else:
            self.overlapping += 1
            if existing != thing:
                if existing in dnaPairings and thing in dnaPairings and (existing,thing) in iupac:
                    # Replace thing with the IUPAC code representing thing and the existing base if
                    # the overlapping reads disagree at this nucleotide.
                    self.col[self.offset] = iupac[(existing,thing)]
                else:
                    # Some sort of conflict between introns/indels/bases/Ns
                    self.col[self.offset] = "X"
                    self.conflicts += 1
        self.advance(1)
    def advance(self, n=1):
        self.offset = (self.offset + n) % self.length
    def log(self, message):
        log(f"@{self.offset:5}: {message}")
    def rewind(self, n):
        if n >= 0:
            raise RuntimeError("Rewinding requires negative offset")
        self.offset += n
        self.log(f"rewinding {n}")
        if self.offset < 0:
            self.offset += self.length
            self.log(f"rewinding {n} around to other side")

class Alignment:
    def __init__(self, splice):
        self.spliceLeft = splice['aIntronStart']
        self.spliceRight = splice['dIntronStart']
        if self.spliceLeft > self.spliceRight:
            (self.spliceLeft, self.spliceRight) = (self.spliceRight, self.spliceLeft)
        self.key = (self.spliceLeft,self.spliceRight)
        self.refStart = self.spliceLeft+1
        self.length = self.spliceRight-self.refStart
        self.chr = splice['chr']
        self.cols = []
        self.docs = []
        self.warnings = []
        self.excluded = 0
        self.repeatMasked = 0
        self.leftOffsetNegative = 0
        self.encompassingExtendsBeyond = 0
        self.addReference(splice['dna'])
        self.cellBarcodes = {}
        self.strand = None
    def baseFilename(self):
        return(f"{self.chr}-{self.key[0]}-{self.key[1]}.aln")
    def determineLeftAndRightSegOffsets(self, doc):
        if doc['strand'] == "+":
            return((doc['seg2Start'] - self.refStart, doc['seg1Start'] - self.refStart))
        return((doc['seg1Start'] - self.refStart, doc['seg2Start'] - self.refStart))
    def join(self, column):
        self.cols.append(column)
    def refAt(self, offset):
        return self.cols[0].col[offset]
    def addReference(self, ref):
        column = Column(self.length)
        self.ref = ref
        if len(self.ref) != self.length:
            raise RuntimeError(f"ref has {len(self.ref)} bp but expecting {self.length}")
        # Add the reference sequence
        for offset in range(self.length):
            base = self.ref[offset]
            upperBase = base.upper()
            column.place(upperBase)
            if base != upperBase:
                self.repeatMasked += 1
        self.join(column)
        self.compressibility = compressionRatio(self.ref)
    def newColumn(self, doc):
        column = Column(self.length)
        column.read = doc['read']
        column.strand = doc['strand']
        column.docId = doc['id']
        column.barcode = doc['barcode']
        column.sample = doc['sample']
        if not column.sample in self.cellBarcodes:
            self.cellBarcodes[column.sample] = {}
        if not column.barcode in self.cellBarcodes[column.sample]:
            self.cellBarcodes[column.sample][column.barcode] = 0
        self.cellBarcodes[column.sample][column.barcode] += 1
        column.addOtherStats(doc)
        return column
    def determineOrderedCigars(self, doc):
        if doc['strand'] == "+":
            return((2,doc['cigar2']),(1,doc['cigar1']))
        return((1,doc['cigar1']),(2,doc['cigar2']))
    def determineReadMates(self, doc):
        if doc['strand'] == "+":
            return((doc['read1'], revcomp(doc['read2'])))
        return((revcomp(doc['read1']), doc['read2']))
    def addEncompassingJunction(self, doc):
        log(f"processing line: {doc['line']}")
        strand = doc['strand']
        cigar1 = doc['cigar1']
        cigar2 = doc['cigar2']
        (read1, read2) = self.determineReadMates(doc)
        log(f"read1: {read1}")
        log(f"read2: {read2}")
        longestReadLen = max(len(read1), len(read2))
        if longestReadLen > self.length:
            log("warning: can't fit read within single cycle")
            self.excluded += 1
            return False
        if len([b for (b,n) in (cigar1+cigar2) if b == "I"]) > 0:
            # TODO: handle insertions
            log("warning: insertion...skipping")
            self.excluded += 1
            return False
        column = self.newColumn(doc)
        column.encompassing = True
        (leftSegOffset, rightSegOffset) = self.determineLeftAndRightSegOffsets(doc)
        log(f"read1 has length {len(read1)}")
        log(f"read2 has length {len(read2)}")
        log(f"ref is {self.length} bp")
        log(f"left offset is {leftSegOffset} and right offset is {rightSegOffset}")
        if leftSegOffset < 0:
            log(f"warning: left offset is negative ({leftSegOffset})")
            self.leftOffsetNegative += 1
            self.excluded += 1
            return False
        if rightSegOffset >= column.length:
            self.encompassingExtendsBeyond += 1
            self.excluded += 1
            return False
        seq = None
        otherSeq = None
        if strand == "+":
            seq = read2
            otherSeq = read1
        else:
            seq = read1
            otherSeq = read2
        seqOffset = 0
        # For an encompassing junction, we start part way through the alignment.
        column.advance(leftSegOffset)
        log(f"starting with {seq} at offset {seqOffset}")
        matches = 0
        misses = 0
        segment = 0
        orderedCigars = self.determineOrderedCigars(doc)
        firstCigar = True
        for c, cigar in orderedCigars:
            cigarDesc = f"cigar{c}: {cigar}"
            for i,(b,n) in enumerate(cigar):
                column.log(f"consuming {cigarDesc} tuple ({b},{n})")
                if b == "S":
                    self.ensureSoftClipAllowed(strand, i, len(cigar), c, True)
                    for j in range(n):
                        seqOffset += 1
                elif b == "M":
                    # Don't allow matches for encompassing reads that extend beyond the circular portion of the reference.
                    if column.offset + n > column.length:
                        log("warning: match beyond putative circular molecule...skipping")
                        self.encompassingExtendsBeyond += 1
                        self.excluded += 1
                        return False
                    for j in range(n):
                        base = seq[seqOffset]
                        if base == self.refAt(column.offset):
                            matches += 1
                        else:
                            misses += 1
                        column.place(base, segment)
                        seqOffset += 1
                elif b == "N":
                    # Don't allow gaps that extend beyond the circular portion of the reference.
                    if column.offset + n > column.length:
                        log("warning: gap beyond putative circular molecule...skipping")
                        self.excluded += 1
                        return False
                    column.placeLargeGap(n)
                elif b == "D":
                    # Deletion
                    if column.offset + n > column.length:
                        log("warning: deletion beyond putative circular molecule...skipping")
                        self.excluded += 1
                        return False
                    for j in range(n):
                        column.place("-")
                else:
                    raise RuntimeError(f"unsupported cigar code {b}")
            if firstCigar:
                firstCigar = False
                segment += 1
                # After we process the first cigar we get set up to process the other one
                column.log(f"moving {rightSegOffset-column.offset} to other segment")
                if column.offset > rightSegOffset:
                    # segment 1 must overlap with the previous chunk.
                    goBack = rightSegOffset - column.offset
                    if goBack < -20:
                        column.log(f"warning: rewinding {goBack} bp")
                    log(f"rewinding")
                    column.rewind(goBack)
                # Forward to segment 1
                log(f"advancing from {column.offset} to {rightSegOffset}")
                log(f"refStart: {self.refStart}, length: {self.length}")
                if column.offset < rightSegOffset:
                    column.advance(rightSegOffset - column.offset)
                # We need to switch reads.
                if seqOffset != len(seq):
                    column.log(f"should switch at end of sequence but at offset {seqOffset} of {len(seq)}")
                    raise ShouldBeAtEndOfSequence("should be switching reads and at end of sequence")
                log("swapping")
                (seq, otherSeq) = (otherSeq, seq)
                seqOffset = 0
                column.log("switching reads")
                # Now we switch to the other cigar string
                column.log(f"starting other segment at offset {seqOffset}")
        if seqOffset != len(seq):
            column.warn("expected to be at end of current seq chunk")
        column.matches = matches
        column.misses = misses
        self.join(column)
        log(f"encompassing read {doc['read']} aligned with {misses} mismatches")
        return True

    def addChimericJunction(self, doc):
        log(f"processing line: {doc['line']}")
        log(f"pFirst: {doc['pFirst']}")
        strand = doc['strand']
        cigar1 = doc['cigar1']
        cigar2 = doc['cigar2']
        (read1, read2) = self.determineReadMates(doc)
        log(f"read1: {read1}")
        log(f"read2: {read2}")
        longestReadLen = max(len(read1), len(read2))
        if longestReadLen > self.length:
            log("warning: can't fit read within single cycle")
            self.excluded += 1
            return
        if len([b for (b,n) in (cigar1+cigar2) if b == "I"]) > 0:
            # TODO: handle insertions
            log("warning: insertion...skipping")
            self.excluded += 1
            return
        column = self.newColumn(doc)
        (leftSegOffset, rightSegOffset) = self.determineLeftAndRightSegOffsets(doc)
        log(f"read1 has length {len(read1)}")
        log(f"read2 has length {len(read2)}")
        log(f"ref is {self.length} bp")
        # If the 'p' gap is in segment 1, then the broken read is read2, otherwise the
        # broken read is read1.
        seq = None
        otherSeq = None
        orderedCigars = self.determineOrderedCigars(doc)
        if doc['pFirst']:
            log("starting with read2")
            seq = read2
            otherSeq = read1
        else:
            log("starting with read1")
            seq = read1
            otherSeq = read2
        seqOffset = 0
        matches = 0
        misses = 0
        segment = 0
        firstCigar = True
        for c, cigar in orderedCigars:
            cigarDesc = f"cigar{c}: {cigar}"
            for i,(b,n) in enumerate(cigar):
                column.log(f"consuming {cigarDesc} tuple ({b},{n})")
                if b == "S":
                    self.ensureSoftClipAllowed(strand, i, len(cigar), c, False)
                    seqOffset += n
                elif b == "M":
                    for j in range(n):
                        base = seq[seqOffset]
                        if base == self.refAt(column.offset):
                            matches += 1
                        else:
                            misses += 1
                        column.place(base, segment)
                        seqOffset += 1
                elif b == "p":
                    if n < 0:
                        # In the case it's negative, we rewind and replace bases with IUPAC codes.
                        column.rewind(n)
                    else:
                        column.advance(n)
                    # When we hit p, we swap sequences
                    (seq, otherSeq) = (otherSeq, seq)
                    segment += 1
                    seqOffset = 0
                    column.log(f"switching reads at p in segment {c}")
                elif b == "N":
                    # Don't allow gaps that extend beyond the circular portion of the reference.
                    if column.offset + n > column.length:
                        log("warning: gap beyond putative circular molecule...skipping")
                        self.excluded += 1
                        return
                    column.placeLargeGap(n)
                elif b == "D":
                    # Deletion
                    for j in range(n):
                        column.place("-")
                else:
                    raise RuntimeError(f"unsupported cigar code {b}")
            if firstCigar:
                firstCigar = False
                segment += 1
                # After we process the first cigar, we get set up to process the other one
                column.log(f"moving {rightSegOffset-column.offset} to segment 1")
                if column.offset > rightSegOffset:
                    # this segment must overlap with the previous chunk.
                    goBack = rightSegOffset - column.offset
                    if goBack < -20:
                        column.log(f"warning: rewinding {goBack} bp")
                    log(f"rewinding")
                    column.rewind(goBack)
                # Forward to other segment
                log(f"advancing")
                if column.offset < rightSegOffset:
                    column.advance(rightSegOffset - column.offset)
                # We need to switch reads.
                if seqOffset != len(seq):
                    column.log(f"should switch at end of sequence but at offset {seqOffset} of {len(seq)}")
                    raise ShouldBeAtEndOfSequence("should be switching reads and at end of sequence")
                log("swapping")
                (seq, otherSeq) = (otherSeq, seq)
                seqOffset = 0
                if doc['pFirst']:
                    column.log("switching reads")
                else:
                    column.log("switching reads again")
                # Now we switch to the other cigar string
                column.log(f"starting other segment at offset {seqOffset}")
        if column.offset != 0:
            raise RuntimeError(f"expected to be at end of reference but at {column.offset}")
        if seqOffset != len(seq):
            column.warn("expected to be at end of current seq chunk")
        column.matches = matches
        column.misses = misses
        self.join(column)
        log(f"overlapping read {doc['read']} aligned with {misses} mismatches")
        log("joining column")

    def ensureSoftClipAllowed(self, strand, cigPos, cigLen, segment, isEncompassing):
        if isEncompassing:
            # TODO: not entirely sure this is correct
            if cigPos == 0 or cigPos+1 == cigLen:
                return
        if strand == "+":
            if segment == 2 and cigPos == 0:
                return
            if segment == 1 and cigPos+1 == cigLen:
                return
        else:
            if segment == 1 and cigPos == 0:
                return
            if segment == 2 and cigPos+1 == cigLen:
                return
        raise SoftClippingNotAllowed(f"soft-clip for strand {strand} in seg {segment} cig pos {cigPos} of {cigLen}")


parser = argparse.ArgumentParser(description='Assemble evidence of cirular RNAs.')
parser.add_argument('-overlapping', type=str, help='pickle containing overlapping junctions', required=True)
parser.add_argument('-encompassing', type=str, help='pickle containing encompassing junctions')
parser.add_argument('-focus', type=str, help='focus on just one alignment')
parser.add_argument('-downsample', type=int, help='If there are too may encompassign reads, limit to this number', default=500)
parser.add_argument('-max-compressibility', type=float, dest='maxZip', help='maximum compressiblity (otherwise skipped)')
parser.add_argument('-min-overlapping', type=int, dest='minOverlapping', help='minimum number of overlapping reads', default=1)
parser.add_argument('-min-encompassing', type=int, dest='minEncompassing', help='minimum number of encompassing reads', default=0)
parser.add_argument('-min-length', type=int, dest='minLen', help='minimum length of circular RNA (bp)', default=0)
parser.add_argument('-limit', type=int, help='limit number of alignments (default = 0, unlimited)', default=0)
args = parser.parse_args()
print(args)

log(f"loading pickle: {args.overlapping}")
overlapping = pickle.load(open(args.overlapping, "rb"))
encompassing = {}
if not args.encompassing is None:
    log(f"loading pickle: {args.encompassing}")
    encompassing = pickle.load(open(args.encompassing, "rb"))['docs']
docs = overlapping['docs']
backSplices = overlapping['backSplices']

pr = cProfile.Profile()
pr.enable()
if not args.focus is None:
    log(f"focusing on {args.focus}")
# Sort by ID.
tupleOffset = 0
softClippingErrors = 0
poorlyAlignedSequenceErrors = 0
alignments = 0
tuples = sorted(list(backSplices.items()), key=lambda x: x[0])
skipped = 0
tried = 0
log(f"number of tuples: {len(tuples)}")
for spliceId, splice in tuples:
    if args.limit > 0 and tried >= args.limit:
        log(f"hit limit of {args.limit} alignments tried")
        break
    if tried % 100 == 0 and tried > 0:
        keepPercent = int(alignments / tried * 1000) / 10
        log(f"tried {tried} candidates, skipped {skipped}, and aligned {alignments} ({keepPercent}%)")
    tried += 1
    log(f"processing splice {spliceId} ({tupleOffset})")
    tupleOffset += 1
    aln = Alignment(splice)
    if aln.length < args.minLen:
        log(f"circular RNA is too short ({aln.length}bp), skipping")
        skipped += 1
        continue
    if (not (args.maxZip is None)) and aln.compressibility > args.maxZip:
        log(f"skipping because compressibility {aln.compressibility} > {args.maxZip}")
        skipped += 1
        continue
    if not args.focus is None:
        if aln.baseFilename() == args.focus:
            log(f"found focus {args.focus}")
        else:
            log(f"not focusing on {aln.baseFilename()}...skipping")
            skipped += 1
            continue
    for docId in splice['docIds']:
        doc = docs[docId]
        log(f"processing doc {docId} {doc['sample']}/{doc['barcode']}")
        try:
            aln.addChimericJunction(doc)
        except SoftClippingNotAllowed as ex:
            softClippingErrors += 1
            log(f"unallowed softclipping for doc {docId}: {ex}")
            continue
        except ShouldBeAtEndOfSequence as ex:
            poorlyAlignedSequenceErrors += 1
            continue
    if len(aln.cols)-1 < args.minOverlapping:
        log(f"skipping alignment because overlapping reads, {len(aln.cols)-1} < {args.minOverlapping}")
        skipped += 1
        continue
    if len(aln.cols) > args.minOverlapping:
        # Look through the encompassing junctions and find ones that support this junction.
        toAdd = []
        for docId, doc in encompassing.items():
            if doc['chr'] != aln.chr:
                continue
            if doc['bound5Prime'] < aln.refStart:
                continue
            if doc['bound3Prime'] > aln.refStart+aln.length:
                continue
            toAdd.append(doc)
        log(f"found {len(toAdd)} encompassing reads for {spliceId}")
        if len(toAdd) < args.minEncompassing:
            log("too few encompassing reads")
            skipped += 1
            continue
        alignments += 1

log(f"encountered {softClippingErrors} soft clipping errors")
log(f"encountered {poorlyAlignedSequenceErrors} poorly aligned sequencing errors")
log(f"found {alignments} alignments")

log("done")
pr.disable()
pr.print_stats()


# END
