
import pysam
import re
import os
import argparse
import sys
import gzip
from collections import defaultdict

class JsRead(object):
    """For junction reads in alignment.bam"""
    def __init__(self, read, known_js):
        self.read = read
        self.js = self.get_js_from_cigartuples()
        self.left_block_li = list()
        self.right_block_li = list()
        self.left_sequence_li = list()
        self.right_sequence_li = list()
        self.mid_sequence_li = list()
        if read.reference_name in known_js.keys():
            self.new_js = [x for x in self.js if x not in known_js[read.reference_name]]
        else:
            self.new_js = self.js


    __slots__ = (
        "read", "js", "new_js", "left_block_li", "right_block_li",
        "left_sequence_li", "right_sequence_li", "mid_sequence_li"
    )


    def get_js_from_cigartuples(self):
        start = self.read.reference_start
        js_li = list()
        for cigar, length in self.read.cigartuples:
            if cigar == 0:
                start = start+length
            elif cigar == 3:
                js_li.append((start, start+length))
                start = start+length
            elif cigar in [4, 1]:
                continue
        block_js = [(self.read.blocks[i][1], self.read.blocks[i + 1][0]) for i in range(len(self.read.blocks) - 1)]
        for js in js_li:
            try:
                assert js in block_js
            except:
                print("Error:", self.read.cigarstring, self.read.qname)
                sys.exit(1) 
        return js_li

    @property
    def name(self):
        return self.read.query_name

    @property
    def is_read1(self):
        return self.read.is_read1

    @property
    def chrom(self):
        return self.read.reference_name

    @property
    def strand(self):
        if self.read.is_reverse:
            return "-"
        else:
            return "+"

    def ref2query(self, site):
        assert self.read.reference_start <= site <= self.read.reference_end
        if site == self.read.reference_end:
            return len(self.read.query_sequence)
        move_site = self.read.reference_start
        query_size = 0

        for cigar, length in self.read.cigar:
            assert cigar <= 6
            if move_site == site:
                continue
            if cigar == 0:
                if site-move_site <= length:  # ref site is in this block
                    query_size = query_size + (site-move_site)
                    move_site += (site-move_site)
                else:
                    move_site += length
                    query_size += length 
            elif cigar in [1, 4, 5]:  # I S H
                query_size += length
            elif cigar in [2, 3, 6]:  # D N P
                move_site += length
        return query_size
    
    def fetch_arm_info(self):
        for js_start, js_end in self.new_js:
            query_start = self.ref2query(js_start)
            query_end = self.ref2query(js_end)
            sequence = self.read.query_sequence.upper()
            mid_sequence = sequence[query_start:query_end]
            self.mid_sequence_li.append(mid_sequence)
            left_blocks = list()
            right_blocks = list()
            for s, e in self.read.blocks:
                if s < js_start:
                    left_blocks.append((s, min(e, js_start)))
                if e > js_end:
                    right_blocks.append((max(s, js_end), e))

            left_block = left_blocks[-1]
            right_block = right_blocks[0]

            left_sequence = sequence[self.ref2query(left_block[0]): self.ref2query(left_block[1])]
            right_sequence = sequence[self.ref2query(right_block[0]): self.ref2query(right_block[1])]

            self.left_block_li.append(left_block)
            self.right_block_li.append(right_block)
            self.left_sequence_li.append(left_sequence)
            self.right_sequence_li.append(right_sequence)


class BamFile(object):
    """Methods for alignment.bam"""
    def __init__(self, f_bam):
        self.bam = pysam.AlignmentFile(f_bam, "rb")

    def iter_new_js_reads(self, known_js, rule):
        for read in self.bam.fetch():
            if read.is_secondary:
                continue
            if read.is_qcfail:
                continue
            if read.is_unmapped:
                continue
            if read.mate_is_unmapped:
                continue
            if read.reference_name is None:
                continue
            if read.reference_start is None:
                continue
            if read.reference_end is None:
                continue
            if read.cigarstring.find("N") == -1:  # No junction
                continue
            nh_tag = read.get_tag("NH") if read.has_tag("NH") else None
            if nh_tag != 1:
                continue

            if rule == "1++,1--,2+-,2-+":
                if read.is_read2:
                    read.is_reverse = not read.is_reverse
                    read.mate_is_reverse = not read.mate_is_reverse
            elif rule == "1+-,1-+,2++,2--":
                if read.is_read1:
                    read.is_reverse = not read.is_reverse
                    read.mate_is_reverse = not read.mate_is_reverse
            elif rule == "+-,-+":
                read.is_reverse = not read.is_reverse
            js_read = JsRead(read, known_js)
            if not js_read.new_js:
                continue
            js_read.fetch_arm_info()
            yield js_read

class JsReadPair(object):
    """For junction pairs in alignment.bam"""

    def __init__(self, name):
        self.name = name
        self.read1 = None
        self.read2 = None
        self.new_js = None
        self.chrom = None
        self.strand = None
        self.donor_blocks =list()
        self.acceptor_blocks = list()
        self.donor_sequence = list()
        self.acceptor_sequence = list()
        self.mid_sequence = list()
        self.is_conflict = None
        self.conflict_info = None

    __slots__ = (
        "name", "read1", "read2", "new_js", "chrom", "strand",
        "donor_blocks", "acceptor_blocks",
        "donor_sequence", "acceptor_sequence", "mid_sequence",
        "is_conflict", "conflict_info"
    )

    def add_read(self, js_read):
        if js_read.is_read1:
            assert self.read1 is None
            self.read1 = js_read
        else:
            assert self.read2 is None
            self.read2 = js_read

    def arm_sequence_diff(self, read1_squence, read2_sequence, pos):
        read1_sequence = list(read1_squence)
        read2_sequence = list(read2_sequence)
        length = min(len(read1_sequence), len(read2_sequence))
        diff_base = 0
        if pos == "left":
            for j in range(length):
                if read1_sequence[-(j+1)] != read2_sequence[-(j+1)]:
                    diff_base += 1
        else:
            for j in range(length):
                if read1_sequence[j] != read2_sequence[j]:
                    diff_base += 1
        return diff_base/length

    def merge_js_info(self):
        
        if self.read1 is None:
            self.new_js = self.read2.new_js
            self.chrom = self.read2.chrom
            self.strand = self.read2.strand
            self.donor_blocks = self.read2.left_block_li
            self.acceptor_blocks = self.read2.right_block_li
            self.donor_sequence = self.read2.left_sequence_li
            self.acceptor_sequence = self.read2.right_sequence_li
            self.mid_sequence = self.read2.mid_sequence_li
        elif self.read2 is None:
            self.new_js = self.read1.new_js
            self.chrom = self.read1.chrom
            self.strand = self.read1.strand
            self.donor_blocks = self.read1.left_block_li
            self.acceptor_blocks = self.read1.right_block_li
            self.donor_sequence = self.read1.left_sequence_li
            self.acceptor_sequence = self.read1.right_sequence_li
            self.mid_sequence = self.read1.mid_sequence_li
        else:
            new_js = list(set(self.read1.new_js + self.read2.new_js))
            self.new_js = new_js
            self.chrom = self.read1.chrom
            self.strand = self.read1.strand
            for i in range(len(new_js)):
                temp_js = new_js[i]
                if temp_js in self.read1.new_js and temp_js in self.read2.new_js:
                    read1_indx = [indx for indx, block in enumerate(self.read1.new_js) if block == temp_js][0]
                    read2_index = [indx for indx, block in enumerate(self.read2.new_js) if block == temp_js][0]
                    left_arm_diff = self.arm_sequence_diff(self.read1.left_sequence_li[read1_indx], self.read2.left_sequence_li[read2_index], "left")
                    right_arm_diff = self.arm_sequence_diff(self.read1.right_sequence_li[read1_indx], self.read2.right_sequence_li[read2_index], "right")
                    if max(left_arm_diff, right_arm_diff) > 0.25:
                        self.conflict_info = "Different arm sequence between R1 and R2"
                    self.donor_sequence.append(self.read1.left_sequence_li[read1_indx])
                    self.acceptor_sequence.append(self.read1.right_sequence_li[read1_indx])
                    self.donor_blocks.append(self.read1.left_block_li[read1_indx])
                    self.acceptor_blocks.append(self.read1.right_block_li[read1_indx])
                    self.mid_sequence.append(self.read1.mid_sequence_li[read1_indx])
                elif temp_js in self.read1.new_js and temp_js not in self.read2.new_js:
                    read1_indx = [indx for indx, block in enumerate(self.read1.new_js) if block == temp_js][0]
                    self.donor_sequence.append(self.read1.left_sequence_li[read1_indx])
                    self.acceptor_sequence.append(self.read1.right_sequence_li[read1_indx])
                    self.donor_blocks.append(self.read1.left_block_li[read1_indx])
                    self.acceptor_blocks.append(self.read1.right_block_li[read1_indx])
                    self.mid_sequence.append(self.read1.mid_sequence_li[read1_indx])
                elif temp_js in self.read2.new_js and temp_js not in self.read1.new_js:
                    read2_indx = [indx for indx, block in enumerate(self.read2.new_js) if block == temp_js][0]
                    self.donor_sequence.append(self.read2.left_sequence_li[read2_indx])
                    self.acceptor_sequence.append(self.read2.right_sequence_li[read2_indx])
                    self.donor_blocks.append(self.read2.left_block_li[read2_indx])
                    self.acceptor_blocks.append(self.read2.right_block_li[read2_indx])
                    self.mid_sequence.append(self.read2.mid_sequence_li[read2_indx])


    def write(self, fo):
        if self.read1 is not None and self.read2 is not None:
            assert self.read1.strand == self.read2.strand 
        for i in range(len(self.new_js)):
            js_start, js_end = self.new_js[i]
            donor_block = "%s-%s" % (self.donor_blocks[i][0], self.donor_blocks[i][1])
            acceptor_block = "%s-%s" % (self.acceptor_blocks[i][0], self.acceptor_blocks[i][1]) 
            donor_sequence = self.donor_sequence[i]
            acceptor_sequence = self.acceptor_sequence[i]
            mid_sequence = self.mid_sequence[i]
            if self.strand == "+":
                outline = [self.chrom, js_start, self.strand, donor_block, donor_sequence, self.chrom, js_end, self.strand, acceptor_block, acceptor_sequence, mid_sequence, self.name, "STAR.Aligned"]
            else:
                outline = [self.chrom, js_end, self.strand, acceptor_block, acceptor_sequence, self.chrom, js_start, self.strand, donor_block, donor_sequence, mid_sequence, self.name, "STAR.Aligned"]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")

def load_known_js(f_known_js):
    known_JS_dict = defaultdict(set)
    if f_known_js is not None:
        print("Loading known JS ...")
        with open(f_known_js) as fi:
            for line in fi:
                data = line.rstrip("\n").split("\t")
                chrom = data[0]
                start = int(data[1])
                end = int(data[2])
                known_JS_dict[chrom].add((start, end))
            all_js_num = sum([len(x) for x in known_JS_dict.values()])
        print("{0} known JS were loaded.".format(all_js_num))
    return known_JS_dict

class ChimericJS(object):
    def __init__(self, name):
        self.name = name
        self.read1_li = list()
        self.read2_li = list()
        self.read1_donor = None
        self.read1_acceptor = None
        self.read2_donor = None
        self.read2_acceptor = None
        self.read1_fastq = None
        self.read2_fastq = None
        self.rule = None
        self.mid_sequence = None

    def add_read(self, read, rule):
        self.rule = rule
        if read.is_read1:
            if rule == "1+-,1-+,2++,2--":
                read.is_reverse = not read.is_reverse
                read.mate_is_reverse = not read.mate_is_reverse
            self.read1_li.append(read)
        if read.is_read2:
            if rule == "1++,1--,2+-,2-+":
                read.is_reverse = not read.is_reverse
                read.mate_is_reverse = not read.mate_is_reverse
            self.read2_li.append(read)

    def get_match_sequence(self, read):
        abs_start = 0
        abs_end = 0
        start = 0
        cigar0_index = 0
        for cigar, length in read.cigar:
            if cigar in [4, 2]: # S
                start += length
            elif cigar in [1]: # I
                start += length
                abs_end += length
            elif cigar in [3, 4, 5]:
                continue
            elif cigar == 0:
                if cigar0_index == 0:
                    abs_start = start
                    abs_end = start
                cigar0_index += 1
                start += length
                abs_end += length
        return read.query_sequence[abs_start:abs_end]

    def fetch_chimeric_from_read(self, read, known_js):
        js_read = JsRead(read, known_js)
        if js_read.new_js:
            if len(js_read.new_js) != 1:
                return None
            assert len(js_read.new_js) == 1
            new_js = js_read.new_js[0]
            js_read.fetch_arm_info()

            if self.rule == "1+-,1-+,2++,2--":
                if read.is_read1:
                    self.read1_donor = (js_read.chrom, new_js[1], js_read.strand, "-".join(list(map(str, js_read.right_block_li[0]))), js_read.right_sequence_li[0])
                    self.read1_acceptor = (js_read.chrom, new_js[0], js_read.strand, "-".join(list(map(str, js_read.left_block_li[0]))), js_read.left_sequence_li[0])
                    self.read1_mid_sequence = js_read.mid_sequence_li[0]
                if read.is_read2:
                    self.read2_donor = (js_read.chrom, new_js[0], js_read.strand, "-".join(list(map(str, js_read.left_block_li[0]))), js_read.left_sequence_li[0])
                    self.read2_acceptor = (js_read.chrom, new_js[1], js_read.strand, "-".join(list(map(str, js_read.right_block_li[0]))), js_read.right_sequence_li[0])
                    self.read2_mid_sequence = js_read.mid_sequence_li[0]
            elif self.rule == "1++,1--,2+-,2-+":
                if read.is_read1:
                    self.read1_donor = (js_read.chrom, new_js[0], js_read.strand, "-".join(list(map(str, js_read.left_block_li[0]))), js_read.left_sequence_li[0])
                    self.read1_acceptor = (js_read.chrom, new_js[1], js_read.strand,"-".join(list(map(str, js_read.right_block_li[0]))), js_read.right_sequence_li[0])
                    self.read1_mid_sequence = js_read.mid_sequence_li[0]
                if read.is_read2:
                    self.read2_donor = (js_read.chrom, new_js[1], js_read.strand, "-".join(list(map(str, js_read.right_block_li[0]))), js_read.right_sequence_li[0])
                    self.read2_acceptor = (js_read.chrom, new_js[0], js_read.strand, "-".join(list(map(str, js_read.left_block_li[0]))), js_read.left_sequence_li[0])
                    self.read2_mid_sequence = js_read.mid_sequence_li[0]

    def get_reverse_seq(self, sequence):
        base_dict = {'n':'N', "a":"T", "t":"A", "c":"G", "g":"C"}
        seq = list(sequence.lower())[::-1]
        rev_seq = [base_dict[base] for base in seq]
        return "".join(rev_seq)

    def get_sequence_relative_pos(self, read, seq):
        if read.is_reverse:
            seq = self.get_reverse_seq(seq)
        if read.is_read1:
            for x in re.finditer(seq, self.read1_fastq):
                start, end = x.span()
                return start+1, end
        else:
            for x in re.finditer(seq, self.read2_fastq):
                start, end = x.span()
                return start+1, end

    def fetch_chimeric_from_read_li(self, read_li):
        part1 = read_li[0]
        part2 = read_li[1]
        part1_strand = "-" if part1.is_reverse else "+"
        part2_strand = "-" if part2.is_reverse else "+"
        part1_block = ["%s-%s" % (b[0], b[1]) for b in part1.get_blocks()]
        part2_block = ["%s-%s" % (b[0], b[1]) for b in part2.get_blocks()]
        part1_sequence = self.get_match_sequence(part1)
        part2_sequence = self.get_match_sequence(part2)

        part1_start, part1_end = self.get_sequence_relative_pos(part1, part1_sequence)
        part2_start, part2_end = self.get_sequence_relative_pos(part2, part2_sequence)
        
        if part1_end < part2_start:
            left_chrom, left_strand, left_block, left_sequence = part1.reference_name, part1_strand, ";".join(part1_block), part1_sequence
            left_site = part1.reference_end if part1_strand == "+" else part1.reference_start
            right_chrom, right_strand, right_block, right_sequence = part2.reference_name, part2_strand, ";".join(part2_block), part2_sequence
            right_site = part2.reference_start if part2_strand == "+" else part2.reference_end
            if part1.is_read1:
                self.read1_mid_sequence = self.read1_fastq[part1_end:(part2_start-1)]
            else:
                self.read2_mid_sequence = self.read2_fastq[part1_end:(part2_start-1)]
        else:
            left_chrom, left_strand, left_block, left_sequence = part2.reference_name, part2_strand, ";".join(part2_block), part2_sequence
            left_site = part2.reference_end if part2_strand == "+" else part2.reference_start
            right_chrom, right_strand, right_block, right_sequence = part1.reference_name, part1_strand, ";".join(part1_block), part1_sequence
            right_site = part1.reference_start if part1_strand == "+" else part1.reference_end
            if part1.is_read1:
                self.read1_mid_sequence = self.read1_fastq[part2_end:(part1_start-1)]
            else:
                self.read2_mid_sequence = self.read2_fastq[part2_end:(part1_start-1)]

        if part1.is_read1:
            self.read1_donor = (left_chrom, left_site, left_strand, left_block, left_sequence)
            self.read1_acceptor = (right_chrom, right_site, right_strand, right_block, right_sequence)
        if part1.is_read2:
            self.read2_donor = (left_chrom, left_site, left_strand, left_block, left_sequence)
            self.read2_acceptor = (right_chrom, right_site, right_strand, right_block, right_sequence)
    
    def write(self, fo):
        if self.read1_donor is not None:
            outline = list(self.read1_donor) + list(self.read1_acceptor) + [self.read1_mid_sequence, self.name, "STAR.chimeric"]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")
        if self.read2_donor is not None:
            outline = list(self.read2_donor) + list(self.read2_acceptor) + [self.read2_mid_sequence, self.name, "STAR.chimeric"]
            outline = list(map(str, outline))
            fo.write("\t".join(outline)+"\n")

def is_overlap(interval1, interval2):
    x1, y1 = interval1
    x2, y2 = interval2
    if y1 <= x2 or y2 <= x1:
        return False
    else:
        return True

def block_overlap(target_blocks, source_blocks):
    target_li = list()
    for target_blcok in target_blocks:
        temp_li = list()
        for source_block in source_blocks:
            if is_overlap(target_blcok, source_block):
                temp_li.append(True)
        if temp_li:
            target_li.append(True)
    if target_li:
        return True
    else:
        return False

def load(f_in):
    with gzip.open(f_in, "rb") as f_i:
        array = []
        for i, line_data in enumerate(f_i):
            line = line_data.rstrip().decode()
            array.append(line)
            if i % 4 == 3:
                yield array
                array = []

def load_chimeric_fastq(f_fastq1, f_fastq2, chimeric_js_dict, rule):
    loader1 = load(f_fastq1)
    loader2 = load(f_fastq2)
    while True:
        try:
            name1, seq1, strand1, qua1 = next(loader1)
            name2, seq2, strand2, qua2 = next(loader2)
            name1 = name1.strip().split(" ")[0].split("@")[1]
            name2 = name2.strip().split(" ")[0].split("@")[1]
            assert name1 == name2
            if name1 in chimeric_js_dict:
                if rule == "1+-,1-+,2++,2--":
                    chimeric_js_dict[name1].read1_fastq = chimeric_js_dict[name1].get_reverse_seq(seq1)
                    chimeric_js_dict[name1].read2_fastq = seq2
                elif rule == "1++,1--,2+-,2-+":
                    chimeric_js_dict[name1].read1_fastq = seq1
                    chimeric_js_dict[name1].read2_fastq = chimeric_js_dict[name1].get_reverse_seq(seq2)
        except StopIteration:
            break
    return chimeric_js_dict

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    
    base_group.add_argument("-b", "--bam", type=str, dest="Bam", metavar="Aligned.sort.uniq.bam", required=True)
    base_group.add_argument("-c", "--chimeric", type=str, dest="Chimeric", metavar="Chimeric.out.sam", required=True)
    base_group.add_argument("-R1", "--read1", type=str, dest="read1", metavar="R1.fastq.gz", required=True)
    base_group.add_argument("-R2", "--read2", type=str, dest="read2", metavar="R2.fastq.gz", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.txt", required=True)
    base_group.add_argument("--known-js", type=str, dest="bed", metavar="ref.intron.bed", required=True)
    base_group.add_argument("--rule", type=str, dest="rule", metavar="1+-,1-+,2++,2--", required=True)

    return(parser.parse_args(args))


def main():
    args = sys.argv[1:]
    args = parse_args(args)

    f_bam = args.Bam
    f_chimeric = args.Chimeric
    f_R1 = args.read1
    f_R2 = args.read2
    f_output = args.output
    f_known_js = args.bed
    rule = args.rule


    known_JS_dict = load_known_js(f_known_js)
    header = ["DonorChrom", "DonorCJS", "DonorStrand", "DonorBlock", "DonorSequence", "AcceptorChrom", "AcceptorCJS", "AcceptorStrand", "AcceptorBlock", "AcceptorSequence", "MidSequence", "ReadName", "Source"]

    print("Loading alignment file ...")
    js_read_dict = dict()
    bam = BamFile(f_bam)
    for indx, js_read in enumerate(bam.iter_new_js_reads(known_JS_dict, rule)):
        if indx % 10000 == 0:
            print("{0} novel junction reads were loaded.".format(indx))
        name = js_read.name
        if name not in js_read_dict.keys():
            js_read_dict[name] = JsReadPair(name)
        js_read_dict[name].add_read(js_read)

    print("Loading chimeric file ...")
    ChimericSam = pysam.AlignmentFile(f_chimeric, "r")
    chimeric_read_dict = dict()
    for read in ChimericSam:
        if read.qname not in chimeric_read_dict:
            chimeric_read_dict[read.qname] = ChimericJS(read.qname)
        chimeric_read_dict[read.qname].add_read(read, rule)

    chimeric_js_dict = load_chimeric_fastq(f_R1, f_R2, chimeric_read_dict, rule)
    for read_name in chimeric_js_dict:
        read_info = chimeric_js_dict[read_name]
        if len(read_info.read1_li) == 1 and len(read_info.read2_li) == 2:
            read1 = read_info.read1_li[0]
            part1 = read_info.read2_li[0]
            part2 = read_info.read2_li[1]
            read1_block = read1.get_blocks()
            part1_blcok = part1.get_blocks()
            part2_blcok = part2.get_blocks()
            if block_overlap(read1_block, part1_blcok) or block_overlap(read1_block, part2_blcok):
                read_info.fetch_chimeric_from_read_li(read_info.read2_li)
            else:
                read_info.fetch_chimeric_from_read(read_info.read1_li[0], known_JS_dict)
                read_info.fetch_chimeric_from_read_li(read_info.read2_li)
        elif len(read_info.read1_li) == 2 and len(read_info.read2_li) == 1:
            read1 = read_info.read2_li[0]
            part1 = read_info.read1_li[0]
            part2 = read_info.read1_li[1]
            read1_block = read1.get_blocks()
            part1_blcok = part1.get_blocks()
            part2_blcok = part2.get_blocks()
            if block_overlap(read1_block, part1_blcok) or block_overlap(read1_block, part2_blcok):
                read_info.fetch_chimeric_from_read_li(read_info.read1_li)
            else:
                read_info.fetch_chimeric_from_read(read_info.read2_li[0], known_JS_dict)
                read_info.fetch_chimeric_from_read_li(read_info.read1_li)
        elif len(read_info.read1_li) == 2 and len(read_info.read2_li) == 2:
            read_info.fetch_chimeric_from_read(read_info.read1_li)
            read_info.fetch_chimeric_from_read(read_info.read2_li)
        elif len(read_info.read1_li) == 1 and len(read_info.read2_li) == 1:
            read_info.fetch_chimeric_from_read(read_info.read1_li[0], known_JS_dict)
            read_info.fetch_chimeric_from_read(read_info.read2_li[0], known_JS_dict)
        else:
            print("error")

    out = open(f_output, "w")
    out.write("\t".join(header)+"\n")
    for read_name in js_read_dict:
        js_read = js_read_dict[read_name]
        js_read.merge_js_info()
        js_read.write(out)
    for read_name in chimeric_js_dict:
        read_info = chimeric_js_dict[read_name]
        read_info.write(out)
    out.close()
if __name__ == "__main__":
    main()
