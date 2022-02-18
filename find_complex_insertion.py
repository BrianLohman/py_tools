#!/usr/bin/env python
# USAGE: find_complex_insertions.py -r [reference fasta] -b [BAM format file (uncompressed)] -o [output file name]
# for use when thinned reads are mapped back to plasmid only reference
# using ONLY the PRIMARY alignment
# using ONLY read with bit FLAG of 0 or 16 (forward or reverse)
# collects the 20 bases upstream and the 20 bases downstream of the insertion as a single string
# ./find_complex_insertions.py -r pk184.fa -s sorted_046_3LTR_pk184.sam -o first_pass.txt
# load modules
import argparse
import sys
from collections import Counter
from inspect import currentframe
from inspect import getframeinfo

from pyfaidx import Fasta

# setup argparse
parser = argparse.ArgumentParser(
    description="Get reference sequence adjacent to reads in filtered BAM"
)
parser.add_argument("-r", "--reference", dest="ref", help="reference fasta")
parser.add_argument("-s", "--SAM", dest="sam", help="filtered SAM file (non-binary")
parser.add_argument("-o", "--output_file", dest="out", help="name of output file")
args = parser.parse_args()

# name input and output
ref = Fasta(args.ref)
sam_file = args.sam
forward_out = open("forward_" + args.out, "w")
reverse_out = open("reverse_" + args.out, "w")
fail_cigars = open("fail_cigars_" + args.out, "w")
log = open(args.out[:-4] + ".log", "w")
fail_log = []

# print the headers for the output file
print(
    "Consensus",
    "Insertion",
    "5_flanking",
    "3_flanking",
    "CIGAR",
    "POS",
    sep="\t",
    file=forward_out,
)
print(
    "Consensus",
    "Insertion",
    "5_flanking",
    "3_flanking",
    "CIGAR",
    "POS",
    sep="\t",
    file=reverse_out,
)

# counters
k = 0
u = 0
sa = 0
sarc = 0
no_ref = 0
m = 0
l = 0
h = 0
complex_cigar = 0
f = 0
test = 0

with open(sam_file, "r") as sam:
    for toks in (line.strip().split("\t") for line in sam):
        # progress counter
        k += 1
        if k % 10000 == 0:
            print(str(k) + " reads processed")

        # skip headers
        if toks[0] in ["@HD", "@SQ", "@PG"]:
            continue

        # use only reads with a flag for foward or reverseq
        if toks[1] in ["0", "16"]:
            pass
        else:
            if toks[1] == "4":
                u += 1
                continue
            if toks[1] == "2048":
                sa += 1
                continue
            if toks[1] == "2064":
                sarc += 1
                continue

        # skip the reads that don't map to the reference (* in RNAME)
        if toks[2] == "*":
            no_ref += 1
            continue

        # set empty variables for each iteration of the loop
        flanking = None
        begin = None
        end = None
        flanking3 = None
        flanking5 = None
        flanking_middle = None
        tail = None
        left = None
        insertion = None
        flanking_looped = None

        # start with simple cigar strings, these are shorter then complex ones
        if toks[5].count("M") == 1 and toks[5].count("S") == 1:

            # FORWARD insertions (mapped region in the right portion of read means HIV LTR in left)
            if toks[5][-1] == "M":

                # find the reference sequence after the LTR
                if int(toks[3]) > 20:
                    m += 1
                    begin = (
                        int(toks[3]) - 1
                    )  # the first base to pull from the reference (SeqIO is 0 based index)
                    end = begin + 20  # the last base to pull from the reference
                    left = begin - 20  # the leftmost base to pull from the 5` side
                    insertion = int(toks[3])  # the insertion site
                    flanking3 = str(ref["verified_pk184"][begin:end])
                    flanking5 = str(ref["verified_pk184"][left:begin])
                    flanking = str(flanking5 + flanking3)

                else:
                    l += 1
                    insertion = int(toks[3])  # the insertion site
                    begin = (
                        int(toks[3]) - 1
                    )  # the first base to pull from the reference (SeqIO is 0 based index)
                    end = begin + 20  # the last base to pull from the reference
                    tail = (20 - begin) * -1
                    flanking3 = str(ref["verified_pk184"][0:end])
                    flanking5 = str(ref["verified_pk184"][tail:])
                    flanking = str(flanking5 + flanking3)

                # check that flanking has been asigned and it is the correct length
                if flanking == None:
                    sys.exit("flanking region not assigned")

                if len(flanking) != 40:
                    print(
                        "FORWARD",
                        flanking,
                        str(insertion),
                        flanking5,
                        flanking3,
                        str(tail),
                        str(begin),
                        str(end),
                    )
                    sys.exit("flanking not correct length")

                print(
                    flanking,
                    str(insertion),
                    flanking5,
                    flanking3,
                    toks[5],
                    toks[3],
                    sep="\t",
                    file=forward_out,
                )
                continue

            # REVERSE insertions (mapped region in the left portion of the read means HIV LTR in the right)
            elif toks[5][-1] == "S":

                # set the insertion site as the last base of the matched region + 1
                insertion = int(toks[3]) + int(toks[5].split("M")[0])

                # find the reference sequence after the LTR
                if insertion > 20 and insertion <= 2405:
                    m += 1
                    begin = (
                        insertion - 1
                    )  # the first base to pull from the reference (SeqIO is 0 based index)
                    end = begin + 20  # the last base to pull from the reference
                    left = begin - 20  # the leftmost base to pull from the 5` side
                    flanking3 = str(ref["verified_pk184"][begin:end])
                    flanking5 = str(ref["verified_pk184"][left:begin])
                    flanking = str(flanking5 + flanking3)

                elif insertion >= 1 and insertion <= 20:
                    l += 1
                    begin = (
                        insertion - 1
                    )  # the first base to pull from the reference (SeqIO is 0 based index)
                    end = begin + 20  # the last base to pull from the reference
                    tail = (20 - begin) * -1
                    flanking3 = str(ref["verified_pk184"][0:end])
                    flanking5 = str(ref["verified_pk184"][tail:])
                    flanking = str(flanking5 + flanking3)

                elif insertion >= 2406:
                    h += 1
                    if (
                        insertion == 2425
                    ):  # special case: the insertion is exactly between the end and the begining of pk184
                        flanking5 = str(ref["verified_pk184"][2405:2425])
                        flanking3 = str(ref["verified_pk184"][0:20])
                        flanking = str(flanking5 + flanking3)

                    else:
                        begin = (
                            insertion - 1
                        )  # the first base to pull from the reference (SeqIO is 0 based index)
                        left = begin - 20
                        flanking3 = str(ref["verified_pk184"][begin:])
                        flanking_looped = str(
                            ref["verified_pk184"][: 20 - len(flanking3)]
                        )
                        flanking5 = str(ref["verified_pk184"][left:begin])
                        flanking = str(flanking5 + flanking3 + flanking_looped)

            else:
                sys.exit("error handling insertion point in reverse oriention")

            # check that flanking has been asigned and it is the correct length
            if flanking == None:
                sys.exit("flanking region not assigned")

            if len(flanking) != 40:
                print(
                    "REVERSE",
                    flanking,
                    str(insertion),
                    flanking5,
                    flanking3,
                    flanking_looped,
                    str(begin),
                    str(left),
                )
                print(toks)
                sys.exit("flanking not correct length")

            print(
                flanking,
                str(insertion),
                flanking5,
                flanking3,
                toks[5],
                toks[3],
                sep="\t",
                file=reverse_out,
            )
            continue

        else:
            complex_cigar += 1

            # FORWARD insertions: the mapped portion of the read is on the right
            # the first part of the cigar string a soft clipeed region and the ajacent match is > 10
            if toks[5][2] == "S" or toks[5][3] == "S":

                # require the leading soft clipped region is at least 20 bp
                if int(toks[5].split("S")[0]) < 20:
                    f += 1
                    fail_log.append(getframeinfo(currentframe()).lineno)
                    print(toks[5], file=fail_cigars)
                    continue

                # require the match be greater than 10 (the 3rd or 4th character will be M)
                elif toks[5].split("S")[1][2] == "M" or toks[5].split("S")[1][3] == "M":
                    if int(toks[3]) > 20:
                        m += 1
                        begin = (
                            int(toks[3]) - 1
                        )  # the first base to pull from the reference (SeqIO is 0 based index)
                        end = begin + 20  # the last base to pull from the reference
                        left = begin - 20  # the leftmost base to pull from the 5` side
                        insertion = int(toks[3])  # the insertion site
                        flanking3 = str(ref["verified_pk184"][begin:end])
                        flanking5 = str(ref["verified_pk184"][left:begin])
                        flanking = str(flanking5 + flanking3)

                    else:
                        l += 1
                        insertion = int(toks[3])  # the insertion site
                        begin = (
                            int(toks[3]) - 1
                        )  # the first base to pull from the reference (SeqIO is 0 based index)
                        end = begin + 20  # the last base to pull from the reference
                        tail = (20 - begin) * -1
                        flanking3 = str(ref["verified_pk184"][0:end])
                        flanking5 = str(ref["verified_pk184"][tail:])
                        flanking = str(flanking5 + flanking3)

                    # check that flanking has been asigned and it is the correct length
                    if flanking == None:
                        sys.exit("flanking region not assigned")

                    if len(flanking) != 40:
                        print(
                            "FORWARD",
                            flanking,
                            str(insertion),
                            flanking5,
                            flanking3,
                            str(tail),
                            str(begin),
                            str(end),
                        )
                        sys.exit("flanking not correct length")

                    print(
                        flanking,
                        str(insertion),
                        flanking5,
                        flanking3,
                        toks[5],
                        toks[3],
                        sep="\t",
                        file=forward_out,
                    )
                    continue

                # 100+ soft clipped bases at the start and a match at the end
                elif toks[5][-1] == "M" and len(toks[5].split("S")[0]) == 3:

                    if int(toks[3]) > 20:
                        m += 1
                        begin = (
                            int(toks[3]) - 1
                        )  # the first base to pull from the reference (SeqIO is 0 based index)
                        end = begin + 20  # the last base to pull from the reference
                        left = begin - 20  # the leftmost base to pull from the 5` side
                        insertion = int(toks[3])  # the insertion site
                        flanking3 = str(ref["verified_pk184"][begin:end])
                        flanking5 = str(ref["verified_pk184"][left:begin])
                        flanking = str(flanking5 + flanking3)

                    else:
                        l += 1
                        insertion = int(toks[3])  # the insertion site
                        begin = (
                            int(toks[3]) - 1
                        )  # the first base to pull from the reference (SeqIO is 0 based index)
                        end = begin + 20  # the last base to pull from the reference
                        tail = (20 - begin) * -1
                        flanking3 = str(ref["verified_pk184"][0:end])
                        flanking5 = str(ref["verified_pk184"][tail:])
                        flanking = str(flanking5 + flanking3)

                    # check that flanking has been asigned and it is the correct length
                    if flanking == None:
                        sys.exit("flanking region not assigned")

                    if len(flanking) != 40:
                        print(
                            "FORWARD",
                            flanking,
                            str(insertion),
                            flanking5,
                            flanking3,
                            str(tail),
                            str(begin),
                            str(end),
                        )
                        sys.exit("flanking not correct length")

                    print(
                        flanking,
                        str(insertion),
                        flanking5,
                        flanking3,
                        toks[5],
                        toks[3],
                        sep="\t",
                        file=forward_out,
                    )
                    continue

                # CIGARs that do not fall into any of the above categories
                # CIGARs with a leading soft clip of 9 or less will end up here
                else:
                    f += 1
                    fail_log.append(getframeinfo(currentframe()).lineno)
                    print(toks[5], file=fail_cigars)
                    continue

            # REVERSE insertions
            elif toks[5][-1] == "S":

                # leading and triling soft clip filters
                # require the tailing soft clip to be greater than 20
                # these statements catch M, I, and D, in the last 3 characters before the terminal S
                if "M" in toks[5][-4:-1]:
                    if int(toks[5][-4:-1].split("M")[-1]) < 20:
                        print(toks[5], file=fail_cigars)
                        fail_log.append(getframeinfo(currentframe()).lineno)
                        continue
                    else:
                        pass

                elif "I" in toks[5][-4:-1]:
                    if int(toks[5][-4:-1].split("I")[-1]) < 20:
                        print(toks[5], file=fail_cigars)
                        fail_log.append(getframeinfo(currentframe()).lineno)
                        continue
                    else:
                        pass

                elif "D" in toks[5][-4:-1]:
                    if int(toks[5][-4:-1].split("D")[-1]) < 20:
                        print(toks[5], file=fail_cigars)
                        fail_log.append(getframeinfo(currentframe()).lineno)
                        continue
                    else:
                        pass

                # if none of the above are present then the terminal soft clip is over 100 bp
                else:
                    if int(toks[5][-4:-1]) < 20:
                        print(toks[5], file=fail_cigars)
                        fail_log.append(getframeinfo(currentframe()).lineno)
                        continue
                    else:
                        pass

                # when there is a leading soft clip, remove reads with >20S
                if toks[5].count("S") > 1:
                    if int(toks[5].split("S")[0]) > 20:
                        print(toks[5], file=fail_cigars)
                        fail_log.append(getframeinfo(currentframe()).lineno)
                        continue
                    else:
                        pass

                # start with POS and adjust for sum of matches and deltions
                insertion = int(toks[3]) + 1  # start with the first mapped maching base

                # get sum of matches
                ma = toks[5].split("M")
                for i, x in enumerate(ma):
                    # possible first (9 or less leading matc)
                    if i == 0 and x.isdigit():
                        insertion = insertion + int(ma[0])
                        continue
                    # last (all cigars should end with S)
                    if x[-1] == "S":
                        continue
                    # when the cigar starts with a soft clip
                    if x[1] == "S" or x[2] == "S":
                        insertion = insertion + int(x.split("S")[-1])
                        continue
                    # the middle of the cigar
                    if "I" in x:
                        insertion = insertion + int(x.split("I")[-1])
                        continue
                    if "D" in x:
                        insertion = insertion + int(x.split("D")[-1])
                        continue

                # if deletions are present add to insertion
                if "D" in toks[5]:
                    d = toks[5].split("D")
                    for i, c in enumerate(d):
                        # the last element in the list (no deletion info present)
                        # print(i, c, insertion)
                        if i + 1 == len(d):
                            continue
                        # the middle of the cigar
                        if c[-2:].isdigit():
                            insertion = insertion + int(c[-2:])
                        else:
                            insertion = insertion + int(c[-1])

                # find the reference sequence after the LTR
                if insertion > 20 and insertion <= 2405:
                    m += 1
                    begin = (
                        insertion - 1
                    )  # the first base to pull from the reference (SeqIO is 0 based index)
                    end = begin + 20  # the last base to pull from the reference
                    left = begin - 20  # the leftmost base to pull from the 5` side
                    flanking3 = str(ref["verified_pk184"][begin:end])
                    flanking5 = str(ref["verified_pk184"][left:begin])
                    flanking = str(flanking5 + flanking3)

                elif insertion >= 1 and insertion <= 20:
                    l += 1
                    begin = (
                        insertion - 1
                    )  # the first base to pull from the reference (SeqIO is 0 based index)
                    end = begin + 20  # the last base to pull from the reference
                    tail = (20 - begin) * -1
                    flanking3 = str(ref["verified_pk184"][0:end])
                    flanking5 = str(ref["verified_pk184"][tail:])
                    flanking = str(flanking5 + flanking3)

                elif insertion >= 2406:
                    h += 1
                    # special case: the insertion is exactly between the end and the begining of pk184
                    if insertion == 2425:
                        flanking5 = str(ref["verified_pk184"][2405:2425])
                        flanking3 = str(ref["verified_pk184"][0:20])
                        flanking = str(flanking5 + flanking3)

                    # When the insertion near the very end of the plasmid
                    elif insertion < 2425:
                        begin = (
                            insertion - 1
                        )  # the first base to pull from the reference (SeqIO is 0 based index)
                        left = begin - 20
                        flanking3 = str(ref["verified_pk184"][begin:])
                        flanking_looped = str(
                            ref["verified_pk184"][: 20 - len(flanking3)]
                        )
                        flanking5 = str(ref["verified_pk184"][left:begin])
                        flanking = str(flanking5 + flanking3 + flanking_looped)

                    # when the insertion point is greater than the length of the plasmid
                    else:
                        begin = (
                            insertion - 2426
                        )  # the first base to pull from the reference (SeqIO is 0 based index)
                        end = begin + 20  # the last base to pull from the reference
                        tail = (20 - begin) * -1
                        flanking3 = str(ref["verified_pk184"][0:end])
                        flanking5 = str(ref["verified_pk184"][tail:])
                        flanking = str(flanking5 + flanking3)

                else:
                    sys.exit("error handling insertion point in reverse oriention")

                # check that flanking has been asigned and it is the correct length
                if flanking == None:
                    sys.exit("flanking region not assigned")

                if len(flanking) != 40:
                    print(
                        "REVERSE",
                        flanking,
                        str(insertion),
                        flanking5,
                        flanking3,
                        flanking_looped,
                        str(begin),
                        str(left),
                    )
                    print(toks)
                    sys.exit("flanking not correct length")

                print(
                    flanking,
                    str(insertion),
                    flanking5,
                    flanking3,
                    toks[5],
                    toks[3],
                    sep="\t",
                    file=reverse_out,
                )
                continue

            else:
                f += 1
                fail_log.append(getframeinfo(currentframe()).lineno)
                print(toks[5], file=fail_cigars)
                continue

print(str(m) + " reads with POS > 20", file=log)
print(str(l) + " reads with POS < 20", file=log)
print(str(h) + " reads with POS < 2404", file=log)
print(str(no_ref) + " reads with * in RNAME", file=log)
print(str(u) + " reads with flag == 4", file=log)
print(str(sa) + " reads with flag == 2048", file=log)
print(str(sarc) + " reads with flag == 2064", file=log)
print(str(complex_cigar) + " reads with complex cigar strings", file=log)
print(str(f) + " failed CIGAR strings", file=log)
print("CIGARs failed on lines:", file=log)
print(Counter(fail_log).keys(), file=log)
print("with n occurances", file=log)
print(Counter(fail_log).values(), file=log)
print("SUCCESS")
