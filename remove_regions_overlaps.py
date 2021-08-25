# remove any overlaps among lines in a bed file used for regions extract
import sys

bed = sys.argv[1]

out = open("compact_" + sys.argv[1], "w")

for i, toks in enumerate(x.strip().split("\t") for x in open(bed)):
    # sanity check: ends must be > than starts
    if int(toks[1]) > int(toks[2]):
        print("start coordinate greater than end coordinate!")
        print(toks)
        break

    if i == 0:
        current_chrom = toks[0]
        current_start = toks[1]
        current_end = toks[2]
        continue

    # when two successive regions are on the same chromsome
    if toks[0] == current_chrom:
        # when the next region does NOT overlap with the current region
        if int(toks[1]) > int(current_end):
            out.write(
                "\t".join([current_chrom, str(current_start), str(current_end) + "\n"])
            )
            current_chrom = toks[0]
            current_start = toks[1]
            current_end = toks[2]

        # when the next region is totally contained in the current region
        if int(toks[1]) > int(current_start) and int(toks[2]) > int(current_end):
            continue

        # when the next region overlaps with the current region, extending the total region
        if int(toks[1]) < int(current_end) and int(toks[2]) > int(current_end):
            out.write("\t".join([current_chrom, str(current_end + 1), toks[2] + "\n"]))
            current_start = current_end + 1
            current_end = toks[2]

    # when two successive regions are NOT on the same chromosome
    if toks[0] != current_chrom:
        out.write(
            "\t".join([current_chrom, str(current_start), str(current_end) + "\n"])
        )
        current_chrom = toks[0]
        current_start = toks[1]
        current_end = toks[2]
