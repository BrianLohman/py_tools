import sys

# compare keys (sample IDs) between two files
# column that contains the key is specified by user

# USAGE: python id_check.py [file A] [file B] [file A key] [file B key ] 

# Assumptions:
# 1. files are tab delimited

file_a = sys.argv[1]
file_b = sys.argv[2]
key_a = int(sys.argv[3])
key_b = int(sys.argv[4])

# initialize number successful matches
n_success = 0

# open the files
file_a = open(sys.argv[1])
file_b = open(sys.argv[2])

# read in the first line
file_a_line = file_a.readline()
file_b_line = file_b.readline()

# initalize line counter
line_no = 0

# loop if either file has reached end of file
while file_a_line != '' or file_b_line != '':

    fields =  file_a_line.strip().split('\t')
    file_a_key = fields[key_a]

    fields =  file_b_line.strip().split('\t')
    file_b_key = fields[key_b]

    if file_a_key != file_b_key:
        print('\t'.join(["Fail", file_a_key, file_b_key]), file=sys.stderr)
    if file_a_key == file_b_key:
        n_success += 1

    # read in next line
    file_a_line = file_a.readline()
    file_b_line = file_b.readline()

    # increment line counter
    line_no += 1
    
print(' '.join([str(n_success), "matches out of", str(line_no)]))
