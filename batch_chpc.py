#!/usr/bin/env python

# Build and submit jobs on CHPC
# STATUS: working: 29 Nov 18
# USAGE: python batch_chpc.py -n [job name] -c [commands] -w [ntasks/node] -m [mem/task]

import sys
import os
import argparse
import re

parser = argparse.ArgumentParser(description='generate and submit jobs to CHPC based on list of commands')
parser.add_argument('-n','--job-name', dest='job', help = 'job base name')
parser.add_argument('-c', '--commands', dest='commands', help = 'list of commands to run, one per line')
parser.add_argument('-w', '--ntasks-node', dest='ntasks', default = 8, help = 'number of tasks per node')
parser.add_argument('-m', '--mem', dest='mem', default = 128000, help = 'minimum memory required per node')
parser.add_argument('-p', '--partition', dest = 'partition', default = "quinlan-shared-kp", help = 'quinlan-kp/quinlan-shared-kp')
parser.add_argument('-t', '--time', dest='time', default = "12:00:00", type = str, help = 'max run time')
args = parser.parse_args()

# count the number of lines in the commands file
f = open(args.commands)
contents = f.readlines()
n_commands = 0
for line in contents:
    stripped_line = line.strip()
    if re.match('^\s*#', stripped_line):
        continue
    if not stripped_line:
	continue
    else:
	n_commands += 1

# loop though commands, printing lines from the commands file at the bottom of the standard header
# while accounting for tasks per node constraints
job_fh = None
line_count = int(1)

for line in open(args.commands):
    if job_fh == None:
        job_fh = args.job + "_" + str(line_count)
        o = open(job_fh + ".job", "w")
        o.write('\n'.join(["#!/bin/bash", "#SBATCH --time="+args.time, "#SBATCH --account=quinlan-kp", "#SBATCH --partition="+args.partition, "#SBATCH --nodes=1", "#SBATCH --ntasks="+str(args.ntasks), "#SBATCH --mem="+str(args.mem), "#SBATCH --job-name="+job_fh, "#SBATCH -o call-"+job_fh+".out", "#SBATCH -e call-"+job_fh+".err", '\n', str(line)]))
	# if there is only one command to run
	if line_count == n_commands:
            o.close()
	    os.system("sbatch "+job_fh+".job")
	    break

	else:
            line_count += 1
	    continue
    
    else:
        if line_count == n_commands:
	    o.write(str(line))
	    o.close()
            os.system("sbatch "+job_fh+".job")
	    break
	
        if line_count % int(args.ntasks) != 0:
	    o.write(str(line))
            line_count += 1
            continue

	if line_count % int(args.ntasks) == 0: 
	    o.write(str(line))
	    o.close()
            os.system("sbatch "+job_fh+".job")
	    job_fh = None
	    line_count += 1
	    continue

