# Build job launchers for quinlan-kp partion of CHPC
import sys
import os

# USAGE: python batch_chpc.py [job name] [commands] [ntasks/node] [mem/task]

# read in job name
job = sys.argv[1]

# define commands, the list of tasks to execute
commands = sys.argv[2]

# number of tasks per node
ntasks = sys.argv[3]

# memory per task
mem = sys.argv[4]

# loop though commands, printing toks at the bottom of the standard header
job_fh = None
line_count = int(1)

for line in open(commands):
    if job_fh == None:
        job_fh = job + "_" + str(line_count)
        o = open(job_fh + ".job", "w")
        o.write('\n'.join(["#!/bin/bash", "#SBATCH --time=12:00:00", "#SBATCH --account=quinlan-kp", "#SBATCH --partition=quinlan-shared-kp", "#SBATCH --nodes=1", "#SBATCH --ntasks="+str(ntasks), "#SBATCH --mem="+str(mem), "#SBATCH --job-name="+job_fh, "#SBATCH -o call-"+job_fh+".out", "#SBATCH -e call-"+job_fh+".err", '\n', str(line)]))
        line_count += 1
        continue

    else:
	if line_count % int(ntasks) != 0:
	    o.write(str(line))
            line_count += 1
            continue

	if line_count % int(ntasks) == 0: 
	    o.write(str(line))
	    o.close()
            os.system("sbatch "+job_fh+".job")
	    job_fh = None
	    line_count += 1
	    continue
