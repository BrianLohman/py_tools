# Build job launchers for quinlan-kp partion of CHPC
import sys
import os

# USAGE: python build_launcher [job name] [commands]

# read in job name
job = sys.argv[1]

# define commands, the list of tasks to execute
commands = sys.argv[2]

# loop though commands, printing toks at the bottom of the standard header
for i, line in enumerate(open(commands)):
    job_fh = job + "_" + str(i)
    o = open(job_fh + ".job", "w")
    o.write('\n'.join(["#!/bin/bash", "#SBATCH --time=12:00:00", "#SBATCH --account=quinlan-kp", "#SBATCH --partition=quinlan-kp", "#SBATCH --nodes=1", "#SBATCH --ntasks=1", "#SBATCH --mem=128000", "#SBATCH --job-name="+job_fh, "#SBATCH -o call-"+job_fh+".out", "#SBATCH -e call-"+job_fh+".err", '\n', str(line)]))
    o.close()
    os.system("sbatch "+job_fh+".job")
