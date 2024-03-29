#!/usr/bin/env python
# Build and submit jobs on CHPC
# STATUS: Working. Updated 18 Feb 2021
# USAGE: python batch_chpc.py -n [job name] -c [commands] -w [ntasks/node] -m [mem/task]
import argparse
import os
import re
import sys

parser = argparse.ArgumentParser(
    description="generate and submit jobs to CHPC based on list of commands"
)
parser.add_argument(
    "-l", "--ml", dest="ml", nargs="+", help="modules to load for each slurm script"
)
parser.add_argument("-n", "--job-name", dest="job", help="job base name")
parser.add_argument(
    "-j", "--commands", dest="commands", help="list of commands to run, one per line"
)
parser.add_argument(
    "-w",
    "--ntasks-node",
    dest="ntasks",
    default=8,
    help="number of tasks per node, default is 8 tasks",
)
parser.add_argument(
    "-m",
    "--mem",
    dest="mem",
    default=32000,
    help="minimum memory required per node, default = 32GB",
)
parser.add_argument(
    "-p",
    "--partition",
    dest="partition",
    default="hci-rw",
    help="hci-rw is only available choice",
)
parser.add_argument(
    "-t",
    "--time",
    dest="time",
    default="12:00:00",
    type=str,
    help="max run time, default is 12 hours",
)
parser.add_argument(
    "-d",
    "--dry-run",
    dest="submit",
    default=True,
    action="store_false",
    help="sbatch job files. or dont.",
)
parser.add_argument(
    "-c",
    "--cpus_task",
    dest="cpus_task",
    default=1,
    help="number of cpus required for each task",
)
parser.add_argument(
    "-x",
    "--wait",
    dest="wait",
    default=False,
    action="store_true",
    help="add wait as last line for background jobs",
)
args = parser.parse_args()

# count the number of lines in the commands file
f = open(args.commands)
contents = f.readlines()
n_commands = 0

for line in contents:
    stripped_line = line.strip()
    if re.match("^\s*#", stripped_line):
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
        o = open(job_fh + ".slurm", "w")
        if (args.ml) is not None:
            o.write(
                "\n".join(
                    [
                        "#!/bin/bash",
                        "#SBATCH --time=" + args.time,
                        "#SBATCH --account=hci-rw",
                        "#SBATCH --partition=" + args.partition,
                        "#SBATCH --nodes=1",
                        "#SBATCH --ntasks=" + str(args.ntasks),
                        "#SBATCH --cpus-per-task=" + str(args.cpus_task),
                        "#SBATCH --mem=" + str(args.mem),
                        "#SBATCH --job-name=" + job_fh,
                        "#SBATCH -o call-" + job_fh + ".out",
                        "#SBATCH -e call-" + job_fh + ".err",
                        "\n",
                        "module load " + " ".join(args.ml),
                        str(line),
                    ]
                )
            )
        else:
            o.write(
                "\n".join(
                    [
                        "#!/bin/bash",
                        "#SBATCH --time=" + args.time,
                        "#SBATCH --account=hci-rw",
                        "#SBATCH --partition=" + args.partition,
                        "#SBATCH --nodes=1",
                        "#SBATCH --ntasks=" + str(args.ntasks),
                        "#SBATCH --cpus-per-task=" + str(args.cpus_task),
                        "#SBATCH --mem=" + str(args.mem),
                        "#SBATCH --job-name=" + job_fh,
                        "#SBATCH -o call-" + job_fh + ".out",
                        "#SBATCH -e call-" + job_fh + ".err",
                        "\n",
                        str(line),
                    ]
                )
            )

        # if there is only one command to run
        if line_count == n_commands:
            o.close()
            if args.submit == True:
                os.system("sbatch " + job_fh + ".slurm")
            break

        else:
            line_count += 1
            continue

    else:
        if line_count == n_commands:
            o.write(str(line))
            if args.wait == True:
                o.write("wait")

            o.close()

            if args.submit == True:
                os.system("sbatch " + job_fh + ".slurm")
            break

        if line_count % int(args.ntasks) != 0:
            o.write(str(line))
            line_count += 1
            continue

    if line_count % int(args.ntasks) == 0:
        o.write(str(line))
        if args.wait == True:
            o.write("wait")

        o.close()

        if args.submit == True:
            os.system("sbatch " + job_fh + ".slurm")

        job_fh = None
        line_count += 1
        continue
