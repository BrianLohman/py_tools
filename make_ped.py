import sys

# USAGE
# python make_ped.py [list-of-crams.txt]

# ASSUMPTIONS
# 1. cram list is sdout of ls operation: SSXXXXX.final.cram, one sample per line

# build dictionaries
bam_to_simons = {}
for line in open("sample_id_dictionary.csv"):
    toks = line.strip().split(",")
    simonsID = toks[0]
    bamID = toks[1]
    bam_to_simons[bamID] = simonsID

simons_to_bam = {}
for line in open("sample_id_dictionary.csv"):
    toks = line.strip().split(",")
    simonsID = toks[0]
    bamID = toks[1]
    simons_to_bam[simonsID] = bamID

# specifiy input/output files
crams = sys.argv[1]
out = open(crams.split("_")[0] + ".ped", "w")

# use the structure of the sample ids to build ped
for line in open(crams):
    toks = line.strip().split(".")
    bamID = toks[0]

    # assign simons sample id
    simonsID = bam_to_simons[bamID]

    # assign family and individual ID
    family = str(simonsID[:-3])
    individualID = bamID

    # assign mother/father
    role = simonsID[6:]

    if role == "mo" or role == "fa":
        mother = "0"
        father = "0"
    else:
        father = simons_to_bam[family + ".fa"]
        mother = simons_to_bam[family + ".mo"]

    # assign phenotype
    if role == "p1":
        pheno = "1"
    else:
        pheno = "0"

    # assign path to cram
    path = (
        "s3://sscwgs-hg38/" + crams.split("_")[0] + "/" + individualID + ".final.cram"
    )

    # assign path to crai
    crai = (
        "s3://sscwgs-hg38/crais/"
        + crams.split("_")[0]
        + "/"
        + individualID
        + ".final.cram.crai"
    )

    # write to file
    out.write(
        "\t".join(
            [family, individualID, father, mother, "0", str(pheno), path, crai, "\n"]
        )
    )
