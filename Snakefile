from csv import DictReader

with open("metadata/SraRunTable.csv") as f_in:
    RUNINFO = list(DictReader(f_in))

rule get_sra_all:
    input: expand("data/{SRR}_{read}.fastq.gz", SRR=[row["Run"] for row in RUNINFO], read=["1", "2"])

rule sra_fqgz:
    output: "data/{SRR}_{read}.fastq.gz"
    input: "data/{SRR}_{read}.fastq"
    shell: "gzip < {input} > {output}"

rule get_sra:
    output: temp("data/{SRR}_{read}.fastq")
    shell: "fastq-dump --split-3 -O $(dirname {output}) {wildcards.SRR}"
