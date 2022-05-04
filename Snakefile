from csv import DictReader

with open("metadata/SraRunTable.csv") as f_in:
    RUNINFO = list(DictReader(f_in))
SAMPLE_MAP = {row["Sample Name"]: row["Run"] for row in RUNINFO}

rule prep_sample_fastq_all:
    input: expand("data/{sample}.{read}.fastq.gz", sample=SAMPLE_MAP.keys(), read=["1", "2"])

rule prep_sample_fastq:
    output: "data/{sample}.{read}.fastq.gz"
    input: lambda w: expand("sra/{SRR}_{{read}}.fastq.gz", SRR=SAMPLE_MAP[w.sample])
    shell: "ln -s ../{input} {output}"

rule get_sra_all:
    input: expand("sra/{SRR}_{read}.fastq.gz", SRR=[row["Run"] for row in RUNINFO], read=["1", "2"])

rule sra_fqgz:
    output: "sra/{SRR}_{read}.fastq.gz"
    input: "sra/{SRR}_{read}.fastq"
    shell: "gzip < {input} > {output}"

rule get_sra:
    output: temp("sra/{SRR}_{read}.fastq")
    shell: "fastq-dump --split-3 -O $(dirname {output}) {wildcards.SRR}"
