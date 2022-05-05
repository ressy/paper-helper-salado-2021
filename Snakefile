from csv import DictReader, DictWriter
from pathlib import Path
from math import ceil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

TABLES = {}
for path in Path("from-paper").glob("*.csv"):
    with open(path) as f_in:
        TABLES[path.stem] = list(DictReader(f_in))

EXTRA = {}
for path in Path("metadata").glob("*.csv"):
    with open(path) as f_in:
        EXTRA[path.stem] = list(DictReader(f_in))

SAMPLE_MAP = {row["Sample Name"]: row["Run"] for row in EXTRA["SraRunTable"]}

### CHIIMP

rule chiimp_run:
    output: "analysis/chiimp/results/results.rds"
    input:
        samples=expand("data/{sample}.{read}.fastq.gz", sample=SAMPLE_MAP.keys(), read=["1", "2"]),
        config="analysis/chiimp/config.yml",
        samples_csv="analysis/chiimp/samples.csv",
        locus_attrs="analysis/chiimp/locus_attrs.csv"
    threads: 12 # separately defined in config.yml
    shell:
        "cd $(dirname {input.config}) && chiimp $(basename {input.config})"

rule chiimp_config:
    output: "analysis/chiimp/config.yml"
    input: "chiimp_config.yml"
    shell: "cp {input} {output}"

rule chiimp_make_samples:
    output: "analysis/chiimp/samples.csv"
    input:
        sra="metadata/SraRunTable.csv",
        loci="analysis/chiimp/locus_attrs.csv"
    run:
        with open(input.loci) as f_in:
            loci = list(DictReader(f_in))
        fields = ["Sample", "Replicate", "Locus", "Filename", "Sex", "Subspecies"]
        rows_out = []
        with open(input[0]) as f_in:
            for row in DictReader(f_in):
                for locus_row in loci:
                    sample, rep = row["Sample Name"].split("_")
                    file_path = Path(".").resolve() / "analysis/trim" / (row["Sample Name"] + ".1.fastq.gz")
                    rows_out.append({
                        "Sample": sample,
                        "Replicate": rep,
                        "Locus": locus_row["Locus"],
                        "Filename": file_path})
        rows_out = sorted(rows_out, key=lambda row: (row["Sample"], row["Replicate"], row["Locus"]))
        with open(output[0], "wt") as f_out:
            writer = DictWriter(
                f_out,
                fieldnames=["Sample", "Replicate", "Locus", "Filename"],
                lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows_out)

rule chiimp_make_locus_attrs:
    output: "analysis/chiimp/locus_attrs.csv"
    input:
        tableS1="from-paper/tableS1.csv",
        extra="metadata/loci.csv"
    params:
        length_buffer_fraction=0.2
    run:
        fields = [
            "Locus", "LengthMin", "LengthMax", "LengthBuffer",
            "Motif", "Primer", "ReversePrimer"]
        def getsize(row):
            txt = row["Allele Size (bp)"]
            parts = [int(x) if x != "" else 0 for x in re.sub("[^-0-9]*", "", txt).split("-")]
            if len(parts) < 2:
                parts = [parts[0] - 10, parts[0] + 10]
            # I *think* the "allele size" they give is without primers.  I
            # think.
            extra = len(row["Forward primer (5' -> 3')"]) + len(row["Reverse primer (5' -> 3')"])
            parts = [p + extra for p in parts]
            return parts
        extra = {}
        with open(input.extra) as f_in:
            for row in DictReader(f_in):
                extra[row["Locus"]] = row
        with open(input.tableS1) as f_in, open(output[0], "wt") as f_out:
            writer = DictWriter(f_out, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            for row in DictReader(f_in):
                if row["Selected"] == "X":
                    parts = getsize(row)
                    # 20% of average of length range, rounded up to nearest 10
                    length_buffer = (parts[1] + parts[0])/2*params.length_buffer_fraction
                    length_buffer = 10*ceil(length_buffer/10)
                    writer.writerow({
                        "Locus": row["Locus name"],
                        "LengthMin": parts[0],
                        "LengthMax": parts[1],
                        "LengthBuffer": length_buffer,
                        "Motif": extra[row["Locus name"]]["Motif"],
                        "Primer": row["Forward primer (5' -> 3')"],
                        "ReversePrimer": row["Reverse primer (5' -> 3')"]})

### MEGASAT

rule megasat:
    output: directory("analysis/megasat/output")
    input:
        data=expand("analysis/megasat/input/{sample}.{read}.fastq", sample=SAMPLE_MAP.keys(), read=["1", "2"]),
        primers="analysis/megasat/primers.tsv"
    params:
        input_dir="analysis/megasat/input",
        min_mismatches=2,
        min_depth=5,
        output_dir=""
    threads: 4
    # From the perl source if six args are given:
    # ($inputPrimers,$mismatches,$m_depth,$max_processors,$dataset,$saveDir)= @ARGV;
    shell:
        """
            mkdir {output}
            perl "MEGASAT/MEGASAT_1.0 for Linux/MEGASAT_Genotype.pl" {input.primers} \
                {params.min_mismatches} {params.min_depth} {threads} \
                {params.input_dir} {output}
        """

rule megasat_input:
    output: "analysis/megasat/input/{sample}.{read}.fastq"
    input: "data/{sample}.{read}.fastq.gz"
    shell: "gunzip < {input} > {output}"

rule megasat_make_primers_file:
    output: "analysis/megasat/primers.tsv"
    input:
        tableS1="from-paper/tableS1.csv",
        extra="metadata/loci.csv"
    run:
        fields = [
            "Locus Name",
            "5' Microsatellite primer",
            "reverse-complement of 3' microsatellite primer",
            "3' flank", "5' flank",
            "repeat_unit sequence",
            "Ratios group"]
        extra = {}
        with open(input.extra) as f_in:
            for row in DictReader(f_in):
                extra[row["Locus"]] = row
        with open(input.tableS1) as f_in, open(output[0], "wt") as f_out:
            writer = DictWriter(f_out, fieldnames=fields, lineterminator="\n", delimiter="\t")
            writer.writeheader()
            for row in DictReader(f_in):
                if row["Selected"] == "X":
                    writer.writerow({
                        "Locus Name": row["Locus name"],
                        "5' Microsatellite primer": row["Forward primer (5' -> 3')"],
                        "reverse-complement of 3' microsatellite primer": row["Reverse primer (5' -> 3')"],
                        "3' flank": extra[row["Locus name"]]["FlankFwd"],
                        "5' flank": extra[row["Locus name"]]["FlankRev"],
                        "repeat_unit sequence": extra[row["Locus name"]]["Motif"],
                        "Ratios group": ""})

### Trim adapters

rule cutadapt_all:
    input: expand("analysis/trim/{sample}.{read}.fastq.gz", sample=SAMPLE_MAP.keys(), read=["1", "2"])

# Near the end of R1 I see:
# GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
# And just before R1:
# ACACTCTTTCCCTACACGACGCTCTTCCGATCT
# Thse match TruSeq adapters here:
# http://tucf-genomics.tufts.edu/documents/protocols/TUCF_Understanding_Illumina_TruSeq_Adapters.pdf
# https://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2013/06/illumina-adapter-sequences_1000000002694-00.pdf
rule cutadapt:
    output:
        report="analysis/trim/{sample}.report.json",
        reads=expand("analysis/trim/{{sample}}.{read}.fastq.gz", read=[1, 2])
    input: expand("data/{{sample}}.{read}.fastq.gz", read=[1, 2])
    params:
        adapter_fwd="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        adapter_rev="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" # or revcmp??
    threads: 4
    shell:
        """
            cutadapt --json {output.report} \
                --cores {threads} \
                -a {params.adapter_fwd} \
                -A {params.adapter_rev} \
                -o {output.reads[0]} -p {output.reads[1]} \
                {input}
        """

### SRA download

rule prep_sample_fastq_all:
    input: expand("data/{sample}.{read}.fastq.gz", sample=SAMPLE_MAP.keys(), read=["1", "2"])

rule prep_sample_fastq:
    output: "data/{sample}.{read}.fastq.gz"
    input: lambda w: expand("sra/{SRR}_{{read}}.fastq.gz", SRR=SAMPLE_MAP[w.sample])
    shell: "ln -s ../{input} {output}"

rule get_sra_all:
    input: expand("sra/{SRR}_{read}.fastq.gz", SRR=SAMPLE_MAP.values(), read=["1", "2"])

rule sra_fqgz:
    output: "sra/{SRR}_{read}.fastq.gz"
    input: "sra/{SRR}_{read}.fastq"
    shell: "gzip < {input} > {output}"

rule get_sra:
    output: temp(expand("sra/{{SRR}}_{read}.fastq", read=[1, 2]))
    params: outdir="sra"
    shell: "fastq-dump --split-3 -O {params.outdir} {wildcards.SRR}"

### Investigating reads by primer

rule primer_split_all:
    input: expand("analysis/primer-split/{thing}/{thing}-{primer}.1.fastq.gz", thing=SAMPLE_MAP.keys(), primer=[row["Locus name"] for row in TABLES["tableS1"]])

rule primer_split:
    output: expand("analysis/primer-split/{{thing}}/{{thing}}-{primer}.1.fastq.gz", primer=[row["Locus name"] for row in TABLES["tableS1"]])
    input:
        reads="analysis/trim/{thing}.1.fastq.gz",
        primers="analysis/primer-split/primers_fwd.fasta"
    threads: 4
    params:
        min_overlap=15,
        # action for matching adapter seq:
        # none means don't actually change the sequences
        action="none"
    shell:
        """
            cutadapt \
                --cores {threads} \
                --overlap {params.min_overlap} \
                --action {params.action} \
                -g file:{input.primers} \
                -o analysis/primer-split/{wildcards.thing}/{wildcards.thing}-{{name}}.1.fastq.gz \
                {input.reads}
        """

rule primer_split_prep_fasta:
    output: "analysis/primer-split/primers_fwd.fasta"
    input: "from-paper/tableS1.csv"
    run:
        with open(input[0]) as f_in, open(output[0], "wt") as f_out:
            for row in DictReader(f_in):
                SeqIO.write(
                    SeqRecord(
                        Seq("^" + row["Forward primer (5' -> 3')"]),
                        id=row["Locus name"],
                        description=""),
                    f_out,
                    "fasta-2line")
