SAMPLES = ['HSM67VF9', 'HSM67VFD', 'HSM67VFJ', 'HSM6XRQB',
           'HSM6XRQI', 'HSM6XRQK', 'HSM6XRQM', 'HSM6XRQO',
           'HSM7CYY7', 'HSM7CYY9', 'HSM7CYYB', 'HSM7CYYD']
RADIUS = ['1', '5', '10']

rule all:
    input: 
        expand("outputs/groot/{sample}_arg90_report.txt", sample = SAMPLES),
        expand("outputs/groot/{sample}_proportion_arg90_derived.csv", sample = SAMPLES),
        "outputs/arg90_matches/arg90_matches_dedup.fna",
        expand("outputs/bcalm/{sample}_r{radius}/cfxA4_AY769933.fna.cdbg_ids.reads.gz.unitigs.gfa", sample = SAMPLES, radius = RADIUS)

rule fastp:
    input:
        R1="inputs/raw/{sample}_R1.fastq.gz",
        R2="inputs/raw/{sample}_R2.fastq.gz",
    output:
        R1="outputs/fastp/{sample}_R1.trim.fq.gz",
        R2="outputs/fastp/{sample}_R2.trim.fq.gz",
        json = "outputs/fastp/{sample}.json",
        html = "outputs/fastp/{sample}.html",
    conda: "envs/env.yml"
    threads: 1
    resources:
        mem_mb=16000
    shell:'''
    fastp --in1 {input.R1} \
      --in2 {input.R2} \
      --out1 {output.R1} \
      --out2 {output.R2} \
      --detect_adapter_for_pe \
      --qualified_quality_phred 4 \
      --length_required 31 --correction \
      --json {output.json} \
      --html {output.html}
    '''

rule remove_host:
# http://seqanswers.com/forums/archive/index.php/t-42552.html
# https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?usp=sharing
    output:
        r1 = 'outputs/bbduk/{sample}_R1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{sample}_R2.nohost.fq.gz',
        human_r1='outputs/bbduk/{sample}_R1.human.fq.gz',
        human_r2='outputs/bbduk/{sample}_R2.human.fq.gz'
    input: 
        r1="outputs/fastp/{sample}_R1.trim.fq.gz",
        r2="outputs/fastp/{sample}_R2.trim.fq.gz",
        human='inputs/host/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz'
    conda: 'envs/env.yml'
    threads: 1
    resources:
        mem_mb=64000
    shell:'''
    bbduk.sh -Xmx64g t=3 in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.human_r1} outm2={output.human_r2} k=31 ref={input.human}
    '''

rule kmer_trim_reads:
    input: 
        'outputs/bbduk/{sample}_R1.nohost.fq.gz',
        'outputs/bbduk/{sample}_R2.nohost.fq.gz'
    output: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    conda: 'envs/env.yml'
    threads: 1
    resources:
        mem_mb=61000
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 60e9 -V - -o {output}
    '''

rule download_groot_db:
    output: directory("arg-annot.90")
    conda: "envs/groot.yml"
    threads: 1
    resources:
        mem_mb=1000
    shell:'''
    groot get -d arg-annot
    '''

rule index_groot_db:
    output: "groot-index"
    input: directory("arg-annot.90")
    conda: "envs/groot.yml"
    threads: 8
    resources:
        mem_mb=64000
    shell:'''
    groot index -m {input} -i {output} -w 100 -p 8
    '''

rule run_groot:
    input: 
        index="groot-index",
        abundtrim="outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output:
        graph=directory("outputs/groot/{sample}_arg90.graph"),
        bam="outputs/groot/{sample}_arg90.bam"
    conda: "envs/groot.yml"
    threads: 2
    resources:
        mem_mb=16000
    shell:'''
    groot align -i {input.index} -f {input.abundtrim} -p 2 -g {output.graph} > {output.bam}
    '''

rule report_groot:
    input: "outputs/groot/{sample}_arg90.bam"
    output: "outputs/groot/{sample}_arg90_report.txt"
    conda: "envs/groot.yml"
    threads: 1
    resources:
        mem_mb=16000
    shell:'''
    groot report --bamFile {input} -c .9 > {output}
    '''

rule combine_report_groot:
    input: expand("outputs/groot/{sample}_arg90_report.txt", sample = SAMPLES)
    output: "outputs/groot/all_arg90_report.txt"
    threads: 1
    resources:
        mem_mb=1000
    shell:'''
    cat {input} > {output}
    '''

rule calc_proportion_arg_reads_per_sample:
    input: bam ="outputs/groot/{sample}_arg90.bam"
    output: "outputs/groot/{sample}_proportion_arg90_derived.csv"
    conda: "envs/samtools.yml"
    threads: 1
    resources:
        mem_mb=2000
    shell:'''
    mapped=$(samtools view -c -F 4 {input.bam})
    proportion=$(echo "scale=10 ; $mapped / 10000000" | bc)
    printf "{wildcards.sample},$proportion" > {output}
    '''
    
rule download_arg_db:
    output: "inputs/arg_db/argannot-args.fna"
    threads: 1
    resources:
        mem_mb=1000
    shell:'''
    wget -O {output} https://github.com/will-rowe/groot/raw/master/db/full-ARG-databases/arg-annot-db/argannot-args.fna
    '''

rule faidx_arg_db:
    input: "inputs/arg_db/argannot-args.fna"
    output: "inputs/arg_db/argannot-args.fna.fai"
    conda: "envs/samtools.yml"
    threads: 1
    resources:
        mem_mb=1000
    shell:'''
    samtools faidx {input}
    '''

rule extract_groot_matches:
    input:
        groot_report="outputs/groot/all_arg90_report.txt",
        arg90="inputs/arg_db/argannot-args.fna",
        arg90_faidx="inputs/arg_db/argannot-args.fna"
    output: "outputs/arg90_matches/arg90_matches.fna"
    threads: 1
    resources:
        mem_mb=4000
    conda: "envs/samtools.yml"
    shell:'''
    samtools faidx {input.arg90} `cut -f1 {input.groot_report}` > {output}
    '''

rule remove_duplicate_groot_matches:
    input: "outputs/arg90_matches/arg90_matches.fna"
    output: "outputs/arg90_matches/arg90_matches_dedup.fna"
    conda: "envs/seqkit.yml"
    threads: 1
    resources:
        mem_mb=4000
    shell:'''
    cat {input} | seqkit rmdup -s -o {output}
    '''

rule spacegraphcats_one_arg:
    input: 
        query = "outputs/arg90_matches/cfxA4_AY769933.fna", 
        conf = "outputs/sgc_conf/{sample}_r{radius}_conf.yml",
        reads = "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output:
        "outputs/sgc_arg_queries_r{radius}/{sample}_k31_r{radius}_search_oh0/cfxA4_AY769933.fna.cdbg_ids.reads.gz",
        "outputs/sgc_arg_queries_r{radius}/{sample}_k31_r{radius}_search_oh0/cfxA4_AY769933.fna.contigs.sig"
    params: outdir = lambda wildcards: "outputs/sgc_arg_queries_r" + wildcards.radius
    conda: "envs/spacegraphcats.yml"
    resources:
        mem_mb = 64000
    threads: 1
    shell:'''
    python -m spacegraphcats run {input.conf} extract_contigs extract_reads --nolock --outdir={params.outdir} --rerun-incomplete 
    '''

rule bcalm_nbhd:
    input: "outputs/sgc_arg_queries_r{radius}/{sample}_k31_r{radius}_search_oh0/cfxA4_AY769933.fna.cdbg_ids.reads.gz"
    output: "outputs/bcalm/{sample}_r{radius}/cfxA4_AY769933.fna.cdbg_ids.reads.gz.unitigs.fa"
    params: "outputs/bcalm/{sample}_r{radius}/cfxA4_AY769933.fna.cdbg_ids.reads.gz"
    resources:
        mem_mb = 4000
    threads: 1
    conda: "envs/spacegraphcats.yml"
    shell:'''
    bcalm -in {input} -out-dir outputs/bcalm/{wildcards.sample}_{wildcards.radius} -kmer-size 31 -abundance-min 1 -out {params}
    '''

rule convert_to_gfa:
    input: "outputs/bcalm/{sample}_r{radius}/cfxA4_AY769933.fna.cdbg_ids.reads.gz.unitigs.fa"
    output: "outputs/bcalm/{sample}_r{radius}/cfxA4_AY769933.fna.cdbg_ids.reads.gz.unitigs.gfa"
    resources:
        mem_mb = 2000
    threads: 1
    conda: "envs/spacegraphcats.yml"
    shell:'''
    python ./scripts/convertToGFA.py {input} {output} 31
    '''
