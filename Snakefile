SAMPLES = ['HSM67VF9', 'HSM67VFD', 'HSM67VFJ', 'HSM6XRQB',
           'HSM6XRQI', 'HSM6XRQK', 'HSM6XRQM', 'HSM6XRQO',
           'HSM7CYY7', 'HSM7CYY9', 'HSM7CYYB', 'HSM7CYYD']

rule all:
    input: expand("outputs/groot/{sample}_arg90.bam", sample = SAMPLES)

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
    output: "arg-annot.90"
    conda: "envs/groot.yml"
    threads: 1
    resources:
        mem_mb=1000
    shell:'''
    groot get -d arg-annot
    '''

rule index_groot_db:
    output: "groot-index"
    input: "arg-annot.90"
    conda: "envs/groot.yml"
    threads: 8
    resources:
        mem_mb=64000
    shell:'''
    groot index -m arg-annot.90 -i groot-index -w 100 -p 8
    '''

rule run_groot:
    input: 
        index="groot-index",
        abundtrim="outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output:
        graph="outputs/groot/{sample}_arg90.graph",
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
    output:
    conda: "envs/groot.yml"
    threads: 1
    resources:
        mem_mb=16000
    shell:'''
    groot report --bamFile {input} -c .9
    '''
