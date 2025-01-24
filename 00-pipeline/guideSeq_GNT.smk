import re
import os
from glob import glob
import pandas as pd
from snakemake.utils import validate


if not os.path.isfile("guideSeq_GNT.yml"):
    raise SystemExit("\n  No config file found in current directory \n")

configfile: "guideSeq_GNT.yml"
#validate(config, "config.schema.yaml")
config_file_path=os.getcwd()+'/'+workflow.configfiles[0]


## Load the sample information files as a TSV file. 'sampleName' column is mandatory
samples = pd.read_table(config["sampleInfo_path"],sep=";").set_index("sampleName", drop=False)
validate(samples, "samples.schema.yaml")


my_list=samples[["Genome","gRNA_sequence","PAM_sequence"]].drop_duplicates()



genomes=my_list["Genome"].tolist()
gRNAs=my_list["gRNA_sequence"].tolist()
PAMs=my_list["PAM_sequence"].tolist()
genomes_unique=list(set(genomes))

#print(samples)
## RUN the pipeline in the project folder.
## snakemake -s ../00-pipeline/guideSeq_GNT.smk -k -j 12 --use-conda --conda-front-end mamba --conda-prefix ../01-envs -n

##########################################################
##########################################################

rule target:
    input: "05-Report/report.html",
     ["03-align/{sample}.UMI.ODN.trimmed.filtered.multi.realigned.bed".format(sample=sample) for sample in samples["sampleName"]]
     #expand("06-offPredict/{gen}_{grna}{pam}.csv", zip, grna=gRNAs, gen=genomes,pam=PAMs),
     #["03-align/{sample}_multi.txt".format(sample=sample) for sample in samples["sampleName"]]
     
     

# Merge index1 and 2 in a new fastq file I3. Easier for demultiplexing.
rule merge_indexes:
    input: I1=config["I1"] , I2=config["I2"]
    output: I3=temp("I3.fastq.gz")
    conda: "../01-envs/env_tools.yml"
    shell: """
        python ../00-pipeline/merge_I1_I2.py {input.I1} {input.I2} {output}
        """

# make a barcode file with I1 & I2 sequences for demultiplexing.
rule make_indexes_fasta:
    input: 
    output: temp("demultiplexing_barcodes.fa")
    run: 
        with open("demultiplexing_barcodes.fa", 'w') as fasta_file:
            for index, row in samples.iterrows():
                sequence = row['index1'] + row['index2']
                fasta_file.write('>' + row['sampleName'] + '\n')
                fasta_file.write(sequence + '\n')


# demultiplexe libraries based on barcodes in the sampleInfo file and I3 file generated before
## 
# check whether an empty library will make the pip to crash (ie wrong index)
##
rule demultiplex_library:
    input: R1=config["R1"], R2=config["R2"], I1=config["I1"], I2=config["I2"],I3=rules.merge_indexes.output.I3, barcodes=rules.make_indexes_fasta.output
    output: 
     R1=["00-demultiplexing/{sample}_R1.fastq.gz".format(sample=sample) for sample in samples["sampleName"]],
     R2=["00-demultiplexing/{sample}_R2.fastq.gz".format(sample=sample) for sample in samples["sampleName"]],
     I3=["00-demultiplexing/{sample}_I3.fastq.gz".format(sample=sample) for sample in samples["sampleName"]]
    conda: "../01-envs/env_tools.yml"
    log: R1= "00-demultiplexing/demultiplexing_R1.log", R2= "00-demultiplexing/demultiplexing_R2.log"
    threads: 6
    shell: """
        
        cutadapt -g ^file:{input.barcodes} -j {threads} -e 0  --action trim --no-indels -o 00-demultiplexing/{{name}}_I3.fastq.gz -p 00-demultiplexing/{{name}}_R1.fastq.gz {input.I3} {input.R1} > {log.R1}
        cutadapt -g ^file:{input.barcodes} -j {threads} -e 0  --action trim --no-indels -o 00-demultiplexing/{{name}}_I3.fastq.gz -p 00-demultiplexing/{{name}}_R2.fastq.gz {input.I3} {input.R2} > {log.R2}
        cutadapt -g ^file:{input.barcodes} -j {threads} -e 0  --action trim --no-indels -o 00-demultiplexing/{{name}}_I3.fastq.gz -p 00-demultiplexing/{{name}}_I1.fastq.gz {input.I3} {input.I1} > {log.R1}
        cutadapt -g ^file:{input.barcodes} -j {threads} -e 0  --action trim --no-indels -o 00-demultiplexing/{{name}}_I3.fastq.gz -p 00-demultiplexing/{{name}}_I2.fastq.gz {input.I3} {input.I2} > {log.R1}
        rm 00-demultiplexing/*unknown* 
        """
        
rule add_UMI:
    input: R1="00-demultiplexing/{sample}_R1.fastq.gz", R2="00-demultiplexing/{sample}_R2.fastq.gz",I3="00-demultiplexing/{sample}_I3.fastq.gz"
    output: R1=temp("00-demultiplexing/{sample}_R1.UMI.fastq.gz"),
     R2=temp("00-demultiplexing/{sample}_R2.UMI.fastq.gz") ,
     I3=temp("00-demultiplexing/{sample}_I3.UMI.fastq.gz")
    threads:6
    conda: "../01-envs/env_tools.yml"
    params: UMI=config["UMI_pattern"] ## bp in 3' of index to considere as UMI
    shell: """
        
        UMI_length=$(expr length {params.UMI})
        
        cutadapt -j {threads} -u -$UMI_length --rename='{{id}}_{{r1.cut_suffix}} {{comment}}' -o {output.I3} -p {output.R1} {input.I3} {input.R1}
        cutadapt -j {threads} -u -$UMI_length --rename='{{id}}_{{r1.cut_suffix}} {{comment}}' -o {output.I3} -p {output.R2} {input.I3} {input.R2}
        """

# remove ODN and discard reads without ODN
rule trim_ODN:
    input: R1=rules.add_UMI.output.R1, R2=rules.add_UMI.output.R2
    output: R1=temp("01-trimming/{sample}_R1.UMI.ODN.fastq.gz"),
     R2=temp("01-trimming/{sample}_R2.UMI.ODN.fastq.gz")
    threads: 6
    log:"01-trimming/{sample}.odn.log"
    conda: "../01-envs/env_tools.yml"
    message: "removing ODN sequence, discard reads without ODN sequence"
    params: R2_leading=lambda wildcards: config[samples["type"][wildcards.sample]][samples["orientation"][wildcards.sample]]["R2_leading"]
    shell: """
        cutadapt -j {threads} -G "{params.R2_leading};max_error_rate=0;rightmost" --discard-untrimmed  -o {output.R1} -p {output.R2} {input.R1} {input.R2} > {log}
        """


# remove leading and trailing ODN sequences
rule trim_reads:
    input: R1=rules.trim_ODN.output.R1, R2=rules.trim_ODN.output.R2
    output: R1=temp("01-trimming/{sample}_R1.UMI.ODN.trimmed.fastq.gz"), 
     R2=temp("01-trimming/{sample}_R2.UMI.ODN.trimmed.fastq.gz")
    threads: 6
    log: "01-trimming/{sample}.trailing.log"
    conda: "../01-envs/env_tools.yml"
    message: "trimming ODN and adaptor sequences in reads"
    params: R2_leading=lambda wildcards: config[samples["type"][wildcards.sample]][samples["orientation"][wildcards.sample]]["R2_leading"], 
     R2_trailing=lambda wildcards: config[samples["type"][wildcards.sample]][samples["orientation"][wildcards.sample]]["R2_trailing"],
     R1_trailing=lambda wildcards: config[samples["type"][wildcards.sample]][samples["orientation"][wildcards.sample]]["R1_trailing"]
    shell: """
        cutadapt -j {threads} -A "{params.R2_trailing};min_overlap=6;max_error_rate=0.1" -a "{params.R1_trailing};min_overlap=6;max_error_rate=0.1" -o {output.R1} -p {output.R2} {input.R1} {input.R2} > {log}
        """



# make a size selection before mapping
rule filter_reads:
    input: R1=rules.trim_reads.output.R1, R2=rules.trim_reads.output.R2
    output:  R1="02-filtering/{sample}_R1.UMI.ODN.trimmed.filtered.fastq.gz", R2="02-filtering/{sample}_R2.UMI.ODN.trimmed.filtered.fastq.gz",
     R1short=temp("02-filtering/{sample}_R1.UMI.ODN.trimmed.tooshort.fastq.gz"), R2short=temp("02-filtering/{sample}_R2.UMI.ODN.trimmed.tooshort.fastq.gz")
    threads: 6
    log: "02-filtering/{sample}.filter.log"
    params: length=config["minLength"]
    conda: "../01-envs/env_tools.yml"
    message: "remove pairs if one mate is shorter than x bp"
    shell: """
        cutadapt -j {threads} --pair-filter=any --minimum-length {params.length} --too-short-output {output.R1short} --too-short-paired-output  {output.R2short} --output {output.R1} --paired-output {output.R2} {input.R1} {input.R2} > {log}
        """




rule make_genome_index:
    input: lambda wildcards: config["genome"][wildcards.genome]["fasta"]
    output: "../02-ressources/{genome}.log"
    threads: 6
    conda: "../01-envs/env_tools.yml"
    params: idx=lambda wildcards: config["genome"][wildcards.genome]["index"]
    shell: """
        bowtie2-build --threads {threads} {input} {params.idx} > {log}
        """




# map reads on the reference genome as pairs
if (config["aligner"]  == "bowtie2" or config["aligner"]  == "Bowtie2") :
    rule alignOnGenome:
        input: R1=rules.filter_reads.output.R1, R2=rules.filter_reads.output.R2
        output: sam=temp("03-align/{sample}.UMI.ODN.trimmed.filtered.sam")
        threads: 6
        log: "03-align/{sample}.UMI.ODN.trimmed.filtered.align.log"
        conda: "../01-envs/env_tools.yml"
        message: "Aligning PE reads on genome with Bowtie2"
        params: index=lambda wildcards: config["genome"][samples["Genome"][wildcards.sample]]["index"],
         minFragLength=config["minFragLength"],
         maxFragLength=config["maxFragLength"]
        shell: """
            bowtie2 -p {threads} --no-unal -I {params.minFragLength} -X {params.maxFragLength} --dovetail --no-mixed --no-discordant --un-conc-gz 03-align/{wildcards.sample}_R%.UMI.ODN.trimmed.unmapped.fastq.gz   -x {params.index} -1 {input.R1} -2 {input.R2} -S {output.sam} 2> {log}
        """
        

elif (config["aligner"]  == "bwa" or config["aligner"]  == "Bwa") :
    rule alignOnGenome:
        input: R1=rules.filter_reads.output.R1, R2=rules.filter_reads.output.R2
        output: sam=temp("03-align/{sample}.UMI.ODN.trimmed.filtered.sam"), unfilterdsam=temp("03-align/{sample}.UMI.ODN.trimmed.unfiltered.sam")
        threads: 6
        log: "03-align/{sample}.UMI.ODN.trimmed.filtered.align.log"
        conda: "../01-envs/env_tools.yml"
        message: "Aligning PE reads on genome with BWA mem"
        params:  index=config["genome"][lambda wildcards:samples["Genome"][wildcards.sample]]["index"]
        shell: """
            #bwa mem -t {threads} {params.index} {input.R1} {input.R2} > {output.sam} 2> {log}
            bwa mem -t {threads} {params.index} {input.R1} {input.R2} > {output.unfilterdsam} 2> {log}
            samtools view -F 0x4 -F 0x8 -F 0x100 -F 0x800 -f 0x2 -b {output.unfilterdsam} > {output.sam}
        """


rule get_R2_multihits_reads:
    input: sam=rules.alignOnGenome.output.sam
    output: list="03-align/{sample}_multi.txt"
    threads: 2
    conda: "../01-envs/env_tools.yml"
    message: "identify R2 reads that have multiple hits with equal alignment score"
    shell: """
        samtools view -f 131 -@ {threads} {input} | awk '
            function abs(v) {{return v < 0 ? -v : v}}
            {{
                as=""; xs=""
                for(i=12; i<=13; i++) {{
                    if($i ~ /^AS:i:/) as=substr($i,6)
                    if($i ~ /^XS:i:/) xs=substr($i,6)
                }}
                if(xs!="") {{
                    if(abs(xs)== abs(as)) {{
                         print $1
                        }}
                    }}
            }}' > {output.list}
    """

rule split_bam:
    input: rnames=rules.get_R2_multihits_reads.output.list, bam = rules.alignOnGenome.output.sam
    output: uniq="03-align/{sample}.UMI.ODN.trimmed.filtered.unique.bam", multi="03-align/{sample}.UMI.ODN.trimmed.filtered.multi.bam"
    threads: 6
    conda: "../01-envs/env_tools.yml"
    shell: """
        samtools view -hb -@ {threads} -N {input.rnames} -o {output.multi} -U {output.uniq} {input.bam}
    """

rule convert_bamtofastq:
    input: rules.split_bam.output.multi
    output: R1="03-align/{sample}.UMI.ODN.trimmed.filtered.multi.R1.fastq.gz", R2="03-align/{sample}.UMI.ODN.trimmed.filtered.multi.R2.fastq.gz"
    threads: 2
    conda: "../01-envs/env_tools.yml"
    shell: """
        samtools fastq -@ {threads} -1 {output.R1} -2 {output.R2}  {input}
    """




rule realign_multi:
    input: R1=rules.convert_bamtofastq.output.R1, R2=rules.convert_bamtofastq.output.R2
    output: sam=temp("03-align/{sample}.UMI.ODN.trimmed.filtered.multi.realigned.sam")
    threads: 6
    log: "03-align/{sample}.UMI.ODN.trimmed.filtered.realign.log"
    conda: "../01-envs/env_tools.yml"
    message: "Aligning PE reads on genome with Bowtie2"
    params: index=lambda wildcards: config["genome"][samples["Genome"][wildcards.sample]]["index"],
     minFragLength=config["minFragLength"],
     maxFragLength=config["maxFragLength"],
     k=config["multi_alignments"]
    shell: """
            bowtie2 -p {threads} -k {params.k} --no-unal -I {params.minFragLength} -X {params.maxFragLength} --dovetail --no-mixed --no-discordant --un-conc-gz 03-align/{wildcards.sample}_R%.UMI.ODN.trimmed.unmapped.fastq.gz   -x {params.index} -1 {input.R1} -2 {input.R2} -S {output.sam} 2> {log}
        """


rule multi_processing:
    input: rules.realign_multi.output.sam
    output: bed="03-align/{sample}.UMI.ODN.trimmed.filtered.multi.realigned.bed", info="03-align/{sample}.UMI.ODN.trimmed.filtered.multi.realigned.info"
    conda: "../01-envs/env_tools.yml"
    shell: """
    
        samtools view -h -f 131 {input} | bedtools bamtobed -i - -ed > {output.bed}
        samtools view -f 131 {input} | awk 'BEGIN{{OFS="\\t"}} {{print $1,$3,$4,$9,$12}}' > {output.info}
    """




#### finish here




# sort alignments by names (required for BEDPE conversion) and position (for viewing)
rule sort_aligned:
    input: sam=rules.alignOnGenome.output.sam
    output: bamPos="03-align/{sample}.UMI.ODN.trimmed.filtered.bam",  bamName="03-align/{sample}.UMI.ODN.trimmed.filtered.sortedName.bam"
    threads: 6
    log:
    conda: "../01-envs/env_tools.yml"
    message: "Sort reads by name"
    shell: """
        samtools sort -@ {threads} {input.sam}  > {output.bamPos}
        samtools index {output.bamPos}
        samtools sort -n -@ {threads} {input.sam} > {output.bamName}
    """


# convert BAM file to BEDPE and identify insertion point, read length
# filter alignment with MAPQ score > threshold in config file
# aggregate reads per fragment size and genomic coordinates
# add a cluster ID to group close IS (distance defined in the config file) 
rule call_IS:
    input: rules.sort_aligned.output.bamName
    output: tmp=temp("04-IScalling/{sample}.pebed"), 
     umi="04-IScalling/{sample}.reads_per_UMI_per_IS.bed"
    threads: 1
    conda: "../01-envs/env_tools.yml"
    log:
    params:  minMAPQ=config["minMAPQ"], UMI=config["UMI_pattern"]
    shell: """
    
        UMI_length=$(expr length {params.UMI})

        bedtools bamtobed -bedpe -mate1 -i {input} > {output.tmp}
        
        # read R2 contains the genome/ODN junction (use $10 of PE bed = mate 2)
        
        ##################################################################
        # count number of reads per UMI and per IS (pos and strand separated)
        ##################################################################

        awk -v umi="$UMI_length" 'BEGIN{{OFS="\\t";FS="\\t"}} ($8>={params.minMAPQ}) && ($1 == $4) {{if($10=="+") print $4,$5,$5,$7,substr($7,length($7)-umi+1,umi),$8,$10,$3-$5; else print $4,$6-1,$6-1,$7,substr($7,length($7)-umi+1,umi),$8,$10,$6-$2}}' {output.tmp} |  sort -k1,1 -k2,3n -k5,5 -k7,7  | bedtools groupby -g 1,2,3,5,7 -c 4,6 -o count_distinct,median > {output.umi}
        
        """

rule correct_UMI:
    input: bed=rules.call_IS.output.umi
    output: "04-IScalling/{sample}.reads_per_UMI_per_IS_corrected.bed"
    conda: "../01-envs/env_R4.3.2.yml"
    threads: 2
    params: hamming=config["UMI_hamming_distance"], filter = config["UMI_filter"],  UMI=config["UMI_pattern"], method = config["UMI_deduplication"]
    shell: """
        Rscript ../00-pipeline/correct_umi.R {input} {params.UMI} {params.filter} {params.hamming} {params.method} {output}
        """


rule Collapse_UMI_IS:
    input: rules.correct_UMI.output
    output: collapse="04-IScalling/{sample}.UMIs_per_IS_in_Cluster.bed"
    threads: 1
    conda: "../01-envs/env_tools.yml"
    log:
    params: window=config["ISbinWindow"], minMAPQ=config["minMAPQ"],minReadsPerUMI=config["minReadsPerUMI"],minUMIPerIS=config["minUMIPerIS"]
    shell: """
        
        ##################################################################
        # aggregate UMI per IS
        ##################################################################
        
        awk 'NR>1 && $6>{params.minReadsPerUMI}' {input} | sort -k1,1 -k2,2n -k3,3n -k5,5 | bedtools groupby -g 1,2,3,5 -c 4,6,4,6,7  -o count_distinct,sum,collapse,collapse,median  | awk 'BEGIN{{OFS="\\t";FS="\\t"}} $5>{params.minUMIPerIS} {{print $1,$2,$3,"ID_"NR,$9,$4,$5,$6,$7,$8}}' |  sort -k1,1 -k2,2n -k3,3n -k5,5 | bedtools cluster -d {params.window} > {output.collapse}
        """
    



rule get_chrom_length:
    input: lambda wildcards: config["genome"][wildcards.genome]["fasta"]
    output: "../02-ressources/{genome}.chrom.length.log"
    threads: 1
    conda: "../01-envs/env_tools.yml"
    params: lambda wildcards: config["genome"][wildcards.genome]["fasta"]+".fai"
    shell: """
        samtools faidx {input}
        """

        
rule get_fasta_around_is:
    input: bed=rules.Collapse_UMI_IS.output.collapse, chr_length=lambda wildcards: config["genome"][samples["Genome"][wildcards.sample]]["fasta"]+".fai"
    output: cluster="04-IScalling/{sample}.cluster_slop.bed", fa="04-IScalling/{sample}.cluster_slop.fa"
    threads: 1
    conda: "../01-envs/env_tools.yml"
    log:
    params: fasta=lambda wildcards: config["genome"][samples["Genome"][wildcards.sample]]["fasta"], 
     slop_size=config["slopSize"]
    shell: """
        #echo -e "#chromosome\tstart\tend\tclusterID\tN_orientations\tmedianMAPQ\tN_IS\tN_UMI\tN_reads" > {output.cluster}
        sort -k11n {input.bed} | bedtools groupby -g 11 -c 1,2,3,4,5,6,7,8 -o distinct,min,max,count_distinct,median,count_distinct,sum,sum |  awk 'BEGIN{{OFS="\\t"}}{{print $2,$3,$4,$1,$7,$6,$5,$8,$9}}' |  bedtools slop -i - -b {params.slop_size} -g {input.chr_length} > {output.cluster} 
        bedtools getfasta -name -fi {params.fasta} -bed {output.cluster} > {output.fa}
        """ 


rule get_stats_fq:
    input: rules.demultiplex_library.output.R1, rules.add_UMI.output.R1, rules.trim_ODN.output.R1, rules.trim_reads.output.R1, rules.filter_reads.output.R1
    output: "05-Report/{sample}.stat"
    conda: '../01-envs/env_tools.yml'
    threads: 6
    shell: """
        seqkit stat -j {threads} -a -T {input} > {output}
        """


### in future, simplify this step to avoid predicting several times the same gRNA
# group by sequence and PAM


rule predict_offtarget_SWOffinder:
    input: 
    output: txt="06-offPredict/{gen}_{grna}{pam}.txt", csv="06-offPredict/{gen}_{grna}{pam}.csv"
    conda: '../01-envs/env_tools.yml'
    params: genome=lambda wildcards: config["genome"][wildcards.gen]["fasta"],
        gRNA=lambda wildcards:{wildcards.grna},
        maxE=config["max_edits_crRNA"],               #Max edits allowed (integer).
        maxM=config["max_edits_crRNA"],               #Max mismatches allowed without bulges (integer).
        maxMB=config["max_edits_crRNA"],             #Max mismatches allowed with bulges (integer).
        maxB=config["SWoffFinder"]["maxB"],               #Max bulges allowed (integer).
        window_size=config["SWoffFinder"]["window_size"], #The window size for choosing the best in a window
        bulges=config["tolerate_bulges"],
        PAM=lambda wildcards:wildcards.pam,
        SWoffFinder=config["SWoffFinder"]["path"]
    threads: 12
    shell: """
        echo {params.gRNA}{params.PAM} > {output.txt}
        
        if [ {params.bulges} = TRUE ] ; 
        then 
            bulge_size={params.maxB}; 
        else 
            bulge_size=0;
        fi

        
        java -cp {params.SWoffFinder}/bin SmithWatermanOffTarget.SmithWatermanOffTargetSearchAlign {params.genome} {output.txt} 06-offPredict/{wildcards.gen} {params.maxE} {params.maxM} {params.maxMB} $(echo $bulge_size) {threads} TRUE {params.window_size} {params.PAM} TRUE
        """


rule report_data:
    input: fasta=rules.get_fasta_around_is.output.fa, cluster=rules.get_fasta_around_is.output.cluster, bed=rules.Collapse_UMI_IS.output.collapse, statfq=rules.get_stats_fq.output,statal=rules.alignOnGenome.log
    output: "05-Report/{sample}.rdata"
    conda: "../01-envs/env_R4.3.2.yml"
    log:
    threads: 4
    params: gRNA_seq=lambda wildcards:samples["gRNA_sequence"][wildcards.sample], 
     gRNA_name=lambda wildcards:samples["gRNA_name"][wildcards.sample],
     PAM=lambda wildcards:samples["PAM_sequence"][wildcards.sample],
     offset=lambda wildcards:samples["Cut_Offset"][wildcards.sample],
     max_edits=config["max_edits_crRNA"],
     bulges=config["tolerate_bulges"]
    shell : """
        Rscript ../00-pipeline/multiple_alignments.R {input.fasta} {input.cluster} {input.bed} {params.gRNA_seq} {params.gRNA_name}  {params.PAM}  {params.offset} {params.max_edits} {params.bulges} {output}
        """
        

rule prepare_annotations:
    input: 
    output: ["../02-ressources/{species}.rds".format(species=species) for species in genomes_unique]
    threads: 6
    conda: "../01-envs/env_R4.3.2.yml"
    params: config=os.path.abspath(workflow.configfiles[0])
    shell: """
        Rscript ../00-pipeline/prepare_annotations.R {params.config}
    """
    

rule annotate_sites:
    input: rdata=rules.report_data.output, species=rules.prepare_annotations.output
    output: "05-Report/{sample}_summary.xlsx"
    conda: "../01-envs/env_R4.3.2.yml"
    threads: 4
    params: species=lambda wildcards: samples["Genome"][wildcards.sample]
    shell: """
        Rscript ../00-pipeline/annotate_cuting_sites.R {params.species} {input.rdata} {output}
        """


rule report:
    input: summaries= ["05-Report/{sample}_summary.xlsx".format(sample=sample) for sample in samples["sampleName"]],
     predictions=expand("06-offPredict/{gen}_{grna}{pam}.csv", zip, grna=gRNAs, gen=genomes,pam=PAMs)
    output: report_rdata= "05-Report/report.rdata"
    conda: "../01-envs/env_R4.3.2.yml"
    threads: 1
    params: sampleInfo=config["sampleInfo_path"], 
     config=os.path.abspath(workflow.configfiles[0]), 
     minUMI_alignments_figure=config["minUMI_alignments_figure"],
     min_predicted_distance=config["min_predicted_distance"],
     max_clusters=config["max_clusters"]
    shell: """
        mkdir -p 05-Report/report-files/;
        Rscript ../00-pipeline/generate_tables.r "{input.summaries}" {params.sampleInfo} {params.config} "{input.predictions}" {params.max_clusters} {params.minUMI_alignments_figure} {params.min_predicted_distance}
        """

rule make_report:
    input:  rules.report.output.report_rdata
    output: report_html="05-Report/report.html"
    conda: "../01-envs/env_R4.3.2.yml"
    params: config_path=config_file_path
    threads: 1
    shell: """
        Rscript ../00-pipeline/publish_report.r {input} {params.config_path} {output}
    """


onsuccess:
    print("Workflow finished, no error")
    shell("cowsay -e \♥♥  You rock !! ")

onerror:
    print("An error occurred")
    shell(" cowsay -e \X\X  Houston, we have a problem !")
