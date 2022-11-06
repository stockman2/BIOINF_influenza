# тестируем bash-скрипты)))
# URL = 'http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/SRR1705851.fastq.gz'
configfile: 'config.yaml'


rule variant_download:
	output:
		'{variant}.fastq.gz'
	params: 
		URL = config["URL_variant"]
	shell:
		'wget -O {output} {params.URL}'

rule variant_unzip:
	input:
		'{variant}.fastq.gz'
	output:
		'{variant}.fastq'
	shell:
		'gzip -d {input} > {output}'

rule bwa_ref_index:
    input:
        "{reference}.fasta"
    output:
        "{reference}.fasta.amb",
        "{reference}.fasta.ann",
        "{reference}.fasta.bwt",
        "{reference}.fasta.pac",
        "{reference}.fasta.sa",
    shell:
        "bwa index {input}"

rule bwa_align:
    input:
        "reference.fasta.amb",
        "reference.fasta.ann",
        "reference.fasta.bwt",
        "reference.fasta.pac",
        "reference.fasta.sa",
        ref = "reference.fasta",
        variant = "{variant}.fastq"
    output:
        "{variant}.bam"
    shell:
        "bwa mem {input.ref} {input.variant} | samtools view -S -b - | samtools sort -o {output}"
        
# rule flagstat:
#	input:
#		'{variant}.bam'
#	output:
#		'{variant}.txt'
#	shell:
#		'samtools flagstat {input} > {output}'
		
rule samtools_index:
	input:
		'{variant}.bam'
	output:
		'{variant}.bam.bai'
	shell:
		'samtools index {input} > {output}'
		
rule mpileup:
	input:
		'{variant}.bam.bai',
		bam = '{variant}.bam',
		ref = 'reference.fasta'
	output:
		'{variant}.mpileup'
	shell:
		'samtools mpileup -f {input.ref} {input.bam} -d 0 > {output}'
		
rule varscan:
	input:
		'{variant}.mpileup'
	output:
		'{variant}.vcf'
	params:
		rate = config["mutation_rate"]
	shell:
		'java -jar "/home/freeman/Workplace/Prac/VarScan.jar" mpileup2snp {input} --min-var-freq {params.rate} --variants --output-vcf 1 > {output}'
		
# rule vcf_parse:
#	input:
#		'{variant}.vcf'
#	output:
#		'{variant}_vcf_parsed.txt'
#	shell:
#		"cat {input} | awk 'NR>24 {print $2, $4, $5, {split($10,a,":"); print a[7]}}' > {output}"

