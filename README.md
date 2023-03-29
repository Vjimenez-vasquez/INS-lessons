# INS-lessons
INS-lessons

# The code 1
```r
# 1 # calidad ##
## fastqc ##
fastqc *.gz -t 12 ; 
mkdir fastqc ;
mv *.html *.zip fastqc/ ; 
ls -lh ;

## ¿phred-score? ##
## ¿qué indican los plots de FASTQC? ##
## ¿qué calidad promedio tienen los reads generados? ##
## ¿qué longitudes promedio tienen los reads generados? ##

# 2 # trimming ##
# http://www.usadellab.org/cms/index.php?page=trimmomatic # 
java -jar trimmomatic-0.39.jar PE 030506923_S6_L001_R1_001.fastq.gz 030506923_S6_L001_R2_001.fastq.gz 030506923_f_paired.fq.gz 030506923_f_unpaired.fq.gz 030506923_r_paired.fq.gz 030506923_r_unpaired.fq.gz LEADING:5 TRAILING:5 MINLEN:140 SLIDINGWINDOW:4:20 ;
java -jar trimmomatic-0.39.jar PE 030510223_S8_L001_R1_001.fastq.gz 030510223_S8_L001_R2_001.fastq.gz 030510223_f_paired.fq.gz 030510223_f_unpaired.fq.gz 030510223_r_paired.fq.gz 030510223_r_unpaired.fq.gz LEADING:5 TRAILING:5 MINLEN:140 SLIDINGWINDOW:4:20 ;
java -jar trimmomatic-0.39.jar PE 030510823_S4_L001_R1_001.fastq.gz 030510823_S4_L001_R2_001.fastq.gz 030510823_f_paired.fq.gz 030510823_f_unpaired.fq.gz 030510823_r_paired.fq.gz 030510823_r_unpaired.fq.gz LEADING:5 TRAILING:5 MINLEN:140 SLIDINGWINDOW:4:20 ;
mkdir trimm ; 
mv *paired.fq.gz trimm/ ; 
rm *unpaired.fq.gz ; 
ls -lh ; 

```

# The code 
```r
## cleaning ##
java -jar trimmomatic.jar PE forward.fastq.gz reverse.fastq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz LEADING:3 TRAILING:3 MINLEN:200 ;
mv output_forward_paired.fq.gz forward_paired.fq.gz ;  
mv output_reverse_paired.fq.gz reverse_paired.fq.gz ;
fastqc forward_paired.fq.gz reverse_paired.fq.gz -t 4 ;
fastq-mcf NexteraPE-PE.fa forward_paired.fq.gz reverse_paired.fq.gz -o uno.gz -o dos.gz -q 20 -x 20 ;
ls -lh ;


## assembly ##
spades -1 uno.fq.gz -2 dos.fq.gz --careful -o bb -t 4 ;
grep ">" bb/scaffolds.fasta ;
ls -lh ;

## order contigs and obtain the pseudo-genome  in contiguator ##

## mapping ##
bwa index PseudoContig.fasta ;

#2# preparar las instrucciones generales#
for r1 in *fq.gz
do
prefix=$(basename $r1 _forward.fq.gz)
r2=${prefix}_reverse.fq.gz

#3# instrucciones para generar el archivo .bam#
bwa mem -t 4 PseudoContig.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 4 -bS -T PseudoContig.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 4 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 4 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 4 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup -@ 4 ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index -@ 4 ${prefix}.bam ;
rm ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;
done ;
ls -lh ;

## mapping quality ##
./qualimap bamqc -bam 120020820.bam -outfile 120020820.pdf

## genome annotation ##
conda activate prokka_env ;
prokka PseudoContig.fasta --outdir 120020820.annotation --cpus 4 --genus Bartonella ;
ls -lh ;

## pan genome analysis ##
#extraer las regiones comunes: PANGENOME#
roary -p 4 -f roary -cd 95 -e -n -v *gff ;
cp roary/core_gene_alignment.aln ;
snp-sites -m -o snps.fasta -c core_gene_alignment.aln ; 
aliview snps.fasta ;
ls -lh ;
```
