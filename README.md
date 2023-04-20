# INS-lessons
INS-lessons

# The code for INS
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
## trimmomatic ##
java -jar trimmomatic-0.39.jar PE -threads 25 030506923_S6_L001_R1_001.fastq.gz 030506923_S6_L001_R2_001.fastq.gz 030506923_f_paired.fq.gz 030506923_f_unpaired.fq.gz 030506923_r_paired.fq.gz 030506923_r_unpaired.fq.gz SLIDINGWINDOW:4:20 MINLEN:100 ;
java -jar trimmomatic-0.39.jar PE -threads 25 030510223_S8_L001_R1_001.fastq.gz 030510223_S8_L001_R2_001.fastq.gz 030510223_f_paired.fq.gz 030510223_f_unpaired.fq.gz 030510223_r_paired.fq.gz 030510223_r_unpaired.fq.gz SLIDINGWINDOW:4:20 MINLEN:100 ;
java -jar trimmomatic-0.39.jar PE -threads 25 030510823_S4_L001_R1_001.fastq.gz 030510823_S4_L001_R2_001.fastq.gz 030510823_f_paired.fq.gz 030510823_f_unpaired.fq.gz 030510823_r_paired.fq.gz 030510823_r_unpaired.fq.gz SLIDINGWINDOW:4:20 MINLEN:100 ;
mkdir trimm ; 
mv *_paired.fq.gz trimm/ ; 
rm *_unpaired.fq.gz ;
cd trimm/ ;
fastqc *.gz -t 12 ; 
mkdir fastqc ;
mv *.html *.zip fastqc/ ; 
ls -lh ;
# ¿qué calidad tienen ahora los reads? #
# ¿cuales son sus longitudes? #

## fastq-mcf (alternative) ##
fastq-mcf NexteraPE-PE.fa 030510223_f_paired.fq.gz 030510223_r_paired.fq.gz -o 030510223_f_clean.fq.gz -o 030510223_r_clean.fq.gz -q 20 -x 20 ;
fastq-mcf NexteraPE-PE.fa 030506923_f_paired.fq.gz 030506923_r_paired.fq.gz -o 030506923_f_clean.fq.gz -o 030506923_r_clean.fq.gz -q 20 -x 20 ;
fastq-mcf NexteraPE-PE.fa 030510823_f_paired.fq.gz 030510823_r_paired.fq.gz -o 030510823_f_clean.fq.gz -o 030510823_r_clean.fq.gz -q 20 -x 20 ;
fastqc *clean.fq.gz -t 12 ; 
mkdir fastqc ;
mv *.html *.zip fastqc/ ; 
ls -lh ; 

# 3 # spades assembly ##
spades -1 030506923_f_paired.fq.gz -2 030506923_r_paired.fq.gz --careful -o 030506923_spades -t 25 ;
spades -1 030510223_f_paired.fq.gz -2 030510223_r_paired.fq.gz --careful -o 030510223_spades -t 25 ;
spades -1 030510823_f_paired.fq.gz -2 030510823_r_paired.fq.gz --careful -o 030510823_spades -t 25 ;
ls -lh ; 

# 4 # cambiar nombres ##
mv 030506923_spades/contigs.fasta 030506923_spades/030506923_contigs.fasta  ; 
mv 030506923_spades/scaffolds.fasta 030506923_spades/030506923_scaffolds.fasta ; 
mv 030506923_spades/030506923_contigs.fasta 030506923_spades/030506923_scaffolds.fasta . ; 

mv 030510823_spades/contigs.fasta 030510823_spades/030510823_contigs.fasta  ; 
mv 030510823_spades/scaffolds.fasta 030510823_spades/030510823_scaffolds.fasta ; 
mv 030510823_spades/030510823_contigs.fasta 030510823_spades/030510823_scaffolds.fasta . ;

mv 030510223_spades/contigs.fasta 030510223_spades/030510223_contigs.fasta  ; 
mv 030510223_spades/scaffolds.fasta 030510223_spades/030510223_scaffolds.fasta ; 
mv 030510223_spades/030510223_contigs.fasta 030510223_spades/030510223_scaffolds.fasta . ;
ls -lh ; 

# ¿cuántos son los contigs resultantes? # 
# ¿qué longitudes presentan? # 

# 5 # preparar la data ##
grep ">" sequences.fasta | wc -l ;
grep ">" sequences.fasta ; 
sed 's/\ |.*//g' sequences.fasta > sequences2.fasta ; 
grep ">" sequences2.fasta ; 

# 6 # añadir IDs a cada contigs de cada muestra ##
# install.packages("seqinr") #
library(seqinr) ;
r <- dir() ;
head <- gsub("_contigs.fasta","",r) ;
a <- 0 ;
for (i in 1:length(head)){ ;
  a <- read.fasta(r[i]) ;
  names(a) <- paste0(rep(head[i],length(a)),"_",names(a)) ;
  write.fasta(a,names(a), file.out=paste0(head[i],".fas")) ;
} ;
q("no") ;

cat *.fas > contigs.fasta

# 7 # be patient ... codes in construction ##

```

# The code for INS
```r
# 8 # identificacion taxonomica de los reads # 
for r1 in *fastq.gz
do
prefix=$(basename $r1 _R1.fastq.gz)
r2=${prefix}_R2.fastq.gz
kraken2 --paired --use-names --gzip-compressed --db /home/veronica/Programas/database/kraken2/k2_viral_20230314 --threads 25 $r1 $r2 --report ${prefix}_report.txt --output ${prefix}_kraken2.out ;
done ;
mkdir kraken_run2 ;
mv *.out *.txt kraken_run2/ ;
cd kraken_run2/ ; 
ls -lh ; 
```

# Don't touch this, don't even look 
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
