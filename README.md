# Sv4_HistoricalStarlings
Scripts and notes related to manuscript entitled "Historical museum samples enable the examination of divergent and parallel evolution during invasion"

Full curated scripts + additional notes will be uploaded shortly.

2021_11_20_notebook_10049: PopGen analysis

2021_11_20_notebook_10040: Outlier analysis

This page contains small vignette's on the following
<li>DArTseq SNP variant calling</li>
<li>Some basic popgen</li>

<h1> SNP variant calling for raw RRS (DArTseq) data</h1>

The following is a summary of some SNP calling pipelines that I have used to process raw DArTseq data, but they may be applied to other tyes of reduced representation sequencing (RRS) raw data.

Variant calling falls into a few primary steps:
<li>Mapping raw reads to a reference genome (or for refernce free variant calling forming the pseudogenome and calling references to this)</li>
<li>Calling the variant SNPs for each sample based on the mapped reads</li>
<li>Filtering the called SNPs</li>

For the first two steps of the process, you can mix and match mapping and variant calling programs.

<h3>Data cleaning: Stacks </h3>

<pre class="r"><code>
RAW_DIR=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/rawdata/batch2_split
OUTPUT_DIR=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/processing/samples_batch
BARCODE_DIR=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/processing/barcodes

for i in 1 2 3 4 5 6;
do
process_radtags -p ${RAW_DIR}/batch_$i/ -o ${OUTPUT_DIR}/batch_$i/ -b ${BARCODE_DIR}/bardcode_named_$i.txt -r -c -q --renz_1 pstI --renz_2 sphI
done
</code></pre>
This should produce files that look like:

Run QC:
<pre class="r"><code>
fastqc *.fq.gz --outdir=fastqc/
</code></pre>

Count reads:
<pre class="r"><code>
zcat *.fq.gz | echo $((`wc -l`/4))
</code></pre>

<h3>BWA (mapping)</h3>

Place your genome somewhere sensible and index it for BWA

<pre class="r"><code>
genome_fa=Sturnus_vulgaris_2.3.1.simp.fasta
bwa index -p bwa/stuv $genome_fa &> bwa/bwa_index.oe
</code></pre>

Aligning with bwa mem:
<pre class="r"><code>
for i in cat sample_names.txt;
do
sample=$i
echo WORKING WITH ${sample}
fq_file=../../samples_batch/merged/${sample}.fq.gz
bam_file=${sample}.bam
bwa_db=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/genome/bwa/stuv
bwa mem -t 16 -M $bwa_db $fq_file | samtools view -bS | samtools sort > $bam_file
done
</code></pre>

check reads mapped successfully: 
<pre class="r"><code>
samtools flagstat SAMPLENAME.bam
</code></pre>

pic of mapping data

Calling variants with gstacks:

<pre class="r"><code>
gstacks -I ./ -M ../historic_populations.txt -O ./
populations -P ./ -M ../historic_populations.txt --vcf
</code></pre>

Removed 0 loci that did not pass sample/population constraints from 412079 loci.
Kept 412079 loci, composed of 28225204 sites; 20 of those sites were filtered, 243589 variant sites remained.
    27345713 genomic sites, of which 846899 were covered by multiple loci (3.1%).
Mean genotyped sites per locus: 68.47bp (stderr 0.00).

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  aw: 11.272 samples per locus; pi: 0.17294; all/variant/polymorphic sites: 15515519/232371/142127; private alleles: 24046
  hist: 4.8985 samples per locus; pi: 0.12483; all/variant/polymorphic sites: 18008520/91625/34830; private alleles: 13284
  mv: 11.015 samples per locus; pi: 0.14832; all/variant/polymorphic sites: 19889591/226170/98305; private alleles: 4740
  mw: 11.146 samples per locus; pi: 0.16776; all/variant/polymorphic sites: 15527851/232281/128802; private alleles: 13279
  nc: 11.442 samples per locus; pi: 0.17365; all/variant/polymorphic sites: 15591084/235304/135530; private alleles: 16176
  or: 11.45 samples per locus; pi: 0.16313; all/variant/polymorphic sites: 14915003/232899/117912; private alleles: 7942

 

Calling SNPS with filtering:
Strict filtering:

populations -P ./ -M ../historic_populations.txt --vcf -r 0.50 --min_maf 0.01 --write_random_snp -t 8
Removed 271920 loci that did not pass sample/population constraints from 416460 loci.
Kept 144540 loci, composed of 10087650 sites; 67754 of those sites were filtered, 77412 variant sites remained.
    9994834 genomic sites, of which 91970 were covered by multiple loci (0.9%).
Mean genotyped sites per locus: 69.79bp (stderr 0.01).

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  aw: 13.444 samples per locus; pi: 0.15815; all/variant/polymorphic sites: 9110939/68361/47379; private alleles: 3605
  hist: 8.0858 samples per locus; pi: 0.15427; all/variant/polymorphic sites: 7780598/14725/8706; private alleles: 1541
  mv: 13.407 samples per locus; pi: 0.14296; all/variant/polymorphic sites: 7740125/63240/34156; private alleles: 1022
  mw: 13.325 samples per locus; pi: 0.15843; all/variant/polymorphic sites: 8904808/67473/44953; private alleles: 2671
  nc: 13.392 samples per locus; pi: 0.16369; all/variant/polymorphic sites: 8775160/71236/48457; private alleles: 3798
  or: 13.452 samples per locus; pi: 0.15581; all/variant/polymorphic sites: 8641133/69988/42823; private alleles: 2169 

lenient filtering:

populations -P ./ -M ../historic_populations.txt --vcf --min_maf 0.01 --write_random_snp -t 8
Removed 0 loci that did not pass sample/population constraints from 416460 loci.
Kept 416460 loci, composed of 28517922 sites; 42232 of those sites were filtered, 109260 variant sites remained.
    27621426 genomic sites, of which 862808 were covered by multiple loci (3.1%).
Mean genotyped sites per locus: 68.46bp (stderr 0.00).

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  aw: 10.178 samples per locus; pi: 0.20525; all/variant/polymorphic sites: 15518001/101327/65371; private alleles: 5356
  hist: 4.4469 samples per locus; pi: 0.16671; all/variant/polymorphic sites: 18101356/35776/14624; private alleles: 3257
  mv: 9.9247 samples per locus; pi: 0.17966; all/variant/polymorphic sites: 20244797/97302/48946; private alleles: 1800
  mw: 10.047 samples per locus; pi: 0.20144; all/variant/polymorphic sites: 15532433/101150/62414; private alleles: 3827
  nc: 10.379 samples per locus; pi: 0.20966; all/variant/polymorphic sites: 15593641/103102/66563; private alleles: 5388
  or: 10.392 samples per locus; pi: 0.19816; all/variant/polymorphic sites: 14918843/101541/59292; private alleles: 3053
Populations is done.

 

mid filtering:

populations -P ./ -M ../historic_populations.txt --vcf -r 0.30 --min_maf 0.01 --write_random_snp -t 8
Removed 242732 loci that did not pass sample/population constraints from 416460 loci.
Kept 173728 loci, composed of 12100028 sites; 60629 of those sites were filtered, 91091 variant sites remained.
    11947646 genomic sites, of which 150246 were covered by multiple loci (1.3%).
Mean genotyped sites per locus: 69.65bp (stderr 0.00).

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  aw: 12.391 samples per locus; pi: 0.16745; all/variant/polymorphic sites: 10485200/79413/53953; private alleles: 4418
  hist: 6.7673 samples per locus; pi: 0.16979; all/variant/polymorphic sites: 10354693/20446/11287; private alleles: 1948
  mv: 12.29 samples per locus; pi: 0.14624; all/variant/polymorphic sites: 9501521/74403/39154; private alleles: 1408
  mw: 12.223 samples per locus; pi: 0.1662; all/variant/polymorphic sites: 10370645/79162/51538; private alleles: 3273
  nc: 12.323 samples per locus; pi: 0.17553; all/variant/polymorphic sites: 10343134/83163/55866; private alleles: 4566
  or: 12.398 samples per locus; pi: 0.16368; all/variant/polymorphic sites: 10030684/81463/49117; private alleles: 2667

 

 

Aligning with bwa aln:
https://github.com/lh3/bwa#type

Illumina single-end reads shorter than ~70bp:

cd /srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/processing/align/bwa_aln_alignment
bwa_db=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/genome/bwa/stuv

for sample in ABGYY ATB1 ATB11 ATB13 ATB14 ATB15 ATB16 ATB17 ATB18 ATB19 ATB3 ATB4 ATB6 ATB7 ATB8 ATB9 BPBB BPYY BRYY BYPP GYRR HIST11 HIST12 HIST15 HIST2 HIST3 HIST4 HIST5 HIST6 HIST8 HIST9 MONK1 MONK11 MONK13 MONK15 MONK17 MONK2 MONK20 MONK22 MONK25 MONK27 MONK28 MONK4 MONK5 MONK6 MONK7 NOR1 NOR10 NOR11 NOR12 NOR13 NOR14 NOR15 NOR16 NOR17 NOR18 NOR20 NOR3 NOR4 NOR8 NOR9 PPYY PRBB PRRR PYBB PYRR PYYY sv057 sv058 sv059 sv060 sv061 sv062 sv063 sv064 sv065 sv066 sv087 sv088 sv089 sv090 sv091 V10181 V10531 V10731A;
do
echo WORKING WITH ${sample}
bwa aln -t 16 -B 5 $bwa_db ../../samples_batch/merged/${sample}.fq.gz > ${sample}.sai
bwa samse $bwa_db ${sample}.sai ../../samples_batch/merged/${sample}.fq.gz > ${sample}.sam
samtools view -bS ${sample}.sam | samtools sort -o ${sample}.bam
done
samtools flagstat ABGYY.bam
samtools flagstat ATB1.bam
b of 5:
ATB18: 2935527/4370005=67.17%  mapped
NOR17: 1787288/2733999=65.37% mapped

check reads mapped successfully: 

samtools flagstat HIST11.bam
samtools flagstat HIST12.bam
samtools flagstat HIST15.bam
samtools flagstat HIST2.bam
samtools flagstat HIST3.bam
samtools flagstat HIST4.bam
samtools flagstat HIST5.bam
samtools flagstat HIST6.bam
samtools flagstat HIST8.bam
samtools flagstat HIST9.bam 
HIST11: 146792/2349968=6.25% mapped
HIST12:275346/3076038=8.95% mapped
HIST8: 156867/2843194=5.52%% mapped
 

module load stacks/2.2 
gstacks -I ./ -M ../historic_populations.txt -O ./
populations -P ./ -M ../historic_populations.txt --vcf
 

Removed 0 loci that did not pass sample/population constraints from 415346 loci.
Kept 415346 loci, composed of 26722792 sites; 29 of those sites were filtered, 242682 variant sites remained.
    24592952 genomic sites, of which 1780173 were covered by multiple loci (7.2%).
Mean genotyped sites per locus: 64.33bp (stderr 0.00).

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  aw: 11.113 samples per locus; pi: 0.17193; all/variant/polymorphic sites: 14706012/231686/138996; private alleles: 23581
  hist: 4.781 samples per locus; pi: 0.12279; all/variant/polymorphic sites: 16028533/90003/33943; private alleles: 16147
  mv: 10.86 samples per locus; pi: 0.14802; all/variant/polymorphic sites: 18114133/225547/96574; private alleles: 4768
  mw: 10.991 samples per locus; pi: 0.16685; all/variant/polymorphic sites: 14708770/231599/126085; private alleles: 13185
  nc: 11.289 samples per locus; pi: 0.17285; all/variant/polymorphic sites: 14792920/234605/132768; private alleles: 15991
  or: 11.303 samples per locus; pi: 0.16254; all/variant/polymorphic sites: 14145933/232123/115637; private alleles: 7963
Populations is done.

 

 

Calling SNPS with filtering:
Strict filtering:

populations -P ./ -M ../historic_populations.txt --vcf -r 0.50 --min_maf 0.01 --write_random_snp -t 8
Removed 272960 loci that did not pass sample/population constraints from 415346 loci.
Kept 142386 loci, composed of 9215837 sites; 65847 of those sites were filtered, 75485 variant sites remained.
    9139366 genomic sites, of which 70462 were covered by multiple loci (0.8%).
Mean genotyped sites per locus: 64.72bp (stderr 0.00).

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  aw: 13.385 samples per locus; pi: 0.1581; all/variant/polymorphic sites: 8312464/66482/45995; private alleles: 3589
  hist: 7.8898 samples per locus; pi: 0.15151; all/variant/polymorphic sites: 6999644/12969/7498; private alleles: 1484
  mv: 13.345 samples per locus; pi: 0.14289; all/variant/polymorphic sites: 7022738/61475/33191; private alleles: 1052
  mw: 13.265 samples per locus; pi: 0.15851; all/variant/polymorphic sites: 8120833/65611/43631; private alleles: 2683
  nc: 13.336 samples per locus; pi: 0.16368; all/variant/polymorphic sites: 7993470/69360/46969; private alleles: 3696
  or: 13.394 samples per locus; pi: 0.15513; all/variant/polymorphic sites: 7877829/68163/41584; private alleles: 2129
Populations is done.

lenient filtering:

populations -P ./ -M ../historic_populations.txt --vcf --min_maf 0.01 --write_random_snp -t 8
Removed 0 loci that did not pass sample/population constraints from 415346 loci.
Kept 415346 loci, composed of 26722792 sites; 39995 of those sites were filtered, 107561 variant sites remained.
    24592952 genomic sites, of which 1780173 were covered by multiple loci (7.2%).
Mean genotyped sites per locus: 64.33bp (stderr 0.00).

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  aw: 10.028 samples per locus; pi: 0.20536; all/variant/polymorphic sites: 14708443/100640/64193; private alleles: 5414
  hist: 4.2583 samples per locus; pi: 0.15146; all/variant/polymorphic sites: 16104957/32828/12171; private alleles: 3377
  mv: 9.7637 samples per locus; pi: 0.18003; all/variant/polymorphic sites: 18119953/96664/48182; private alleles: 1791
  mw: 9.9017 samples per locus; pi: 0.20159; all/variant/polymorphic sites: 14711609/100408/61356; private alleles: 3862
  nc: 10.238 samples per locus; pi: 0.20998; all/variant/polymorphic sites: 14794677/102375/65561; private alleles: 5512
  or: 10.248 samples per locus; pi: 0.19904; all/variant/polymorphic sites: 14148468/100835/58591; private alleles: 3145
Populations is done.

 

mid filtering:

populations -P ./ -M ../historic_populations.txt --vcf -r 0.30 --min_maf 0.01 --write_random_snp -t 8
Removed 244064 loci that did not pass sample/population constraints from 415346 loci.
Kept 171282 loci, composed of 11070243 sites; 58528 of those sites were filtered, 89015 variant sites remained.
    10933716 genomic sites, of which 126562 were covered by multiple loci (1.2%).
Mean genotyped sites per locus: 64.63bp (stderr 0.00).

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  aw: 12.298 samples per locus; pi: 0.16771; all/variant/polymorphic sites: 9613384/77796/52617; private alleles: 4382
  hist: 6.6219 samples per locus; pi: 0.15842; all/variant/polymorphic sites: 9379639/18103/9472; private alleles: 1892
  mv: 12.195 samples per locus; pi: 0.14685; all/variant/polymorphic sites: 8682424/72833/38252; private alleles: 1459
  mw: 12.131 samples per locus; pi: 0.16676; all/variant/polymorphic sites: 9510690/77565/50254; private alleles: 3271
  nc: 12.234 samples per locus; pi: 0.17594; all/variant/polymorphic sites: 9476561/81578/54547; private alleles: 4600
  or: 12.319 samples per locus; pi: 0.16407; all/variant/polymorphic sites: 9198872/79815/47971; private alleles: 2686
Populations is done.

 

 

 

Bowtie + GATK
IF I RERUN REDO BAM FILES THEY WERE COPIED OVR BY BWA

https://www.biostars.org/p/81924/

Create a Bowtie genome database:

module load bowtie/2.3.5.1
cd /srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/genome
bowtie2-build -f Sturnus_vulgaris_2.3.1.simp.fasta Sturnus_vulgaris_2.3.1.simp
And GATK

module load java/8u121
module load samtools/1.10
module load picard/2.18.26
java -Xmx10g -jar /apps/picard/2.18.26/picard.jar CreateSequenceDictionary R=Sturnus_vulgaris_2.3.1.simp.fasta O=Sturnus_vulgaris_2.3.1.simp.dict
samtools faidx Sturnus_vulgaris_2.3.1.simp.fasta
 

module load gatk/4.1.0.0
cd /srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/processing/align/bowtie_GATK_alignment
export _JAVA_OPTIONS="-Xmx128g"
for i in ABGYY ATB1 ATB11 ATB13 ATB14 ATB15 ATB16 ATB17 ATB18 ATB19 ATB3 ATB4 ATB6 ATB7 ATB8 ATB9 BPBB BPYY BRYY BYPP GYRR HIST11 HIST12 HIST15 HIST2 HIST3 HIST4 HIST5 HIST6 HIST8 HIST9 MONK1 MONK11 MONK13 MONK15 MONK17 MONK2 MONK20 MONK22 MONK25 MONK27 MONK28 MONK4 MONK5 MONK6 MONK7 NOR1 NOR10 NOR11 NOR12 NOR13 NOR14 NOR15 NOR16 NOR17 NOR18 NOR20 NOR3 NOR4 NOR8 NOR9 PPYY PRBB PRRR PYBB PYRR PYYY sv057 sv058 sv059 sv060 sv061 sv062 sv063 sv064 sv065 sv066 sv087 sv088 sv089 sv090 sv091 V10181 V10531 V10731A;
do
echo working with $i;
GENOME=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/genome/Sturnus_vulgaris_2.3.1.simp
READS=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/processing/samples_batch/merged
bowtie2 -p 10 --phred33 --very-sensitive-local -x ${GENOME} -I 149 -X 900 --rg-id "$i" --rg SM:"$i" -U ${READS}/"$i".fq.gz -S "$i".sam
samtools view -bS ${i}.sam | samtools sort -o ${i}.bam
java -Xmx128g -jar /apps/picard/2.18.26/picard.jar MarkDuplicates INPUT="$i".bam OUTPUT="$i"_mark.bam METRICS_FILE="$i"_sort.metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
java -Xmx48g -jar /apps/picard/2.18.26/picard.jar BuildBamIndex I="$i"_mark.bam
gatk --java-options "-Xmx120G" HaplotypeCaller -R ${GENOME}.fasta -I ${i}_mark.bam -O ${i}.g.vcf -ERC GVCF --native-pair-hmm-threads 16;
done
 

GENOME=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/genome/Sturnus_vulgaris_2.3.1.simp.fasta
gatk --java-options "-Xmx120G" CombineGVCFs \
--reference $GENOME \
--variant ABGYY.g.vcf \
--variant ATB11.g.vcf \
--variant ATB13.g.vcf \
--variant ATB14.g.vcf \
--variant ATB15.g.vcf \
--variant ATB16.g.vcf \
--variant ATB17.g.vcf \
--variant ATB18.g.vcf \
--variant ATB19.g.vcf \
--variant ATB1.g.vcf \
--variant ATB3.g.vcf \
--variant ATB4.g.vcf \
--variant ATB6.g.vcf \
--variant ATB7.g.vcf \
--variant ATB8.g.vcf \
--variant ATB9.g.vcf \
--variant BPBB.g.vcf \
--variant BPYY.g.vcf \
--variant BRYY.g.vcf \
--variant BYPP.g.vcf \
--variant GYRR.g.vcf \
--variant HIST11.g.vcf \
--variant HIST12.g.vcf \
--variant HIST15.g.vcf \
--variant HIST2.g.vcf \
--variant HIST3.g.vcf \
--variant HIST4.g.vcf \
--variant HIST5.g.vcf \
--variant HIST6.g.vcf \
--variant HIST8.g.vcf \
--variant HIST9.g.vcf \
--variant MONK11.g.vcf \
--variant MONK13.g.vcf \
--variant MONK15.g.vcf \
--variant MONK17.g.vcf \
--variant MONK1.g.vcf \
--variant MONK20.g.vcf \
--variant MONK22.g.vcf \
--variant MONK25.g.vcf \
--variant MONK27.g.vcf \
--variant MONK28.g.vcf \
--variant MONK2.g.vcf \
--variant MONK4.g.vcf \
--variant MONK5.g.vcf \
--variant MONK6.g.vcf \
--variant MONK7.g.vcf \
--variant NOR10.g.vcf \
--variant NOR11.g.vcf \
--variant NOR12.g.vcf \
--variant NOR13.g.vcf \
--variant NOR14.g.vcf \
--variant NOR15.g.vcf \
--variant NOR16.g.vcf \
--variant NOR17.g.vcf \
--variant NOR18.g.vcf \
--variant NOR1.g.vcf \
--variant NOR20.g.vcf \
--variant NOR3.g.vcf \
--variant NOR4.g.vcf \
--variant NOR8.g.vcf \
--variant NOR9.g.vcf \
--variant PPYY.g.vcf \
--variant PRBB.g.vcf \
--variant PRRR.g.vcf \
--variant PYBB.g.vcf \
--variant PYRR.g.vcf \
--variant PYYY.g.vcf \
--variant sv057.g.vcf \
--variant sv058.g.vcf \
--variant sv059.g.vcf \
--variant sv060.g.vcf \
--variant sv061.g.vcf \
--variant sv062.g.vcf \
--variant sv063.g.vcf \
--variant sv064.g.vcf \
--variant sv065.g.vcf \
--variant sv066.g.vcf \
--variant sv087.g.vcf \
--variant sv088.g.vcf \
--variant sv089.g.vcf \
--variant sv090.g.vcf \
--variant sv091.g.vcf \
--variant V10181.g.vcf \
--variant V10531.g.vcf \
--variant V10731A.g.vcf \
--output historical_GATK_variantscombined.vcf
 

gatk --java-options "-Xmx120g" GenotypeGVCFs \
-R $GENOME \
-V historical_GATK_variantscombined.vcf \
-O historical_GATK_variantsgenotyped.vcf
 

Only Hist:

GENOME=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/genome/Sturnus_vulgaris_2.3.1.simp.fasta
gatk --java-options "-Xmx120G" CombineGVCFs \
--reference $GENOME \
--variant HIST11.g.vcf \
--variant HIST12.g.vcf \
--variant HIST15.g.vcf \
--variant HIST2.g.vcf \
--variant HIST3.g.vcf \
--variant HIST4.g.vcf \
--variant HIST5.g.vcf \
--variant HIST6.g.vcf \
--variant HIST8.g.vcf \
--variant HIST9.g.vcf \
--output historical_GATK_variantscombined_histonly.vcf
 

gatk --java-options "-Xmx120g" GenotypeGVCFs \
-R $GENOME \
-V historical_GATK_variantscombined_histonly.vcf \
-O historical_GATK_variantsgenotyped_histonly.vcf

 

 

 

 

BWA-mem and GATK
module load gatk/4.1.0.0
cd /srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/processing/align/
cp bwa_mem_alignment/*.bam bwamem_GATK_alignment/
cd bwamem_GATK_alignment
export _JAVA_OPTIONS="-Xmx128g"
module load samtools/1.10
module load java/8u121
module load gatk/4.1.9.0
for i in ABGYY ATB1 ATB11 ATB13 ATB14 ATB15 ATB16 ATB17 ATB18 ATB19 ATB3 ATB4 ATB6 ATB7 ATB8 ATB9 BPBB BPYY BRYY BYPP GYRR HIST11 HIST12 HIST15 HIST2 HIST3 HIST4 HIST5 HIST6 HIST8 HIST9 MONK1 MONK11 MONK13 MONK15 MONK17 MONK2 MONK20 MONK22 MONK25 MONK27 MONK28 MONK4 MONK5 MONK6 MONK7 NOR1 NOR10 NOR11 NOR12 NOR13 NOR14 NOR15 NOR16 NOR17 NOR18 NOR20 NOR3 NOR4 NOR8 NOR9 PPYY PRBB PRRR PYBB PYRR PYYY sv057 sv058 sv059 sv060 sv061 sv062 sv063 sv064 sv065 sv066 sv087 sv088 sv089 sv090 sv091 V10181 V10531 V10731A;
do
echo working with $i;
GENOME=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/genome/Sturnus_vulgaris_2.3.1.simp
samtools sort ${i}.bam -o ${i}_sort.bam
java -Xmx48g -jar /apps/picard/2.18.26/picard.jar AddOrReplaceReadGroups INPUT=${i}_sort.bam OUTPUT=${i}_sortRG.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=######## RGSM=${i}_sort_mark
java -Xmx128g -jar /apps/picard/2.18.26/picard.jar MarkDuplicates INPUT=${i}_sortRG.bam OUTPUT=${i}_sort_mark.bam METRICS_FILE=${i}_sort.metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
java -Xmx48g -jar /apps/picard/2.18.26/picard.jar BuildBamIndex I=${i}_sort_mark.bam
gatk --java-options "-Xmx4g" HaplotypeCaller -R ${GENOME}.fasta -I ${i}_sort_mark.bam -O ${i}.g.vcf -ERC GVCF --native-pair-hmm-threads 16
done
for i in NOR1;
do
echo working with $i;
GENOME=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/genome/Sturnus_vulgaris_2.3.1.simp
samtools sort ${i}.bam -o ${i}_sort.bam
java -Xmx48g -jar /apps/picard/2.18.26/picard.jar AddOrReplaceReadGroups INPUT=${i}_sort.bam OUTPUT=${i}_sortRG.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=######## RGSM=${i}_sort_mark
java -Xmx128g -jar /apps/picard/2.18.26/picard.jar MarkDuplicates INPUT=${i}_sortRG.bam OUTPUT=${i}_sort_mark.bam METRICS_FILE=${i}_sort.metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
java -Xmx48g -jar /apps/picard/2.18.26/picard.jar BuildBamIndex I=${i}_sort_mark.bam
gatk --java-options "-Xmx4g" HaplotypeCaller -R ${GENOME}.fasta -I ${i}_sort_mark.bam -O ${i}.g.vcf -ERC GVCF --native-pair-hmm-threads 16
done
 
GENOME=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/genome/Sturnus_vulgaris_2.3.1.simp.fasta
gatk --java-options "-Xmx120G" CombineGVCFs \
--reference $GENOME \
--variant ABGYY.g.vcf \
--variant ATB11.g.vcf \
--variant ATB13.g.vcf \
--variant ATB14.g.vcf \
--variant ATB15.g.vcf \
--variant ATB16.g.vcf \
--variant ATB17.g.vcf \
--variant ATB18.g.vcf \
--variant ATB19.g.vcf \
--variant ATB1.g.vcf \
--variant ATB3.g.vcf \
--variant ATB4.g.vcf \
--variant ATB6.g.vcf \
--variant ATB7.g.vcf \
--variant ATB8.g.vcf \
--variant ATB9.g.vcf \
--variant BPBB.g.vcf \
--variant BPYY.g.vcf \
--variant BRYY.g.vcf \
--variant BYPP.g.vcf \
--variant GYRR.g.vcf \
--variant HIST11.g.vcf \
--variant HIST12.g.vcf \
--variant HIST15.g.vcf \
--variant HIST2.g.vcf \
--variant HIST3.g.vcf \
--variant HIST4.g.vcf \
--variant HIST5.g.vcf \
--variant HIST6.g.vcf \
--variant HIST8.g.vcf \
--variant HIST9.g.vcf \
--variant MONK11.g.vcf \
--variant MONK13.g.vcf \
--variant MONK15.g.vcf \
--variant MONK17.g.vcf \
--variant MONK1.g.vcf \
--variant MONK20.g.vcf \
--variant MONK22.g.vcf \
--variant MONK25.g.vcf \
--variant MONK27.g.vcf \
--variant MONK28.g.vcf \
--variant MONK2.g.vcf \
--variant MONK4.g.vcf \
--variant MONK5.g.vcf \
--variant MONK6.g.vcf \
--variant MONK7.g.vcf \
--variant NOR10.g.vcf \
--variant NOR11.g.vcf \
--variant NOR12.g.vcf \
--variant NOR13.g.vcf \
--variant NOR14.g.vcf \
--variant NOR15.g.vcf \
--variant NOR16.g.vcf \
--variant NOR17.g.vcf \
--variant NOR18.g.vcf \
--variant NOR1.g.vcf \
--variant NOR20.g.vcf \
--variant NOR3.g.vcf \
--variant NOR4.g.vcf \
--variant NOR8.g.vcf \
--variant NOR9.g.vcf \
--variant PPYY.g.vcf \
--variant PRBB.g.vcf \
--variant PRRR.g.vcf \
--variant PYBB.g.vcf \
--variant PYRR.g.vcf \
--variant PYYY.g.vcf \
--variant sv057.g.vcf \
--variant sv058.g.vcf \
--variant sv059.g.vcf \
--variant sv060.g.vcf \
--variant sv061.g.vcf \
--variant sv062.g.vcf \
--variant sv063.g.vcf \
--variant sv064.g.vcf \
--variant sv065.g.vcf \
--variant sv066.g.vcf \
--variant sv087.g.vcf \
--variant sv088.g.vcf \
--variant sv089.g.vcf \
--variant sv090.g.vcf \
--variant sv091.g.vcf \
--variant V10181.g.vcf \
--variant V10531.g.vcf \
--variant V10731A.g.vcf \
--output historical_BWAmem_GATK_variantscombined.vcf
 

gatk --java-options "-Xmx120g" GenotypeGVCFs \
-R $GENOME \
-V historical_BWAmem_GATK_variantscombined.vcf \
-O historical_BWAmem_GATK_variantsgenotyped.vcf
 

 

grep -v "#" historical_BWAmem_GATK_variantsgenotyped.vcf | wc -l
grep -v "#" populations.snps.vcf | wc -l
