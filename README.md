# Sv4_HistoricalStarlings
Scripts and notes related to manuscript entitled "Historical museum samples enable the examination of divergent and parallel evolution during invasion"

2021_11_20_notebook_10049: PopGen analysis

2021_11_20_notebook_10040: Outlier analysis

This page contains small vignette's on the following
<li>DArTseq SNP variant calling</li>


<h2> SNP variant calling for raw RRS (DArTseq) data</h2>

The following is a summary of some SNP calling pipelines that I have used to process raw DArTseq data, but they may be applied to other types of reduced representation sequencing (RRS) raw data.

Variant calling falls into a few primary steps:
<li>Mapping raw reads to a reference genome (or for refernce free variant calling forming the pseudogenome and calling references to this)</li>
<li>Calling the variant SNPs for each sample based on the mapped reads</li>
<li>Filtering the called SNPs, however I won't cover this here as it is quite project specific.</li><p>


For the first two steps of the process, you can mix and match mapping and variant calling programs, though some are best paired. The below code covers:
* Mapping with [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) mem and aln, then variant calling with Stacks [Gstacks](https://catchenlab.life.illinois.edu/stacks/comp/gstacks.php).
* Mapping with [bowtie](http://bowtie-bio.sourceforge.net/index.shtml), then variant calling with [GATK](https://gatk.broadinstitute.org/hc/en-us).

Ulitmately I used BWA aln as my mapping software of choice, as this is [recommended](https://github.com/lh3/bwa#type) for Illumina single-end reads shorter than ~70bp, which is good for DArT-seq.

## Software

Below is a list of the modules used (and versions).

<pre class="r"><code>module load stacks/2.2
module load fastqc/0.11.8
module load bwa/0.7.17
module load bowtie/2.3.5.1
module load java/8u121
module load samtools/1.10
module load picard/2.18.26
module load gatk/4.1.0.0
</code></pre>


## Data cleaning: Stacks

<pre class="r"><code>RAW_DIR=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/rawdata/batch2_split
OUTPUT_DIR=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/processing/samples_batch
BARCODE_DIR=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/processing/barcodes

process_radtags -p ${RAW_DIR}/batch_$i/ -o ${OUTPUT_DIR}/batch_$i/ -b ${BARCODE_DIR}/bardcode_named_$i.txt -r -c -q --renz_1 pstI --renz_2 sphI
</code></pre>


Run QC:
<pre class="r"><code>fastqc *.fq.gz --outdir=fastqc/
</code></pre>

Count reads:
<pre class="r"><code>zcat *.fq.gz | echo $((`wc -l`/4))
</code></pre>

## Mapping: BWA aln and mem

Place your genome somewhere sensible and index it for BWA.

<pre class="r"><code>genome_fa=Sturnus_vulgaris_2.3.1.simp.fasta
bwa index -p bwa/stuv $genome_fa &> bwa/bwa_index.oe
bwa_db=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/genome/bwa/stuv
</code></pre>

<h3>Aligning with bwa mem</h3>
Here I have my sample names as a list within the file 'sample_names.txt'.
Aligning with bwa mem:
<pre class="r"><code>for i in cat sample_names.txt;
do
sample=$i
echo WORKING WITH ${sample}
fq_file=../../samples_batch/merged/${sample}.fq.gz
bam_file=${sample}.bam
bwa_db=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/genome/bwa/stuv
bwa mem -t 16 -M $bwa_db $fq_file | samtools view -bS | samtools sort > $bam_file
done
</code></pre>

Check reads mapped successfully: 
<pre class="r"><code>
samtools flagstat SAMPLENAME.bam #sub out SAMPLENAME with a few sample names
</code></pre>

## Aligning with bwa aln

<pre class="r"><code>for sample in cat sample_names.txt;
do
echo WORKING WITH ${sample}
bwa aln -t 16 -B 5 $bwa_db ../../samples_batch/merged/${sample}.fq.gz > ${sample}.sai
bwa samse $bwa_db ${sample}.sai ../../samples_batch/merged/${sample}.fq.gz > ${sample}.sam
samtools view -bS ${sample}.sam | samtools sort -o ${sample}.bam
done
</code></pre>

Check reads mapped successfully

<pre class="r"><code>samtools flagstat ABGYY.bam
samtools flagstat ATB1.bam
</code></pre>

<i>ATB18: 2935527/4370005=67.17%  mapped<br>
NOR17: 1787288/2733999=65.37% mapped</i>

As you can see, my historical samples had much lower mapping rates:

<pre class="r"><code>samtools flagstat HIST11.bam
samtools flagstat HIST12.bam
samtools flagstat HIST8.bam
</code></pre>

<i>HIST11: 146792/2349968=6.25% mapped<br>
HIST12:275346/3076038=8.95% mapped<br>
HIST8: 156867/2843194=5.52%% mapped<br></i>


## Variant Calling: gstacks

For BWA mem and aln, the output is a [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) file. Next, I called variants using gstacks. This is not a necessity; once you have the BAM files you can use different variant calling methods.

<pre class="r"><code>gstacks -I ./ -M ../historic_populations.txt -O ./
</code></pre>

Now the variant calling is done, we can produce a [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file using Stacks [Populations](https://catchenlab.life.illinois.edu/stacks/comp/populations.php)
 
<pre class="r"><code>populations -P ./ -M ../historic_populations.txt --vcf
</code></pre>

Populations is nice because it gives you some summary info at the end of the variant calling process:

<i>Removed 0 loci that did not pass sample/population constraints from 412079 loci.<br>
Kept 412079 loci, composed of 28225204 sites; 20 of those sites were filtered, 243589 variant sites remained.<br><br>
    27345713 genomic sites, of which 846899 were covered by multiple loci (3.1%).<br>
Mean genotyped sites per locus: 68.47bp (stderr 0.00).<br>

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  aw: 11.272 samples per locus; pi: 0.17294; all/variant/polymorphic sites: 15515519/232371/142127; private alleles: 24046
  hist: 4.8985 samples per locus; pi: 0.12483; all/variant/polymorphic sites: 18008520/91625/34830; private alleles: 13284
  mv: 11.015 samples per locus; pi: 0.14832; all/variant/polymorphic sites: 19889591/226170/98305; private alleles: 4740
  mw: 11.146 samples per locus; pi: 0.16776; all/variant/polymorphic sites: 15527851/232281/128802; private alleles: 13279
  nc: 11.442 samples per locus; pi: 0.17365; all/variant/polymorphic sites: 15591084/235304/135530; private alleles: 16176
  or: 11.45 samples per locus; pi: 0.16313; all/variant/polymorphic sites: 14915003/232899/117912; private alleles: 7942</i>

You can even do some preliminary filtering on Stacks: 

<pre class="r"><code>populations -P ./ -M ../historic_populations.txt --vcf -r 0.50 --min_maf 0.01 --write_random_snp -t 8
</code></pre>

Though I prefer the filtering options on VCFtools personally, as there are more options. The flag 'write_random_snp' or 'write-single-snp' is useful to run here through as it keep only one SNP per locus.
 

## Alternate mapping and variant calling process: Bowtie + GATK

You will find that this process has very similar steps to the above BWA and Stacks processes. The overall workflow of SNP variant calling is usually quite similar.

Create a Bowtie genome database:

<pre class="r"><code>bowtie2-build -f Sturnus_vulgaris_2.3.1.simp.fasta Sturnus_vulgaris_2.3.1.simp
</code></pre>

You will also need to create a sequence dictionary for GATK.

<pre class="r"><code>java -Xmx10g -jar /apps/picard/2.18.26/picard.jar CreateSequenceDictionary R=Sturnus_vulgaris_2.3.1.simp.fasta O=Sturnus_vulgaris_2.3.1.simp.dict
samtools faidx Sturnus_vulgaris_2.3.1.simp.fasta
</code></pre>

The below covers a few steps: bowtie2 mapping; samtools processed a SAM into a sorted BAM file; picard marks duplicates; picard index; GATK haplotype caller. Each of these steps is run independently for each sample.

<pre class="r"><code>
for i in cat sample_names.txt;
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
</code></pre>

Finally, we combine the individually called haplotypes into one VCF file.

<pre class="r"><code>GENOME=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv4_Historic/genome/Sturnus_vulgaris_2.3.1.simp.fasta
gatk --java-options "-Xmx120G" CombineGVCFs \
--reference $GENOME \
--variant SAMPLE_A.g.vcf \
...
--variant SAMPLE_Z.g.vcf \
--output historical_GATK_variantscombined.vcf
</code></pre>

And finally perform joint genotyping on the samples

<pre class="r"><code>gatk --java-options "-Xmx120g" GenotypeGVCFs \
-R $GENOME \
-V historical_GATK_variantscombined.vcf \
-O historical_GATK_variantsgenotyped.vcf
</code></pre>

GATK doesn't produce any summaries like populations, however we can make a quick count of the number of SNPs using the below:

<pre class="r">grep -v "#" historical_BWAmem_GATK_variantsgenotyped.vcf | wc -l
</code></pre>
