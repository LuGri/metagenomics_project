#This script should be run in an empty directory

#------------------#
#-01 Simulate Data-#
#------------------#

# Activate environment
conda activate omics 

mkdir -p data && cd data

#download genomes
ncbi-genome-download all --assembly-accessions GCA_000204565.1,GCA_000023385.1,GCA_000025525.1 -F fasta -s genbank
#GCA_000204565.1 Clostridium botulinum
#GCA_000023385.1 Geobacillus sp.
#GCA_000025525.1 #Methanocaldococcus sp.

mkdir -p ncbi
mv genbank ncbi

#concatenate all files into one multifasta file
find ncbi -name '*fna.gz' | xargs zcat > selected_ncbi_genomes.fasta

#abundance file pulled from github repo

# Simulate reads from selected genomes and abundance file
iss generate --genomes selected_ncbi_genomes.fasta --abundance_file selected_ncbi_genomes_abundance.txt --n_reads 2000000 --model MiSeq --output ncbi_lc --cpus 8

#Quality Control

# Set read names as variables
cd ..
IN=data/ncbi_lc_R1.fastq
IN2=data/ncbi_lc_R2.fastq

# Basic statistics
seqstats $IN
seqstats $IN2

# More elaborate statistics
seqtk fqchk $IN
seqtk fqchk $IN2

#FastQC
mkdir -p 01_fastqc_results
fastqc --outdir 01_fastqc_results $IN $IN2

#Preprocessing with bbduk
bbduk.sh in1=$IN in2=$IN2 out1=$IN.bbduk.fq out2=$IN2.bbduk.fq qtrim=r trimq=10 maq=10 minlen=100

#fastqc again
fastqc --outdir 01_fastqc_results $IN.bbduk.fq $IN2.bbduk.fq

#-----------------------#
#-02 Taxonomic Profiling#
#-----------------------#

#set symbolic links
ln -s data/ncbi_lc_R1.fastq.bbduk.fq read1.fq
ln -s data/ncbi_lc_R2.fastq.bbduk.fq read2.fq

# Prepare KrakenUniq database: This was already done, please don't do it again
#DBDIR=/mirror/kraken_db
#mkdir -p $DBDIR
# Download the taxonomy
#krakenuniq-download --db $DBDIR taxonomy
# Download complete bacterial and archaeal genomes from RefSeq (and mask low-complexity regions)
#krakenuniq-download --db $DBDIR --threads 10 --dust refseq/bacteria refseq/archaea
# Build database
#krakenuniq-build --db $DBDIR --kmer-len 31 --threads 10 --taxids-for-genomes --taxids-for-sequences

# Classify reads with KrakenUniq
DBDIR=/mirror/kraken_db
krakenuniq --db $DBDIR --threads 10 --report-file kraken_taxonomy_profile.txt --paired read1.fq read2.fq  > kraken_read_classification.tsv


# Classify reads with Metaxa
conda activate metaxa
OUT=metaxa_out
metaxa2 -1 read1.fq -2 read2.fq -g ssu --mode metagenome --plus T --cpu 8 --megablast T -o $OUT

#Filter for long alignments
CUTOFF=240
cat $OUT.taxonomy.txt | awk -F'\t' -v cutoff="$CUTOFF" '$4 > cutoff {print $0}' > $OUT.taxonomy.filtered.txt

#Visualization
conda activate omics
python /apps/scripts/lineage2taxid.cpython-36.pyc --mail valid-mail-address-is-required-for-ncbi-lookups --db /apps/etetoolkit/taxa.sqlite -i $OUT.taxonomy.filtered.txt -c2 > $OUT.taxonomy.filtered.taxids.txt
ktImportTaxonomy $OUT.taxonomy.filtered.taxids.txt -o $OUT.taxonomy.filtered.krona.html -q 1 -t 2 -s 5

#Classify reads with mOTUs
OUT=motus_taxonomy_profile.txt
motus profile -f read1.fq -r read2.fq -t 8 > $OUT

# Get only lines where abundance is > 0
cat $OUT | awk -F'\t' 'BEGIN {OFS = FS} /^#/ {print $0} {if ($2 > 0) print $0}' > $OUT.filtered



# Profile metagenome with MetaPhlAn
conda activate metaphlan 
OUT=metaphlan_taxonomy_profile.txt
metaphlan2.py read1.fq,read2.fq --bowtie2out metagenome.bowtie2.bz2 --nproc 8 --input_type fastq > $OUT

# Visualize with Krona
metaphlan2krona.py -p $OUT -k metaphlan_krona.out
ktImportText -o $OUT.krona.html metaphlan_krona.out

#organize outputs
mkdir 02_profiling
mkdir 02_profiling/metaxa
mv metaxa_out.* 02_profiling/metaxa

mkdir 02_profiling/krakenuniq
mv kraken_* 02_profiling/krakenuniq

mkdir 02_profiling/motus
mv motus_taxonomy_profile* 02_profiling/motus

mkdir 02_profiling/metaphlan
mv meta* 02_profiling/metaphlan/

#-------------#
#-03 Assembly-#
#-------------#

conda activate omics

# Assemble the metagenome using the Megahit assembler
OUT=megahit_out
megahit -1 read1.fq -2 read2.fq -o megahit_out

#Evaluate assembly with quast
IN=megahit_out/final.contigs.fa
quast -1 read1.fq -2 read2.fq $IN

#Calculate contig coverage

IN=megahit_out/final.contigs.fa

# Align reads back to contigs/scaffolds with the short read mapper bbmap/bbwrap.sh
bbwrap.sh ref=$IN in=read1.fq in2=read2.fq out=aln.sam.gz kfilter=22 subfilter=15 maxindel=80
# Output per contig coverage to cov.txt with bbmap/pileup.sh
pileup.sh in=aln.sam.gz out=cov.txt

# Get important contigs stats (length/coverage/GC)
cat cov.txt | awk -F'\t' 'BEGIN {OFS = FS} {print $1, $2, $3, $9}' | sed 's/ [^\t]*//' > assembly_stats.txt

# Shorter and low-coverage contigs are often of lower quality, so it's better to get rid of them.
# Filter by length (e.g. 500bp or 1kb) and coverage (e.g. 3x)
cat assembly_stats.txt | awk -F'\t' 'BEGIN {OFS = FS} {if ($2 > 3 && $3 > 500) print $0}' > assembly_stats.filtered.txt
cat assembly_stats.filtered.txt | cut -f1 | grep -v '^#ID' > assembly_stats.filtered.contignames
# seqtk subseq can filter a fasta file based on a list of sequence names
seqtk subseq $IN assembly_stats.filtered.contignames > $IN.filtered.fasta
# Compare the original and the filtered fasta files
seqstats $IN
seqstats $IN.filtered.fasta

#organize output
mkdir 03_assembly
mv megahit_out quast_results ref assembly* aln.sam.gz cov.txt 03_assembly/

#------------#
#-04 Binning-#
#------------#

IN=03_assembly/megahit_out/final.contigs.fa.filtered.fasta 

#Metabat

OUT=metabat_out/bin

# Use existing mapping files and convert them to bam-files
SAM=03_assembly/aln.sam.gz
zcat $SAM | samtools view -S -b > tmp.bam
samtools sort tmp.bam -o aln.sorted.bam && rm tmp.bam
samtools index aln.sorted.bam
mkdir -p bamdir
mv *bam *bam.bai bamdir

jgi_summarize_bam_contig_depths --outputDepth depth.txt --pairedContigs paired.txt --minContigLength 1000 --minContigDepth 2 bamdir/*.bam

metabat -i $IN -a depth.txt -o $OUT -v --saveTNF saved.tnf --saveDistance saved.dist

#checkM
conda activate checkm
IN=metabat_out/
OUT=checkm_results
checkm lineage_wf -t 20 -x fa $IN $OUT > checkm_summary.txt

#organize output
mkdir 04_binning
mv bamdir checkm_results checkm_summary.txt depth.txt maxbin_out metabat_out paired.txt 04_binning/
