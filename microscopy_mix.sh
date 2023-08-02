#BSUB -q normal
#BSUB -J 96
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 54

cd `pwd`
mkdir tmp
tmpdir=tmp
sample_name=$(basename `pwd`)
dropseq_root=/share/home/guoguoji/tools/Drop-seq_tools-2.5.1

##fastq --> bam
java -jar ${dropseq_root}/3rdParty/picard/picard.jar FastqToSam F1=H_R1.fastq.gz F2=H_R2.fastq.gz \
O=H.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name TMP_DIR=$tmpdir

######################Paired################################

#####-------------  Cell Barcode -------------
#####add RT barcode
#########tag barcodes1
${dropseq_root}/TagBamWithReadSequenceExtended \
BASE_RANGE=1-6:22-27:43-56:84-100 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=true DISCARD_READ=false TAG_NAME=BC NUM_BASES_BELOW_QUALITY=4 \
INPUT=H.bam OUTPUT=$tmpdir/H1.bam COMPRESSION_LEVEL=5

${dropseq_root}/TagBamWithReadSequenceExtended \
BASE_RANGE=1-4 BASE_QUALITY=10 BARCODED_READ=2 TAG_BARCODED_READ=false DISCARD_READ=true TAG_NAME=CQ NUM_BASES_BELOW_QUALITY=1 \
INPUT=$tmpdir/H1.bam OUTPUT=$tmpdir/H2.bam COMPRESSION_LEVEL=5  

#########FilterBAM
${dropseq_root}/FilterBam TAG_REJECT=XQ INPUT=$tmpdir/H2.bam OUTPUT=$tmpdir/H3.bam 

##############corrected bam for one mismatch
python /share/home/guoguoji/tools/index_microscopy/index_microscopy_correctBC_drop.py /share/home/guoguoji/tools/index_microscopy/ $tmpdir/H3.bam filtered.bam \
echo "correct sam files done"
#######################


java -Xmx100g -jar /share/home/guoguoji/tools/Drop-seq_tools-2.5.1/3rdParty/picard/picard.jar \
SamToFastq INPUT=filtered.bam FASTQ=R1.fastq  READ1_TRIM=83

# ## PolyATrimmer

cutadapt -a A{10} -j 0 -O 10 --minimum-length=15 -o R1_trim.fastq \
R1.fastq
/share/home/guoguoji/tools/seqtk/seqtk seq -Ar R1_trim.fastq > R1_polyA_trim_reverse.fastq 

# Alignment STAR
/share/home/guoguoji/tools/STAR-2.5.2a/source/STAR \
--genomeDir /share/home/guoguoji/tools/Human-Mouse-Merged-featurecounts-use/GenomeDir/ \
--readFilesIn R1_polyA_trim_reverse.fastq \
--outFileNamePrefix star \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 \
--limitOutSJcollapsed 5000000 \
--outSAMtype BAM Unsorted

# --clip5pNbases 12 \

## MergeBamAlignment
java -Xmx100g -jar ${dropseq_root}/3rdParty/picard/picard.jar MergeBamAlignment \
REFERENCE_SEQUENCE=/share/home/guoguoji/tools/Human-Mouse-Merged-featurecounts-use/GRCh38_GRCm39.dna.primary_assembly.fa \
UNMAPPED_BAM=filtered.bam \
ALIGNED_BAM=starAligned.out.bam \
OUTPUT=merged.bam \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false

 ${dropseq_root}/TagReadWithGeneFunction I=merged.bam \
 O=star_gene_exon_tagged.bam \
 ANNOTATIONS_FILE=/share/home/guoguoji/tools/Human-Mouse-Merged-featurecounts-use/GRCh38_GRCm39.105.sorted.gtf

 # FilterBAM
 ${dropseq_root}/FilterBam I=star_gene_exon_tagged.bam O=mouse.bam REF_SOFT_MATCHED_RETAINED=GRCm39
 ${dropseq_root}/FilterBam I=star_gene_exon_tagged.bam O=human.bam REF_SOFT_MATCHED_RETAINED=hg38


 ${dropseq_root}/DigitalExpression -m 8g \
 I=mouse.bam \
 CELL_BARCODE_TAG=XC \
 MOLECULAR_BARCODE_TAG=XM \
 O=mouse_dge.txt.gz \
 SUMMARY=mouse_dge.summary.txt \
 NUM_CORE_BARCODES=20000 \
 LOCUS_FUNCTION_LIST=INTRONIC \
 TMP_DIR=.

 ${dropseq_root}/DigitalExpression -m 8g \
 I=human.bam \
 CELL_BARCODE_TAG=XC \
 MOLECULAR_BARCODE_TAG=XM \
 O=human_dge.txt.gz \
 SUMMARY=human_dge.summary.txt \
 NUM_CORE_BARCODES=20000 \
 LOCUS_FUNCTION_LIST=INTRONIC \
 TMP_DIR=.

 ${dropseq_root}/DigitalExpression -m 8g \
 I=star_gene_exon_tagged.bam \
 CELL_BARCODE_TAG=XC \
 MOLECULAR_BARCODE_TAG=XM \
 O=_dge.txt.gz \
 SUMMARY=_dge.summary.txt \
 NUM_CORE_BARCODES=20000 \
 LOCUS_FUNCTION_LIST=INTRONIC \
 TMP_DIR=.

 #  Read Summary
 ${dropseq_root}/BamTagHistogram I=star_gene_exon_tagged.bam O=out_cell_readcounts.txt TAG=XC #&& rm -rf $tmpdir
 samtools view star_gene_exon_tagged.bam |awk '{for (i=1; i<=NF; ++i) {if($i ~ "^XF:Z:"){print $i}}}' >type.txt
 gzip type.txt
