#!/bin/bash
umask 0002
## Authors: Logan Voegtly <logan.j.voegtly.ctr@mail.mil> & Kyle Long <kyle.a.long8.ctr@mail.mil> & Gregory Rice <gregory.k.rice.ctr@mail.mil> & Katie Arnold <catherine.e.arnold13.civ@mail.mil>
## Last modified: August 5, 2021
version=1.3.0
## Update: 1.3.0 Includes support for ARTICv4 primer scheme
## Update: 1.2.0 Update how coverage is calculated for summary stats and optimize parameters; incorporate nomerge
## Update: 1.1.0 Add support for QIAseq DIRECT primers
## Update: 1.0.1 Bugfixes: add umask; cd to directory before checking for files; add touch file to indicate analysis is done
## Update: 1.0.0 Consolidated scripts for easier updating in the future; consolidate to use a single parser for summary stats; Option to not run align_trim for probes or shotgun data

sample=$1
threads=$2
# ARTIC or YouSeq primer sets which need to be subtracted
protocol=$3
# Only do R1 analysis; SE=Single end PE=paired end
end_type=$4
# Skips the merging data step
nomerge=$5

cd /mnt/data/

if [[ -d $sample ]]; then
	cd $sample
fi


USAGE="
VAIW_hcov.sh $version
USAGE: docker run -v /path/to/data:/mnt/data/ bdrdgenomics/viral_amplicon_illumina_workflow /usr/local/bin/VAIW/VAIW_hcov.sh <sample prefix> <# threads> <protocol [ie. ARTIC or YouSeq]> <SE [single end] or PE [paired end]> \n
\tSamples are expected to in either /path/to/data/<sample name> or /path/to/data directory \n
\tSamples are expected to be paired data with name format <sample name>_R1.fastq.gz <sample name>_R2.fastq.gz\n
"

# ARTIC V4 primers
if [[ $protocol == *ARTICv4* ]]; then
	reference_dir=/usr/local/bin/artic-ncov2019/primer_schemes/nCoV-2019/V4
	reference="$reference_dir"/SARS-CoV-2.reference.fasta
	primer_bed="$reference_dir"/SARS-CoV-2.scheme.bed
	reference_gff="$reference_dir"/MN908947.gff3
	reference_acc=MN908947_3
	min_length=80
# ARTIC V3 primers
elif [[ $protocol == *ARTIC* ]]; then
	reference_dir=/usr/local/bin/artic-ncov2019/primer_schemes/nCoV-2019/V3
	reference="$reference_dir"/nCoV-2019.reference.fasta
	primer_bed="$reference_dir"/nCoV-2019.bed
	reference_gff="$reference_dir"/MN908947.gff3
	reference_acc=MN908947_3
	min_length=80
# If protocol contains YouSeq
# NOTE: This has been updated to use the same reference accession number for ARTIC and YouSeq
elif [[ $protocol == *YouSeq* ]]; then
	reference_dir=/usr/local/bin/hCoV_YouSeq/
	reference="$reference_dir"/MN908947.3.fasta
	primer_bed="$reference_dir"/youseq_covid-19_primers_bed_NC_045512_v2.0.bed
	reference_gff="$reference_dir"/MN908947.3.gff3
	reference_acc=MN908947_3
	min_length=80
elif [[ $protocol == *QIAseq* ]]; then
	reference_dir=/usr/local/bin/hCoV_QIAseq/
	reference="$reference_dir"/MN908947.3.fasta
	primer_bed="$reference_dir"/QIAseqDIRECTSARSCoV2primersfinal_edited.bed
	reference_gff="$reference_dir"/MN908947.3.gff3
	reference_acc=MN908947_3
	min_length=80
	max_indel=300
# Probe based data or shotgun data (will not run align_trim)
elif [[ $protocol == *probes* || $protocol == *shotgun* ]]; then
	reference_dir=/usr/local/bin/hCoV_YouSeq/
	reference="$reference_dir"/MN908947.3.fasta
	primer_bed="null"
	reference_gff="$reference_dir"/MN908947.3.gff3
	reference_acc=MN908947_3
	min_length=50
else
	echo "<protocol> $protocol is not ARTIC or YouSeq"
	echo -e $USAGE
	exit 3
fi

# If end type is paired end
if [ $end_type == "PE" ]; then
	read1="$sample"_R1.fastq.gz
	read2="$sample"_R2.fastq.gz
	trim_read1="$sample"_trimmed_R1.fastq.gz
	trim_read2="$sample"_trimmed_R2.fastq.gz
	merged_reads="$sample"_trimmed_merged.fastq.gz
	unmerged_read1="$sample"_trimmed_merged_R1.fastq.gz
	unmerged_read2="$sample"_trimmed_merged_R2.fastq.gz
	# Check to ensure R1 and R2 files are present
	if [[ ! -f $read1 && ! -f $read2 ]]; then
		echo "<end type> $end_type"
		echo "<sample name> could not find $read1 or $read2"
		echo -e $USAGE
		exit 2
	fi
# If protocol is single end (R1 data or MinIon data?)
elif [ $end_type == "SE" ]; then
	read1="$sample"_R1.fastq.gz
	trim_read1="$sample"_trimmed_R1.fastq.gz
	bbmerge_stats="null"
	# Check if R1 reads are present
	if ! [[ -f $read1 ]]; then
		read1="$sample".fastq.gz
		# Check if Read 1 does not have the R1 designation
		if [[ ! -f $read1 ]]; then
			echo "<end type> $end_type"
			echo "<sample name> could not find $read1"
			echo -e $USAGE
			exit 2
		fi
	fi
else
	echo "<end type> $end_type is not PE or SE"
	echo -e $USAGE
	exit 4
fi


log="$sample".log

# Alignment Filtering
min_mapping_quality=30

# SNV and Consensus Calling
min_coverage=10
min_frequency=0.3



if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "<threads> $threads is not defined or not a number"
    echo -e $USAGE
    exit 1
fi

if ! [[ "$max_indel" =~ ^[0-9]+$ ]]; then
    max_indel=500
fi
if ! [[ "$min_length" =~ ^[0-9]+$ ]]; then
    min_length=50
fi

echo VAIW_hcov.sh $version > $log
echo Start processing $sample >> $log
echo Start processing $sample

# 	Run filtering with bbduk
echo Running bbduk.sh on $sample >> $log
echo Running bbduk.sh on $sample
suffix="_trimmed"
bbduk_stats="$sample""$suffix"_bbduk.stats.txt
# If paired end run bbduk and bbmerge
if [ $end_type == "PE" ]; then
	bbduk.sh in1=$read1 in2=$read2 out1="$trim_read1" out2="$trim_read2" minlength=$min_length trimpolya=15 qtrim=r trimq=20 maq=20 ow 1>>$log 2> $bbduk_stats
	if ! [[ -f $trim_read1 && -f $trim_read2 ]]; then
		echo "Error: Trimmed read files $trim_read1 and $trim_read2 do not exist." >> $log
		exit 2
	fi

	# Perform merging if nomerge is not set
	if ! [ $nomerge == "nomerge" ]; then

		# Run BBmerge
		echo Running bbmerge.sh on $sample >> $log
		echo Running bbmerge.sh on $sample
		suffix="$suffix"_merged
		bbmerge_stats="$sample""$suffix"_bbmerge.stats.txt
		bbmerge.sh in1="$trim_read1" in2="$trim_read2" out="$merged_reads" outu1="$unmerged_read1" outu2="$unmerged_read2" mininsert=$min_length ow 1>> $log 2> $bbmerge_stats

		# Check if merged reads exist
		if ! [[ -f $merged_reads ]]; then
			echo "Error: Merged read file $merged_reads does not exist." >> $log
			exit 2
		fi
	fi

# If Single End only run bbduk
else
	bbduk.sh in=$read1 out="$trim_read1" minlength=50 trimpolya=15 qtrim=r trimq=20 maq=20 ow 1>>$log 2> $bbduk_stats
	if ! [[ -f $trim_read1  ]]; then
		echo "Error: Trimmed read file $trim_read1 do not exist." >> $log
		exit 2
	# Set the variable merged reads to point to the single end trim reads for the rest of the analysis
	else
		merged_reads="$trim_read1"
	fi
fi

# BBmap
echo Running bbmap alignment on $sample >> $log
echo Running bbmap alignment on $sample
suffix="$suffix"_bbmap
bam_out="$sample""$suffix".bam
bbmap_covstats="$sample""$suffix".covstats.txt

if [ $nomerge == "nomerge" ]; then
	bbmap.sh in1="$trim_read1" in2="$trim_read2" t=$threads ref=$reference basecov=$sample$suffix.basecov.txt nodisk=t out=$bam_out interleaved=f mappedonly=t covstats=$bbmap_covstats statsfile="$sample""$suffix".stats.txt maxindel=$max_indel dellenfilter=$max_indel inslenfilter=$max_indel local=t sam=1.3 ow >> $log
else
	bbmap.sh in="$merged_reads" t=$threads ref=$reference basecov=$sample$suffix.basecov.txt nodisk=t out=$bam_out interleaved=f mappedonly=t covstats=$bbmap_covstats statsfile="$sample""$suffix".stats.txt maxindel=$max_indel dellenfilter=$max_indel inslenfilter=$max_indel local=t sam=1.3 ow >> $log
fi

if ! [[ -f $bam_out ]]; then
	echo "Error: Bam file $bam_out does not exist." >> $log
	exit 2
fi

# Sort bam file
bbmap_bam_stats="$bam_out".stats.txt
samtools stats  -@ $threads "$bam_out" | grep ^SN | cut -f 2- > $bbmap_bam_stats

bam_in=$bam_out
bam_out="$sample""$suffix".sorted.bam
samtools sort -@ $threads $bam_in -o $bam_out 2>> $log
samtools stats -@ $threads "$bam_out" | grep ^SN | cut -f 2- > "$bam_out".stats.txt
if ! [[ -f $bam_out ]]; then
	echo "Error: Bam file $bam_out does not exist." >> $log
	exit 2
fi

if [[ ! $protocol == *probes* && ! $protocol == *shotgun* ]]; then
	# Align Trim
	echo Running align_trim on $sample >> $log
	echo Running align_trim on $sample
	suffix="$suffix"_primertrimmed.rg
	bam_in=$bam_out
	bam_out="$sample""$suffix".sorted.bam
	align_trim --no-read-groups "$primer_bed"  < "$bam_in" 2> "$sample""$suffix".alignreport.er.txt | samtools sort -T "$sample""$suffix" - -o $bam_out
	samtools stats -@ $threads "$sample""$suffix".sorted.bam | grep ^SN | cut -f 2- > "$sample""$suffix".sorted.bam.stats.txt
	if ! [[ -f $bam_out ]]; then
		echo "Error: Bam file $bam_out does not exist." >> $log
		exit 2
	fi
fi

# Filtering
echo Running reformat to filter low quality alignments on $sample >> $log
echo Running reformat to filter low quality alignments on on $sample
suffix="$suffix"_filtered
bam_in=$bam_out
bam_out="$sample""$suffix".sorted.bam
reformat.sh in=$bam_in out=stdout.bam mappedonly=t sam=1.3 ow=t ref=$reference minmapq=$min_mapping_quality 2>> "$sample""$suffix"_reformat.stats.txt | samtools sort -T "$sample""$suffix" - -o $bam_out
bbmap_bam_stats="$sample""$suffix".sorted.bam.stats.txt
samtools stats -@ $threads "$sample""$suffix".sorted.bam | grep ^SN | cut -f 2- > "$sample""$suffix".sorted.bam.stats.txt


# Mpileup for iVar
echo Running samtools mpileup on $sample >> $log
echo Running samtools mpileup on $sample
mpileup_file="$sample""$suffix".vcf
samtools mpileup -A -d 50000 -B -Q 0 --reference $reference $bam_out > $mpileup_file 2>> $log

# iVar consensus calling
echo Running ivar consensus on $sample >> $log
echo Running ivar consensus on $sample
ivar consensus -m "$min_coverage" -t "$min_frequency" -p "$sample".consensus -n N < $mpileup_file >> $log
consensus_stats="$sample".consensus.stats.txt
stats.sh in="$sample".consensus.fa out=$consensus_stats format=2 overwrite=t

# iVar
echo Running ivar variants on $sample >> $log
echo Running ivar variants on $sample
# ivar variant calling
ivar_snv_prefix="$sample"_ivar.snv
ivar_snv_file="$ivar_snv_prefix".tsv
ivar_vcf_file="$ivar_snv_prefix".vcf
ivar variants -m "$min_coverage" -t "$min_frequency" -p $ivar_snv_prefix -r $reference -g "$reference_gff" < $mpileup_file >> $log
snv_passed="$ivar_snv_prefix".pass.tsv
head -n 1 $ivar_snv_file > $snv_passed
grep "TRUE" $ivar_snv_file >> $snv_passed && sed -i "s/^/$sample\t/" $snv_passed
python /usr/local/bin/VAIW/ivar_variants_to_vcf.py --pass_only $ivar_snv_file $ivar_vcf_file


# Parsing sample statistics for summarystats.tsv
echo Parsing summary statistics >> $log
echo Parsing summary statistics
python /usr/local/bin/VAIW/parse_stats_v2.py "$sample" "$bbduk_stats" "$bbmerge_stats" "$bbmap_bam_stats" "$consensus_stats" "$bbmap_covstats" "$mpileup_file" "$snv_passed" "$reference_gff" "$primer_bed" > "$sample".summarystats.tsv

# Other variant analysis


# iVar allvariants
echo "Running ivar variants on $sample for allvariants" >> $log
echo "Running ivar variants on $sample for allvariants"
# ivar variant calling
ivar_allvariants_snv_prefix="$sample"_ivar_allvariants.snv
ivar_allvariants_snv_file="$ivar_allvariants_snv_prefix".tsv
ivar_allvariants_vcf_file="$ivar_allvariants_snv_prefix".vcf
ivar variants -m "$min_coverage" -t .02 -p $ivar_allvariants_snv_prefix -r $reference -g "$reference_gff" < $mpileup_file >> $log
ivar_allvariants_snv_passed="$ivar_allvariants_snv_prefix".pass.tsv
head -n 1 $ivar_allvariants_snv_file > $ivar_allvariants_snv_passed
grep "TRUE" $ivar_allvariants_snv_file >> $ivar_allvariants_snv_passed && sed -i "s/^/$sample\t/" $ivar_allvariants_snv_passed
python /usr/local/bin/VAIW/ivar_variants_to_vcf.py --pass_only $ivar_allvariants_snv_file $ivar_allvariants_vcf_file

#gzip mpileup_file to reduce file size
gzip $mpileup_file


# lowfreq
echo Running lofreq indelqual and call on $sample >> $log
echo Running lofreq indelqual on call $sample
suffix="$suffix"_lofreq_indelqual
bam_in=$bam_out
bam_out="$sample""$suffix".sorted.bam
lofreq_vcf_file="$sample""$suffix".vcf
lofreq_filtered_vcf_file="$sample""$suffix"_filtered.vcf
lofreq indelqual --dindel --ref=$reference $bam_in | samtools sort -T "$sample""$suffix" - -o $bam_out
samtools index -@ $threads $bam_out
lofreq call-parallel --force-overwrite --min-cov $min_coverage --pp-threads $threads --ref $reference --call-indels --out $lofreq_vcf_file $bam_out
lofreq filter --in $lofreq_vcf_file --out $lofreq_filtered_vcf_file --af-min 0.02


# snpEff all variants
echo Running snpEff eff on allvariants vcf $sample >> $log
echo Running snpEff eff on allvariants vcf $sample
snpeff_allvariants_prefix="$sample"_snpeff_allvariants
snpeff_allvariants_vcf_file="$snpeff_allvariants_prefix".vcf
snpeff_allvariants_html_stats="$snpeff_allvariants_prefix"_stats.html
snpeff_allvariants_csv_stats="$snpeff_allvariants_prefix"_stats.csv
snpEff eff -csvStats $snpeff_allvariants_csv_stats -htmlStats $snpeff_allvariants_html_stats $reference_acc $ivar_allvariants_vcf_file > $snpeff_allvariants_vcf_file


# snpEff on lowfreq
echo Running snpEff eff on lofreq vcf $sample >> $log
echo Running snpEff eff on lofreq vcf $sample
snpeff_lofreq_prefix="$sample"_snpeff_lofreq
snpeff_lofreq_vcf_file="$snpeff_lofreq_prefix".vcf
snpeff_lofreq_html_stats="$snpeff_lofreq_prefix"_stats.html
snpeff_lofreq_csv_stats="$snpeff_lofreq_prefix"_stats.csv
snpEff eff -csvStats $snpeff_lofreq_csv_stats -htmlStats $snpeff_lofreq_html_stats $reference_acc $lofreq_filtered_vcf_file > $snpeff_lofreq_vcf_file

# running vcf comparisons
echo Running vcf_comparison.sh on $sample >> $log
echo Running vcf_comparison.sh on $sample
/usr/local/bin/VAIW/vcf_comparison.sh $sample


echo Done processing $sample >> $log
echo Done processing $sample
touch "$sample".VAIW.done

# Go back to the parent directory
if [[ -d ../$sample ]]; then
	cd ..
fi
