#!/bin/bash
umask 0002
# Specific for Docker image 2.0
source /usr/local/etc/profile.d/conda.sh

## Authors: Logan Voegtly <logan.j.voegtly.ctr@mail.mil> & Kyle Long <kyle.a.long8.ctr@mail.mil> & Gregory Rice <gregory.k.rice.ctr@mail.mil> & Katie Arnold <catherine.e.arnold13.civ@mail.mil>
## Last modified: December 7, 2021
version=v2.0.1
## Update: 2.0.1 Bugfixes
## Update: 2.0.0 Rewrite of the script to be more generic for references and use named arguments instead of positional arguments
## Update: 1.2.0 Update how coverage is calculated for summary stats and optimize parameters; incorporate nomerge
## Update: 1.1.0 Add support for QIAseq DIRECT primers
## Update: 1.0.1 Bugfixes: add umask; cd to directory before checking for files; add touch file to indicate analysis is done
## Update: 1.0.0 Consolidated scripts for easier updating in the future; consolidate to use a single parser for summary stats; Option to not run align_trim for probes or shotgun data

USAGE="VAIW.sh $version
USAGE:
docker run -v /path/to/data:/mnt/data/ bdrdgenomics/viral_amplicon_illumina_workflow /usr/local/bin/VAIW/VAIW.sh -1 <read 1> -2 <read 2> -p <protocol> -o <output directory>

Input files are expected to be fastq or fastq.gz format. The reads are expected to be in the provided output directory.

# File Input:
-1, --r1                    Read 1 file name
-2, --r2                    Read 2 file name
    --se                    Single end Read file, not compatible with paired read set
-o, --sample-directory      Path to output directory, VAIW expects reads are in this directory
-n, --sample-name           Sample name, to be used as prefix for output files

# Analysis modifications
    --nomerge               Run without merging paired reads

# Built in Reference and Library type
-r, --reference             Reference genome, current default: SARS-CoV-2 (others TBD)
-p, --protocol              Protocol used to generate library, ie. ARTIC, ARTICv4, QIAseq, YouSeq, shotgun, probes

# User provided Reference:
# Functionality to be added at a later date

# Read Processing:
-q, --read-quality          Minimum read quality, default: 20
    --min_length            Minimum read length, default: 50

# Read Mapping:
    --max-indel             Maximum indel length when doing read mapping, default: 500
    --min-mapping-quality   Minimum quality when filtering bam file, default: 30

# SNV and consensus calling:
    --min-coverage          Minimum read coverage for consensus calling and SNV analysis, default: 10
    --min-frequency         Minimum frequency to call consensus or SNV, default: 0.3
    --lf-min-coverage       Minimum read coverage for low frequency SNV analysis, default: 15
    --lf-min-frequency      Minimum frequency for low frequency SNV analysis, default=0.02

# System parameters
-t, --threads               Number of threads/CPU to use, default: 4
-m, --memory                Amount of memory to be used in m or g (ie. 100m, 1g), default: 1g
"

# Function to print to stdout and log
log=/dev/null
print_log () {
  echo -e "$1" >> "$log"
  echo -e "$1"
}

run_cmd () {
  command=$1; continue_on_error=$2
  print_log "Running: $command"
  if ! eval "$command"; then
    if [[ $continue_on_error == "continue" ]]; then
      print_log "Warning: Command failed. Continuing with analysis."
      return
    fi
    sleep 2
    print_log "\nError: Initial run had error rerunning: $command\n"
    if ! eval "$command"; then
      print_log "\nError: command failed after a second attempt. Stopping program.\n"
      exit 4
    fi
  fi
}

get_bam_stats () {
#  bam_file=$1; stats_file=$2
  run_cmd "samtools stats -@ $threads $1 | grep ^SN | cut -f 2- > $2"
}

# Special case for running ivar in its environment
run_ivar_cmd () {
  run_cmd "conda activate ivar; $1; conda deactivate"
}

# Special case for running lofreq in its environment
run_lofreq_cmd () {
  run_cmd "conda activate lofreq; $1; conda deactivate"
}

is_int () {
  identifier=$1; integer=$2
  if ! [[ "$integer" =~ ^[0-9]+$ ]]; then
    print_log "\nError: $identifier $integer is not defined or not an integer!\n"
    print_log "$USAGE"
    exit 1
fi
}

is_float () {
  identifier=$1; float=$2
  if ! [[ "$float" =~ ^0?\.[0-9]+$ ]]; then
    print_log "\nError: $identifier $float is not defined or not a float!\n"
    print_log "$USAGE"
    exit 1
fi
}

is_memory () {
  identifier=$1; memory=$2
  if ! [[ "$memory" =~ ^[0-9]+[mMgG]$ ]]; then
    print_log "\nError: $identifier $memory is not defined or not a number!\n"
    print_log "$USAGE"
    exit 1
fi
}

# TODO: Integrate User provided reference files: fasta genome, gff file, primer file

# gets all the args to print to log
command="$0 $*"

# Read QC
read_quality=20
min_length=50

# Alignment Filtering
max_indel=500
min_mapping_quality=30

# SNV and Consensus Calling
min_coverage=10
min_frequency=0.3
lf_min_coverage=15
lf_min_frequency=0.02

# System preferences
threads=4
memory="1g"

reference="SARS-CoV-2"
nomerge=0

sample_directory=/mnt/data/
# TODO: implement qscore, min_length, and other read quality metrics
# TODO: Implement variables for variant calling parameters

# Adapted from https://www.golinuxcloud.com/beginners-guide-to-use-script-arguments-in-bash-with-examples/
# Getting options from the command line
while [ -n "$1" ];do
   case "$1" in
        -h|--help)
          echo -e "$USAGE"
          exit;;
        -1|--r1)
          shift
          read1="$1";;
        -2|--r2)
          shift
          read2="$1";;
        --se)
          shift
          se_read="$1";;
        -o|--sample-directory)
          shift
          sample_directory="$1";;
        -n|--sample-name)
          shift
          sample="$1";;
        --nomerge)
          nomerge="nomerge";;
        -r|--reference)
          shift
          reference="$1";;
        -p|--protocol)
          shift
          protocol="$1";;
        -q|--read_quality)
          shift
          is_int "read_quality" "$1"
          read_quality="$1";;
        --min-length)
          shift
          is_int "min_length" "$1"
          min_length="$1";;
        --max-indel)
          shift
          is_int "max_indel" "$1"
          max_indel="$1";;
        --min-mapping-quality)
          shift
          is_int "min_mapping_quality" "$1"
          min_mapping_quality="$1";;
        --min-coverage)
          shift
          is_int "min_coverage" "$1"
          min_coverage="$1";;
        --min-frequency)
          shift
          is_float "min_frequency" "$1"
          min_frequency="$1";;
        --lf-min-coverage)
          shift
          is_int "lf_min_coverage" "$1"
          lf_min_coverage="$1";;
        --lf-min-frequency)
          shift
          is_float "lf_min_frequency" "$1"
          lf_min_frequency="$1";;
        -t|--threads)
          shift
          is_int "threads" "$1"
          threads="$1";;
        -m|--memory)
          shift
          is_memory "Memory" "$1"
          memory="$1";;
        *)
       echo "Error: Incorrect input provided $1"
       echo "Command: $command"
       echo -e "$USAGE"
       exit
   esac
shift
done

# Initialize log file
if [[ ! $sample ]]; then
  echo -e "$USAGE"
  echo "Command: $command"
  exit 2
fi


# Move into directory with samples
if ! cd $sample_directory; then
  echo "Error: no such directory $sample_directory."
  echo "Command: $command"
  echo "$USAGE"
  exit 2
fi

# Set up the log file
log="$sample".log
echo "" > "$log"
print_log "VAIW.sh $version"
print_log "Command: $command"
print_log "\n[ Start processing $sample ]"

print_log "\n[ Analyzing data in $sample_directory ]"

print_log "\n[ Setting Up Reference Data ]"
# Find reference files
if [[ $reference == "SARS-CoV-2" ]]; then
  print_log "Using reference: $reference"
  reference_dir=/usr/local/bin/VAIW/references/SARS_CoV_2
  reference_fasta="$reference_dir"/SARS_CoV_2.reference.fasta
  reference_gff="$reference_dir"/SARS_CoV_2.reference.gff3
  reference_acc=MN908947.3

  # ARTIC V4.1 primers
  if [[ $protocol == *ARTICv4.1* ]]; then
    primer_bed="$reference_dir"/SARS_CoV_2.primers_ARTICv4.1.bed
    min_length=80
  # ARTIC V4 primers
  if [[ $protocol == *ARTICv4* ]]; then
    primer_bed="$reference_dir"/SARS_CoV_2.primers_ARTICv4.bed
    min_length=80
  # ARTIC V3 primers
  elif [[ $protocol == *ARTIC* ]]; then
    primer_bed="$reference_dir"/SARS_CoV_2.primers_ARTICv3.bed
    min_length=80
  # If protocol contains YouSeq
  # NOTE: This has been updated to use the same reference accession number for ARTIC and YouSeq
  elif [[ $protocol == *YouSeq* ]]; then
    primer_bed="$reference_dir"/SARS_CoV_2.primers_YouSeqv2.bed
    min_length=80
  elif [[ $protocol == *QIAseq* ]]; then
    primer_bed="$reference_dir"/SARS_CoV_2.primers_QIAseqDIRECTv1.bed
    min_length=80
    max_indel=300
  # Probe based data or shotgun data (will not run align_trim)
  elif [[ $protocol == *probes* || $protocol == *shotgun* ]]; then
    primer_bed="null"
    min_length=50
  else
    print_log "Error: Protocol $protocol is not yet installed into VAIW Docker for reference $reference."
    print_log "$USAGE"
    exit 3
  fi
else
  print_log "Error: Reference $reference is not yet installed into VAIW Docker. Added functionality is planned for future releases."
  print_log "$USAGE"
  exit 3
fi
print_log "Using primer file: $primer_bed"

# TODO loop through reference file for accession numbers

# Check if required reference files exist
if [[ $primer_bed != "null" ]]; then
  if [[ ! -f $reference_fasta || ! -f $reference_gff || ! -f $primer_bed ]]; then
    print_log "Error: Could not find all required reference files:
    fasta: $reference_fasta,
    gff: $reference_gff,
    primer bed: $primer_bed"
    print_log "$USAGE"
    exit 2
  fi
elif [[ ! -f $reference_fasta || ! -f $reference_gff ]]; then
   print_log "Error: could not find all required reference files:
    fasta: $reference_fasta,
    gff: $reference_gff"
  print_log "$USAGE"
  exit 2
fi
# End checking reference files
print_log "\n[ Setting up input reads ]"
# If single end data is not provided (assumes paired read)
if [[ ! -f $se_read ]]; then
  print_log "Using paired end reads: $read1 and $read2"
  end_type="PE"
  # Set up names for downstream analysis
  trim_read1="$sample"_trimmed_R1.fastq.gz
  trim_read2="$sample"_trimmed_R2.fastq.gz
  merged_reads="$sample"_trimmed_merged.fastq.gz
  unmerged_read1="$sample"_trimmed_merged_R1.fastq.gz
  unmerged_read2="$sample"_trimmed_merged_R2.fastq.gz
  # Check to ensure R1 and R2 files are present
  if [[ ! -f $read1 && ! -f $read2 ]]; then
    print_log "Error: could not find $read1 or $read2 in sample directory $sample_directory"
    print_log "$USAGE"
    exit 2
  fi
# If provided single end reads
else
  print_log "Using single end reads: $se_read"
  end_type="SE"
  read1=$se_read
  trim_read1="$sample"_trimmed_SE.fastq.gz
  bbmerge_stats="null"
  # Check if R1 reads are present
  if ! [[ -f $read1 ]]; then
    print_log "Error: could not find $se_read in sample directory $sample_directory"
    print_log "$USAGE"
    exit 2
  fi
fi


# Start Main

# 	Run filtering with bbduk
print_log "\n[ Read QC with bbduk ]"
suffix="_trimmed"
bbduk_stats="$sample""$suffix"_bbduk.stats.txt
# If paired end run bbduk and bbmerge
if [ $end_type == "PE" ]; then
  # Get parameters for bbduk
  cmd_args=( in1="$read1" in2="$read2" out1="$trim_read1" out2="$trim_read2" minlength="$min_length" trimpolya=15 qtrim=r trimq="$read_quality" maq="$read_quality" ow t="$threads" -Xmx"$memory" )
  # Run bbduk using run command to attempt running the command twice if it should fail the first time
  run_cmd "bbduk.sh ${cmd_args[*]} 1>> $log 2> $bbduk_stats"
  # Check to ensure the trimmed reads are present
  if ! [[ -f $trim_read1 && -f $trim_read2 ]]; then
    print_log "Error: Could not find Bbduk output trimmed read files $trim_read1 and $trim_read2."
    exit 2
  fi

  # Perform merging if nomerge is not set
  if ! [ $nomerge == "nomerge" ]; then

    # Run BBmerge
    print_log "\n[ Merging reads with bbmerge ]"
    suffix="$suffix"_merged
    bbmerge_stats="$sample""$suffix"_bbmerge.stats.txt
    cmd_args=(in1="$trim_read1" in2="$trim_read2" out="$merged_reads" outu1="$unmerged_read1" outu2="$unmerged_read2" mininsert="$min_length" ow t="$threads" -Xmx"$memory")
    run_cmd "bbmerge.sh ${cmd_args[*]} 1>> $log 2> $bbmerge_stats"

    # Check if merged reads exist
    if ! [[ -f $merged_reads ]]; then
      print_log "Error: Could not find bbmerge output merged read file $merged_reads."
      exit 2
    fi
  fi

# If Single End only run bbduk
else
  cmd_args=( in="$read1" out="$trim_read1" minlength="$min_length" trimpolya=15 qtrim=r trimq="$read_quality" maq="$read_quality" ow t="$threads" -Xmx"$memory" )
  run_cmd "bbduk.sh ${cmd_args[*]} 1>> $log 2> $bbduk_stats"

  # Check if trimmed single end read exists
  if ! [[ -f $trim_read1  ]]; then
    print_log "Error: Could not find Bbduk output trimmed read file $trim_read1."
    exit 2
  # Set the variable merged reads to point to the single end trim reads for the rest of the analysis
  else
    merged_reads="$trim_read1"
  fi
fi

# BBmap
print_log "\n[ Aligning reads to reference with bbmap ]"
suffix="$suffix"_bbmap
bam_out="$sample""$suffix".bam
bbmap_covstats="$sample""$suffix".covstats.txt
bbmap_log="$sample""$suffix".bbmap.log

if [ $nomerge == "nomerge" ]; then
  cmd_args=( in1="$trim_read1" in2="$trim_read2" )
else
  cmd_args=(in="$merged_reads")
fi
  cmd_args+=(t="$threads" -Xmx"$memory" ref="$reference_fasta" basecov="$sample""$suffix".basecov.txt nodisk=t out="$bam_out" )
  cmd_args+=(interleaved=f mappedonly=t covstats="$bbmap_covstats" statsfile="$sample""$suffix".stats.txt maxindel="$max_indel" )
  cmd_args+=(dellenfilter="$max_indel" inslenfilter="$max_indel" local=t sam=1.3 ow )
  run_cmd "bbmap.sh ${cmd_args[*]} 1>> $log 2>$bbmap_log"

if ! [[ -f $bam_out ]]; then
  print_log "Error: Could not find bbmap output bam file $bam_out."
  exit 2
fi

# Sort bam file
bbmap_bam_stats="$bam_out".stats.txt
get_bam_stats "$bam_out" "$bbmap_bam_stats"

bam_in=$bam_out
bam_out="$sample""$suffix".sorted.bam
run_cmd "samtools sort -@ $threads $bam_in -o $bam_out 2>> $log"
bam_stats="$bam_out".stats.txt
get_bam_stats "$bam_out" "$bam_stats"

if ! [[ -f $bam_out ]]; then
  print "Error: Could not find samtools sort bam file $bam_out."
  exit 2
fi

# Do not run align_trim on non-amplicon samples (Yes I know this is designed for Amplicon data)
if [[ ! $protocol == *probes* && ! $protocol == *shotgun* ]]; then
  # Align Trim
  print_log "\n[ Trimming primer sequences with align_trim ]"
  suffix="$suffix"_primertrimmed.rg
  bam_in=$bam_out
  bam_out="$sample""$suffix".sorted.bam
  run_cmd "align_trim --no-read-groups $primer_bed  < $bam_in 2> $sample$suffix.alignreport.er.txt | samtools sort -T $sample$suffix - -o $bam_out"
  bam_stats="$sample""$suffix".sorted.bam.stats.txt
  get_bam_stats "$bam_out" "$bam_stats"
  if ! [[ -f $bam_out ]]; then
    print_log "Error: Could not find align_trim output bam file $bam_out."
    exit 2
  fi
fi

# Filtering
print_log "\n[ Final filtering of alignment using reformat.sh ]"
suffix="$suffix"_filtered
bam_in=$bam_out
bam_out="$sample""$suffix".sorted.bam
reformat_stats="$sample""$suffix"_reformat.stats.txt
cmd_args=(in=$bam_in out=stdout.bam mappedonly=t sam=1.3 ow=t ref=$reference_fasta minmapq=$min_mapping_quality)
run_cmd "reformat.sh ${cmd_args[*]} 2>> $reformat_stats | samtools sort -T $sample$suffix - -o $bam_out"
bam_stats="$sample""$suffix".sorted.bam.stats.txt
get_bam_stats "$bam_out" "$bam_stats"



# Mpileup for iVar
print_log "\n[ Samtools mpileup ]"
mpileup_file="$sample""$suffix".vcf
run_cmd "samtools mpileup -A -d 50000 -B -Q 0 --reference $reference_fasta $bam_out > $mpileup_file 2>> $log"

# iVar consensus calling
print_log "\n[ Calling consensus with ivar ]"
run_ivar_cmd "ivar consensus -m $min_coverage -t $min_frequency -p $sample.consensus -n N < $mpileup_file >> $log  && touch consensus.done &"

# iVar
print_log "\n[ Calling variants with ivar ]"
# ivar variant calling
ivar_snv_prefix="$sample"_ivar.snv
ivar_snv_file="$ivar_snv_prefix".tsv
ivar_vcf_file="$ivar_snv_prefix".vcf
run_ivar_cmd "ivar variants -m $min_coverage -t $min_frequency -p $ivar_snv_prefix -r $reference_fasta -g $reference_gff < $mpileup_file >> $log && touch ivar_snv.done &"

# Wait for ivar consensus and snv calls to be done (run in background)
sleep_count=0
until [[ -f consensus.done && -f ivar_snv.done ]]
do
  sleep 1
  let "sleep_count+=1"
  # After 5 minutes report status
  if [[ "$sleep_count" -eq 300 ]]; then
    print_log "Waiting on ivar consensus and snv calling."
    run_cmd "ls -l $ivar_snv_file $sample.consensus.fa"
  fi
done

consensus_stats="$sample".consensus.stats.txt
run_cmd "stats.sh in=$sample.consensus.fa out=$consensus_stats format=2 overwrite=t"

snv_passed="$ivar_snv_prefix".pass.tsv
run_cmd "head -n 1 $ivar_snv_file > $snv_passed"
run_cmd "grep 'TRUE' $ivar_snv_file >> $snv_passed && sed -i 's/^/$sample\t/' $snv_passed" "continue"
run_cmd "python /usr/local/bin/VAIW/ivar_variants_to_vcf.py --pass_only $ivar_snv_file $ivar_vcf_file" "continue"


# Parsing sample statistics for summarystats.tsv
print_log "\n[ Generating summarystats file ]"
run_cmd "python /usr/local/bin/VAIW/parse_stats_v2.py '$sample' '$bbduk_stats' '$bbmerge_stats' '$bbmap_bam_stats' '$consensus_stats' '$bbmap_covstats' '$mpileup_file' '$snv_passed' '$reference_gff' '$primer_bed' > '$sample'.summarystats.tsv" "continue"

# Other variant analysis

# iVar allvariants
print_log "\n[ Calling low frequency variants with ivar ]"
# ivar variant calling
ivar_allvariants_snv_prefix="$sample"_ivar_allvariants.snv
ivar_allvariants_snv_file="$ivar_allvariants_snv_prefix".tsv
ivar_allvariants_vcf_file="$ivar_allvariants_snv_prefix".vcf
run_ivar_cmd "ivar variants -m $lf_min_coverage -t $lf_min_frequency -p $ivar_allvariants_snv_prefix -r $reference_fasta -g $reference_gff < $mpileup_file >> $log"
ivar_allvariants_snv_passed="$ivar_allvariants_snv_prefix".pass.tsv
run_cmd "head -n 1 $ivar_allvariants_snv_file > $ivar_allvariants_snv_passed"
run_cmd "grep 'TRUE' $ivar_allvariants_snv_file >> $ivar_allvariants_snv_passed && sed -i 's/^/$sample\t/' $ivar_allvariants_snv_passed" "continue"
run_cmd "python /usr/local/bin/VAIW/ivar_variants_to_vcf.py --pass_only $ivar_allvariants_snv_file $ivar_allvariants_vcf_file"

#gzip mpileup_file to reduce file size
run_cmd "gzip -f $mpileup_file" "continue"


# lowfreq
print_log "\n[ Calling low frequency variants with lofreq ]"
suffix="$suffix"_lofreq_indelqual
bam_in=$bam_out
bam_out="$sample""$suffix".sorted.bam
lofreq_vcf_file="$sample""$suffix".vcf
lofreq_filtered_vcf_file="$sample""$suffix"_filtered.vcf
run_lofreq_cmd "lofreq indelqual --dindel --ref=$reference_fasta $bam_in | samtools sort -T $sample$suffix - -o $bam_out"
run_cmd "samtools index -@ $threads $bam_out"
run_lofreq_cmd "lofreq call-parallel --force-overwrite --min-cov $lf_min_coverage --pp-threads $threads --ref $reference_fasta --call-indels --out $lofreq_vcf_file $bam_out"
run_lofreq_cmd "lofreq filter --in $lofreq_vcf_file --out $lofreq_filtered_vcf_file --af-min $lf_min_frequency"


# snpEff all variants
print_log "\n[ Converting ivar low frequency SNVs from nt to aa using snpeff ]"
snpeff_allvariants_prefix="$sample"_snpeff_allvariants
snpeff_allvariants_vcf_file="$snpeff_allvariants_prefix".vcf
snpeff_allvariants_html_stats="$snpeff_allvariants_prefix"_stats.html
snpeff_allvariants_csv_stats="$snpeff_allvariants_prefix"_stats.csv
run_cmd "snpEff eff -csvStats $snpeff_allvariants_csv_stats -htmlStats $snpeff_allvariants_html_stats $reference_acc $ivar_allvariants_vcf_file > $snpeff_allvariants_vcf_file"


# snpEff on lowfreq
print_log "\n[ Converting lofreq low frequency SNVs from nt to aa using snpeff ]"
snpeff_lofreq_prefix="$sample"_snpeff_lofreq
snpeff_lofreq_vcf_file="$snpeff_lofreq_prefix".vcf
snpeff_lofreq_html_stats="$snpeff_lofreq_prefix"_stats.html
snpeff_lofreq_csv_stats="$snpeff_lofreq_prefix"_stats.csv
run_cmd "snpEff eff -csvStats $snpeff_lofreq_csv_stats -htmlStats $snpeff_lofreq_html_stats $reference_acc $lofreq_filtered_vcf_file > $snpeff_lofreq_vcf_file"

# running vcf comparisons
print_log "\n[ Comparing low frequency variant calls ]"
vdir=vcf_comparison_dir

run_cmd "bgzip -f '$sample'_snpeff_lofreq.vcf"
run_cmd "bgzip -f '$sample'_snpeff_allvariants.vcf"
run_cmd "bcftools index '$sample'_snpeff_lofreq.vcf.gz"
run_cmd "bcftools index '$sample'_snpeff_allvariants.vcf.gz"
run_cmd "bcftools isec -c none -p ${vdir} '$sample'_snpeff_lofreq.vcf.gz '$sample'_snpeff_allvariants.vcf.gz" 
run_cmd "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AF\t%QUAL\t%DP4\t%ANN\n' ${vdir}/0002.vcf > ${vdir}/'$sample'_comparison_lofreq.tsv"
run_cmd "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%REF_DP][\t%REF_QUAL][\t%ALT_DP][\t%ALT_QUAL][\t%ALT_FREQ][\t%ANN][\t%SAMPLE=%GT]\n' ${vdir}/0003.vcf > ${vdir}/'$sample'_comparison_iVar.tsv"
print_log "running awk"
awk -v sample="$sample" -i inplace -F'\t' -vOFS='\t' ' { gsub("," , "\t" , $8); gsub(",","\t",$9);gsub(",","\t",$10);gsub(",","\t",$11); print sample "\t" $0}' ${vdir}/"$sample"_comparison_lofreq.tsv 
awk -v sample="$sample" -i inplace '{print sample "\t" $0}' ${vdir}/"$sample"_comparison_iVar.tsv 
run_cmd "echo -e 'SAMPLE\tREF_GENOME\tPOS\tREF\tALT\tDP\tAF\tQUAL\tDP_REF_FW\tDP_REF_RV\tDP_ALT_FW\tDP_ALT_RV\tANN\tREF_GENOME\tPOS\tREF\tALT\tREF_DP\tREF_QUAL\tALT_DP\tALT_QUAL\tALT_FREQ\tANN\tSAMPLE' > ${vdir}/'$sample'_comparison_combined.tsv"
run_cmd "paste -d' ' ${vdir}/'$sample'_comparison_lofreq.tsv ${vdir}/'$sample'_comparison_iVar.tsv >> ${vdir}/'$sample'_comparison_combined.tsv"


#run_cmd "/usr/local/bin/VAIW/vcf_comparison.sh $sample"


print_log "\n[ Done processing $sample ]"
touch "$sample".VAIW.done

