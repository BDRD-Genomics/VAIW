## New in v2:
### Version 2 of VAIW includes:
* Rewrite of the script to be more generic for different references
* Ability to analyze user provided references (fasta, gff) and amplicons (primer bed file)
* Use named arguments instead of positional arguments
* Added many arguments for more fine-tuned analysis
* Performed a major update to align_trim to better find primer ends
* More detailed logging
* More robust: it would randomly fail when running bbtools when it failed to allocate memory

## About VAIW
#### This is a workflow to analyze viral genomes prepared with an amplicon based library protocol and sequenced on the Illumina platform. 
#### VAIW natively supports [ARTIC](https://artic.network/ncov-2019) (v3, v4, and v4.1), [QIASeq DIRECT](https://www.qiagen.com/us/products/next-generation-sequencing/rna-sequencing/qiaseq-direct-sars-cov-2-kits/), and  [YouSeq](https://youseq.com) (v2) amplicon protocols for SARS-CoV-2.

#### The user is able to provide references (fasta, gff) and amplicons (primer bed file).

#### The code is available on GitHub [BDRD-Genomics/VAIW](https://github.com/BDRD-Genomics/VAIW) 

### The docker container is available on DockerHub [bdrdgenomics/viral_amplicon_illumina_workflow](https://hub.docker.com/r/bdrdgenomics/viral_amplicon_illumina_workflow)

# USAGE
```
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
# Note: All accessions need to match across user provided files
    --reference-fasta       User provided reference fasta file; Accession needs be the only thing in the headers
    --reference-gff         User provided reference gff file
    --primer-bed            User provided primer bed file; This is not required for shotgun or probe protocols

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
```

# Summary of workflows:
1. Read QC with bbduk
2. Merge paired reads with bbmerge
3. Read mapping of only joined reads with bbmap
4. Clip amplicon primer sequences with align_trim (from ARTIC workflow)
5. Samtools mpileup
6. Consensus calling with iVar consensus
7. SNV calling with iVar variants
8. Parse through all data and generate a summary stats file
9. Additional SNV analysis with lofreq and snpeff

# How to run Docker container: 

## Run the docker container on your data:
```
docker run --rm --memory $memory --cpus $threads -u "$user_id":"$group_id" -v "$host_data_dir":"$container_data_dir" bdrdgenomics/viral_amplicon_illumina_workflow /usr/local/bin/VAIW/VAIW.sh -1 $read1 -2 $read2 -o "$container_data_dir" -n $sample -t $threads -m $memory -p $protocol -r "SARS-CoV-2"
```
 It is recommended to run with the --rm flag to remove the container when it is done.

## Set environmental variables:

### Get user and group ID for running
```
user_id=$(id -u ${USER})
group_id=$(id -g ${USER})
```
### Set the number of threads and memory for the docker container
```
threads=<number of threads>
memory=<memory in m or g>
```
### set the full path to the directory with the read data on the host and in the container
```
host_data_dir=<full path to read data on the host>
container_data_dir=<full path to read data in container>
```
### Set the sample name and read names
```
sample=<sample name>
read1=<read 1 name>
read2=<read 2 name>
```

**Note sample is used as the prefix on files generated**

### Set the protocol type (ARTICv4.1, ARTICv4, ARTIC, QIAseq_DIRECT, or YouSeq) and reference ("SARS-CoV-2")
```
protocol=<protocol>
reference="SARS-CoV-2"
```

### Set the read end type PE (paired end) or SE (single end)
```
end_type=<end_type>
```

# Contact:
Regina Z. Cer MSc.*: ``` regina.z.cer.civ[at]health.mil```


\* Genomics and Bioinformatics Department, Biological Defense Research Directorate, Naval Medical Research Command

### Citation:
Parts of this workflow were adapted from [Utah DoH ARTIC/Illumina Bioinformatic Workflow](https://github.com/CDCgov/SARS-CoV-2_Sequencing/tree/master/protocols/BFX-UT_ARTIC_Illumina) workflow and the [ARTIC workflow](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html)


# Disclaimers: 

This work was supported/funded by work unit number WUN A1417 and P0013_20_AH_01.01

The views expressed in this software reflect the results of research conducted by the author and do not necessarily reflect the official policy or position of the Department of the Navy, Department of Defense, nor the United States Government.

This work was prepared as part of their official duties. Title 17 U.S.C. 105 provides that `copyright protection under this title is not available for any work of the United States Government.' Title 17 U.S.C. 101 defines a U.S. Government work as work prepared by a military service member or employee of the U.S. Government as part of that person's official duties.
