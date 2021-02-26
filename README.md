#### This is a workflow to analyze viral genomes prepared with an amplicon based library protocol and sequenced on the Illumina platform. This currently supports [ARTIC](https://artic.network/ncov-2019) (v3) and [YouSeq](https://youseq.com) (v2) amplicon protocols.

#### Parts of this workflow were adapted from [Utah DoH ARTIC/Illumina Bioinformatic Workflow](https://github.com/CDCgov/SARS-CoV-2_Sequencing/tree/master/protocols/BFX-UT_ARTIC_Illumina) workflow and the [ARTIC workflow](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html)

#### The code is available on GitHub [BDRD-Genomics/VAIW](https://github.com/BDRD-Genomics/VAIW) 

### The docker container is available on DockerHub [bdrdgenomics/viral_amplicon_illumina_workflow](https://hub.docker.com/repository/docker/bdrdgenomics/viral_amplicon_illumina_workflow)

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

## Environmental variables:

### Get user and group ID for running
```
user_id=$(id -u ${USER})
group_id=$(id -g ${USER})
```
### Set the number of threads
```
threads=<number of threads>
```
### set the full path to where the read data
```
host_data_dir=<full path to read data>
```
### Set the sample name
```
sample=<sample name>
```
Samples are expected to be in either ```$host_data_dir/$sample``` or ```$host_data_dir``` directory 
Samples are expected to be paired data with name format:
``` "$sample"_R1.fastq.gz and "$sample"_R2.fastq.gz```

**Note the "_" between sample and R1.fastq.gz and R2.fastq.gz**

### Set the protocol type (ARTIC or YouSeq)
```
protocol=<protocol>
```

### Set the read end type PE (paired end) or SE (single end)
```
end_type=<end_type>
```

### It is recommended to run with the --rm flag to remove the container when it is done.

# To run the docker container on your data:
```
docker run --rm --cpus $threads -u "$user_id":"$group_id" -v $host_data_dir:/mnt/data bdrdgenomics/viral_amplicon_illumina_workflow /usr/bin/local/bin/VAIW/VAIW_hcov.sh $sample $threads $protocol $end_type
```

# Contact:
Logan Voegtly M.Bin.*: ``` logan.j.voegtly.ctr[at]mail.mil```

Kyle Long M.S.*: ```kyle.a.long8.ctr[at]mail.mil```

\* Leidos Inc, Genomics and Bioinformatics Department, Biological Defense Research Directorate, Naval Medical Research Center

# Disclaimers: 

This work was supported/funded by work unit number WUN A1417 and P0013_20_AH_01.01

The views expressed in this software reflect the results of research conducted by the author and do not necessarily reflect the official policy or position of the Department of the Navy, Department of Defense, nor the United States Government.

Kimberly Bishop-Lilly is a federal employee of the United States government. This work was prepared as part of her official duties. Title 17 U.S.C. 105 provides that `copyright protection under this title is not available for any work of the United States Government.' Title 17 U.S.C. 101 defines a U.S. Government work as work prepared by a military service member or employee of the U.S. Government as part of that person's official duties.