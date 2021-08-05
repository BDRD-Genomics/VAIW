# Updated May 25, 2021 by Logan Voegtly
# Now the default script for generating summary stats for SE and PE read sets
# Fixed incorrect reporting of which CDS region the low coverage occurs
# Recalculated

import sys
import re
import os
import os.path
from os import path

#file arguments
sample = sys.argv[1]
bbduk_file = sys.argv[2]
bbmerge_file = sys.argv[3]
bbmap_file = sys.argv[4]
consensus_file = sys.argv[5]
covstats_file = sys.argv[6]
vcf_file = sys.argv[7]
snv_file = sys.argv[8]
gff_file = sys.argv[9]
bed_file = sys.argv[10]

#read bbduk file for number of raw/trimmed reads, and calculate the average raw/trimmed read length
if path.exists(bbduk_file) == True:
	output = os.popen('grep "Input:" ' + bbduk_file).read()
	reads = output.split("\t")[1].split(" ")[0]
	bases = output.split("\t")[3].split(" ")[0]
	avg_read_length = str(round(float(bases)/float(reads),2))

	output = os.popen('grep "Result:" ' + bbduk_file).read()
	trimmed_reads = output.split("\t")[1].split(" ")[0]
	trimmed_bases = output.split("\t")[2].split(" ")[0]
	trimmed_avg_read_length = str(round(float(trimmed_bases)/float(trimmed_reads),2))
else:
	avg_read_length = str(0)
	reads = str(0)
	trimmed_reads = str(0)
	trimmed_avg_read_length = str(0)


#read bbmerge file for the number of joined reads and the average insert size
if path.exists(bbmerge_file) == True:
	output = os.popen('grep "Joined:" ' + bbmerge_file).read()
	joined_reads = output.split("\t")[1].split(" ")[0]
	available_reads = joined_reads

	output = os.popen('grep "Avg Insert:" ' + bbmerge_file).read()
	avg_insert_size = output.split("\t")[1].split(" ")[0].strip()
else:
	avg_insert_size = 'N/A'
	joined_reads = 'N/A'
	available_reads = int(reads)


#read bbmap file for the number of joined reads mapped and then calculate the total percentage mapped
if path.exists(bbmap_file) == True:
	output = os.popen('grep "reads mapped:" ' + bbmap_file).read()
	reads_mapped = output.split("\t")[1].split(" ")[0].strip()
	percentage_reads_mapped = str(round(((float(reads_mapped)/float(available_reads))*100),2))
else:
	reads_mapped = str(0)
	percentage_reads_mapped = str(0)

if path.exists(consensus_file) == True:
	output = os.popen('grep "Main genome contig sequence total:" ' + consensus_file).read()
	consensus_length = output.split("\t")[1].split(" ")[0].strip()
	# LJV updated to include a count of Ns in the sequence
	output = os.popen('grep "Main genome scaffold sequence total:" ' + consensus_file).read()
	scaffold_length = int(output.split("\t")[1].split(" ")[0].strip())
	number_of_ns = str(scaffold_length - int(consensus_length))

else:
	consensus_length = str(0)


#read covstats file for average coverage
# if path.exists(covstats_file) == True:
# 	for line in open(covstats_file):
# 		if line.startswith("#ID"):
# 			continue
# 		else:
# 			avg_cov = line.split("\t")[1]
# else:
# 	avg_cov = str(0)

#this section collects minimum coverage areas
#collect primer regions from bed file
if path.exists(bed_file) == True:
	with open(bed_file) as bed:
		first = str(bed.readline() [1:]).split("\t")[2]
		last = (bed.readlines() [-1:])[0].split("\t")[1]
		first_pos = int(first) + 1
		last_pos = int(last)
		primer_positions = range(int(first_pos), int(last_pos)+1)
else:
	primer_positions = range(0, int(consensus_length))
	first_pos=1

#get positions and coverages from the vcf file
vcf_collection = []
pos_and_pos_cov = []
if path.exists(vcf_file) == True:
	for line in open(vcf_file):
		#for each line, collect the coverage and position
		position = int(line.split("\t")[1])
		position_cov =  int(line.split("\t")[3])
		#mpileup writes deletions as "*", so these will be subtracted from the coverage count
		# LJV removed due to real deletions being called 0 coverage
# 		del_count = int(line.split("\t")[4].count('*'))
# 		position_cov = position_cov - del_count
		pos_and_pos_cov = [position, position_cov]
		vcf_collection.append(pos_and_pos_cov)
else:
	first_low = str(0)
	min_cov_position = str(0)

#makes sure all primer covered positions are included in the list
count = 0
for i in primer_positions:
	try:
		if vcf_collection[count][0] != first_pos+count:
			vcf_collection.insert(count, [first_pos+count, 0])
	except IndexError as error:
		vcf_collection.insert(count, [first_pos+count, 0])
	count = count + 1

#sort list by coverage
collection_sorted_by_cov = sorted(vcf_collection,key=lambda l:l[1])

#collect positions that are the lowest coverage
lowest_cov_positions = []
count = 0
for i in collection_sorted_by_cov:
	if count == 0:
		first_low = i
		lowest_cov_positions.append(i)
	else:
		if i[1] == first_low[1]:
			lowest_cov_positions.append(i)
	count = count + 1

min_cov_position = ""
CDS_regions = []
primer_regions = []
full_line = ""
x=""
elementStart=""
elementStop=""

for i in lowest_cov_positions:
	if min_cov_position == "":
		min_cov_position = " /" + str(i[0])
	else:
		last_min_cov_position = min_cov_position.rsplit(" ", 1)[1]
		if "/" in last_min_cov_position: #if the last min cov position has a slash in it, that means it is the start of a range (i.e., it is position 100 in a range of minimum coverage from 100-110)
			last_min_cov_position_number = last_min_cov_position.split("/")[1]
			if i[0] == int(last_min_cov_position_number) + 1: #if current position is equal to 1 more than the last position added as minimum coverage, then this is a range and it is written as such; if not, it is written as the start of a new range
				min_cov_position = str(min_cov_position) + "- " + str(i[0])
			else:
				min_cov_position = str(min_cov_position) + ", /" + str(i[0])

		else: #no slash, so the last min cov position is not the start of a range
			last_min_cov_position_number = min_cov_position.rsplit(" ", 1)[1]
			if i[0] == int(last_min_cov_position_number) + 1:
				min_cov_position = str(min_cov_position.rsplit(" ",1)[0]) + "- " + str(i[0])
			else:
				min_cov_position = str(min_cov_position) + ", /" + str(i[0])

pattern = re.compile('(?:\-)+')
min_cov_position = re.sub(pattern,'-',min_cov_position)
min_cov_position = min_cov_position.replace(" ","").replace("/","").replace(",",", ")
min_cov_list = min_cov_position.split(", ")

#this section reformats the bed file to a list that puts the left and right sides of a primer in a single string
bed_count = 0
bed_list = []
combine = ""
side = ""
if path.exists(bed_file) == True:
	with open(bed_file) as bed_file:
		for line in bed_file:
			try:
				bed_count = bed_count + 1
				primer_name = line.split("\t")[3]
				if bed_count > 2 and side == "right":
					if "LEFT" in primer_name:
						combine = left.replace("\n","\t") + right
						bed_list.append(combine)
						combine = ""
						side = ""
				if "LEFT" in primer_name:
					left = line
					side = "left"
				else:
					right = line
					side = "right"
			except IndexError as error:
				pass

#this section matches each minimum coverage region to a CDS if applicable, and determines which primers it belongs too
for element in min_cov_list:
	CDS_regions = []
	primer_regions = []
	if "-" in element: #if it has a "-" in it, this means it is a range.
		elementStart = element.split("-")[0]
		elementStop = element.split("-")[1]
		low_cov_range = range(int(elementStart), int(elementStop)) #creates a range of a minimum coverage region to loop through
		for row in open(gff_file): #find out if low coverage position is in CDS
			try:
				if row.split("\t")[2]=="CDS":
					CDSstart = int(row.split("\t")[3])
					CDSstop = int(row.split("\t")[4])
					CDSname = row.split("\t")[8].split(";")[1].split("=")[1]
					for x in low_cov_range:
						if x >= CDSstart and x <= CDSstop:
							CDS_regions.append(CDSname)
			except IndexError as error:
				pass
		unique_CDS_regions = set(CDS_regions)
		unique_CDS_regions = list(unique_CDS_regions)
		if unique_CDS_regions == []:
			unique_CDS_regions = ["Noncoding"]


		for row in bed_list: #find out which primer responsible for this low coverage position
			try:
				primerStart = int(row.split("\t")[1])
				primerStop = int(row.split("\t")[8])
				primerName = row.split("\t")[3]
				for x in low_cov_range:
					if x >= primerStart and x <= primerStop:
						primer_regions.append(primerName.split("_L")[0])
			except IndexError as error:
				pass
		unique_primer_regions = set(primer_regions)
		unique_primer_regions = list(unique_primer_regions)
		if unique_primer_regions == []:
			unique_primer_regions = ["No primer"]
		full_line = full_line + elementStart + "-" + elementStop + " (" + ", ".join(unique_CDS_regions) + "; " + ", ".join(unique_primer_regions) + "), "

	else: #the minimum coverage does not have a "-" which means it is not a range, so it is treated as a singular position
		x = int(element)
		for row in open(gff_file):
			try:
				if row.split("\t")[2]=="CDS":
					CDSstart = int(row.split("\t")[3])
					CDSstop = int(row.split("\t")[4])
					CDSname = row.split("\t")[8].split(";")[1].split("=")[1]
					if x >= CDSstart and x <= CDSstop:
						CDS_regions.append(CDSname)
			except IndexError as error:
				pass
		unique_CDS_regions = set(CDS_regions)
		unique_CDS_regions = list(unique_CDS_regions)
		if unique_CDS_regions == []:
			unique_CDS_regions = ["Noncoding"]

		for row in bed_list:
			try:
				primerStart = int(row.split("\t")[1])
				primerStop = int(row.split("\t")[8])
				primerName = row.split("\t")[3]
				if x >= primerStart and x <= primerStop:
					primer_regions.append(primerName.split("_L")[0])
			except IndexError as error:
				pass
		unique_primer_regions = set(primer_regions)
		unique_primer_regions = list(unique_primer_regions)
		if unique_primer_regions == []:
			unique_primer_regions = ["No primer"]
		full_line = full_line + str(x) + " (" + ", ".join(unique_CDS_regions) + "; " + ", ".join(unique_primer_regions) + "), "

#get rid of last comma and space
full_line = full_line[:-2]

#read snv file file snv and synonymous snv count
snv_count = 0
snv_line_count = 0
synon_count = 0
if path.exists(snv_file) == True:
	for line in open(snv_file):
		snv_line_count += 1
		if snv_line_count == 1:
			continue
		else:
			snv_count += 1
			ref_aa = line.split("\t")[17]
			alt_aa = line.split("\t")[19].strip()
			if ref_aa != alt_aa:
				synon_count +=1
else:
	snv_count = str(0)
	synon_count = str(0)

# New calculation of average coverage
if avg_insert_size == "N/A":
	read_length = float(trimmed_avg_read_length)
else:
	read_length = float(avg_insert_size)
# No dividing by zero
if consensus_length == '0':
	avg_cov = '0'
else:
	# Calculate average coverage and only print 2 decimal points
	avg_cov = '%.2f' % ((read_length * int(reads_mapped))/int(consensus_length))

print("SAMPLE NAME\tNUMBER OF RAW READS\tAVERAGE RAW READ LENGTH\tNUMBER OF TRIMMED READS\tAVERAGE TRIMMED READ LENGTH\tNUMBER OF JOINED READS\tAVERAGE INSERT SIZE\tNUMBER OF JOINED READS MAPPED\tPERCENTAGE OF JOINED READS MAPPED\tCONSENSUS LENGTH\tAVERAGE COVERAGE\tMINIMUM COVERAGE OUTSIDE OF ENDS\tMINIMUM COVERAGE POSITIONS\tTOTAL NUMBER OF SNV\tTOTAL NUMBER OF NON-SYNONYMOUS SNVS\tNUMBER OF Ns IN CONSENSUS")
print(sample + "\t" + reads + "\t" + avg_read_length + "\t" + trimmed_reads + "\t" + trimmed_avg_read_length + "\t" + joined_reads + "\t" + avg_insert_size + "\t" + reads_mapped + "\t" + percentage_reads_mapped + "\t" + consensus_length + "\t" + str(avg_cov) + "\t" + str(first_low[1]) + "\t" + full_line + "\t" + str(snv_count) + "\t" + str(synon_count)+ "\t" + number_of_ns)
