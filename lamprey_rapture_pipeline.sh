################################################################################
################################################################################
################################################################################
# Lamprey Rapture Pipeline
# 10-8-18
# Seth Smith
#
# This is the pipeline used to genotype the lamprey rapture data.
#
# I've tried to format most steps in the following way: Fist, write a command using 
# ls, grep, and awk -F that writes a list of library names or sample names. Then,
# pipe this list into a while loop that echos commands into a list. The list of commands
# serves as documentation of how each library/sample was processed. Finally, run the
# list of commands using parallel -j.  
#
################################################################################
################################################################################
################################################################################
#
# The lamprey mitome was accessed from GenBank: U11880.1
# The lamprey genome was accessed from GenBank: /genomes/all/GCA/002/833/325/GCA_002833325.1_Pmar_germline_1.0
# Downloaded on 10-8-18.  Scaffolds are not on LGs/chromosomes
# scaffolds were renamed, mitome and genome were catted, and fasta was normalized with picard
# java -jar ~/bin/picard/picard.jar NormalizeFasta I=lamprey_g_m.fa O=pmarinus_genomic.fa &

################################################################################
################################################################################
################################################################################

# Make a full set of genome indexes
samtools faidx pmarinus_genomic.fa &
java -jar ~/bin/picard/picard.jar CreateSequenceDictionary R=pmarinus_genomic.fa O=pmarinus_genomic.dict &
bwa index pmarinus_genomic.fa &

################################################################################
################################################################################
################################################################################

# 1) Let's rename these fastq files. 
# Novogene mixed up SL3 and SL5 so we need to correct this.
# Also, let's make the names a little more human readable....
mv SL1_USPD16089468-D701-AK1680_HNT7MBBXX_L3_1.fq.gz SL1_1.fq.gz &
mv SL1_USPD16089468-D701-AK1680_HNT7MBBXX_L3_2.fq.gz SL1_2.fq.gz &
mv SL2_USPD16089468-D702-AK1681_HNT7MBBXX_L3_1.fq.gz SL2_1.fq.gz &
mv SL2_USPD16089468-D702-AK1681_HNT7MBBXX_L3_2.fq.gz SL2_2.fq.gz &
mv SL3_USPD16089468-D705-AK1543_HNT7MBBXX_L3_1.fq.gz SL5_1.fq.gz &
mv SL3_USPD16089468-D705-AK1543_HNT7MBBXX_L3_2.fq.gz SL5_2.fq.gz &
mv SL4_USPD16089468-D704-AK1780_HNT7MBBXX_L3_1.fq.gz SL4_1.fq.gz &
mv SL4_USPD16089468-D704-AK1780_HNT7MBBXX_L3_2.fq.gz SL4_2.fq.gz &
mv SL5_USPD16089468-D703-AK1682_HNT7MBBXX_L3_1.fq.gz SL3_1.fq.gz &
mv SL5_USPD16089468-D703-AK1682_HNT7MBBXX_L3_2.fq.gz SL3_2.fq.gz &

################################################################################
################################################################################
################################################################################

# Run FastQC on everything
mkdir /scratch/seth/lamprey/fastqc
ls | grep -v "2.fq.gz" | awk -F "_1.fq.gz" '{print $1}' | sort | uniq | while read -r LINE; do
    echo "fastqc "$LINE"_1.fq.gz -o /scratch/seth/lamprey/fastqc" >> fastqc.txt
    echo "fastqc "$LINE"_2.fq.gz -o /scratch/seth/lamprey/fastqc" >> fastqc.txt
done ;
parallel < fastqc.txt

################################################################################
################################################################################
################################################################################
#flip reads into proper orientation

mkdir barcodes
#transfer stacks barcode files (with unique sample names) to the barcodes directory

ls | grep fq | grep -v "_2.fq.gz" | awk -F "_1.fq.gz" '{print $1}' | uniq | sort | while read -r LINE; do
    echo "perl ~/scripts/bRAD_flip_trim.pl /scratch/seth/lamprey/barcodes/"$LINE".barcodes <(gzip -d -c "$LINE"_1.fq.gz) <(gzip -d -c "$LINE"_2.fq.gz) "$LINE"_flipped.1.fq "$LINE"_flipped.2.fq" >> flip_commands.txt
done &
parallel < flip_commands.txt

################################################################################
################################################################################
################################################################################
# compress files - might require a lot of temporary storage
ls | grep flipped | grep -v ".2.fq" | awk -F ".1.fq" '{print $1}' | uniq | sort | while read -r LINE; do
    gzip "$LINE".1.fq &
    gzip "$LINE".2.fq &
done &

################################################################################
################################################################################
################################################################################

# Demultiplex using process_radtags
# source activate stacks

mkdir /scratch/seth/lamprey/demult
mkdir /scratch/seth/lamprey/logs

#
ls | grep flipped | grep -v ".2.fq.gz" | awk -F "_flipped.1.fq.gz" '{print $1}' | sort | uniq | while read -r PLATES; do
    echo "process_radtags -1 ./"$PLATES"_flipped.1.fq.gz -2 ./"$PLATES"_flipped.2.fq.gz -i gzfastq -y gzfastq -o /scratch/seth/lamprey/demult/ -b /scratch/seth/lamprey/barcodes/"$PLATES".barcodes --inline_null -e sbfI --barcode_dist_1 1 --retain_header" >> process_commands.txt &
done ;
 
less -S process_commands.txt
#
parallel < process_commands.txt &

################################################################################
################################################################################
################################################################################

# 5) remove clonal reads

cd /scratch/seth/lamprey/demult
mkdir /scratch/seth/lamprey/clone

ls | grep -v ".2.fq.gz" | awk -F ".1.fq.gz" '{print $1}' | sort | uniq | while read -r LINE; do
    echo "clone_filter -1 ./"$LINE".1.fq.gz -2 ./"$LINE".2.fq.gz -i gzfastq -o /scratch/seth/lamprey/clone 2> "$LINE".log &" >> clone_filter_commands.txt
done ;
#
less -S clone_filter_commands.txt
#
parallel -j 16 < clone_filter_commands.txt

################################################################################
################################################################################
################################################################################

# 6) Trim reads, remove adapters, toss reads that are too short to effectivley map

#You'll need a file called "adapters.fa" in your home directory

mkdir trim
cd /scratch/seth/lamprey/trim


for file in $(ls /scratch/seth/lamprey/clone | grep fq); do ln -s /scratch/seth/lamprey/clone/$file ;
done ;
cp ~/adapters.fa ./
ls | grep -v ".2.fq.gz" | grep -v log | grep ".1.fq.gz" | awk -F ".1.1.fq.gz" '{print $1}' | sort | uniq | while read -r LINE; do
    echo "trimmomatic PE -threads 1 -phred33 "$LINE".1.1.fq.gz "$LINE".2.2.fq.gz "$LINE"_paired.1.fq.gz "$LINE"_unpaired.1.fq.gz "$LINE"_paired.2.fq.gz "$LINE"_unpaired.2.fq.gz ILLUMINACLIP:adapters.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50 2> "$LINE".trim.log" >> trimmomatic_commands.txt
done ;

less -S trimmomatic_commands.txt 

parallel -j 16 < trimmomatic_commands.txt &

################################################################################
################################################################################
################################################################################

# 7) Map to the reference genome using BWA-mem

cd /scratch/seth/lamprey
mkdir map
cd /scratch/seth/lamprey/map

##Generate sample list and pipe into a command writing while-loop
ls /scratch/seth/lamprey/trim | grep ".1.fq.gz" | grep -v "unpaired" | grep "paired" | awk -F "_paired" '{print $1}' | sort | uniq | while read -r LINE; do
    echo "bwa mem -t 1 -v 1 -T 30 -h 5 -M -R '@RG\tID:"$LINE"\tSM:"$LINE"\tPL:ILLUMINA\tLB:LB1' /scratch/seth/Genomes/lamprey/pmarinus_genomic.fa <(gzip -d -c /scratch/seth/lamprey/trim/"$LINE"_paired.1.fq.gz) <(gzip -d -c /scratch/seth/lamprey/trim/"$LINE"_paired.2.fq.gz) | samtools view -Sb - > "$LINE".pe.bam" >> BWA_Commands.txt;
done ;

### Convert single quotes to double quotes because of quoting funkiness above 
cat BWA_Commands.txt | sed -e 's|['\'']|\"|g' > BWA_Commands.list

#" (do not run or worry about this. Its a place holder for the extra quote in the sed command)

rm BWA_Commands.txt
###Maintain 30 concurrent jobs
parallel -j 30 < BWA_Commands.list &

################################################################################
################################################################################
################################################################################

#sort the bam files
cd /scratch/seth/lamprey/map

ls | grep bam | awk -F "." '{print $1}' | sort | uniq | while read -r LINE; do
    echo "samtools sort "$LINE".pe.bam -o "$LINE".pe.sorted.bam" >> Sort_Commands.txt
done ;

parallel < Sort_Commands.txt &

################################################################################
################################################################################
################################################################################

# Remove reads that are not properly paired and with mq less than 30

ls | grep sorted.bam | awk -F "." '{print $1}' | sort | uniq | while read -r LINE; do
    echo "samtools view -q 30 -h -f2 -F2308  "$LINE".pe.sorted.bam | samtools view -Sb > "$LINE".pe.sorted.filtered.bam" >> read_filter.txt
done ;



################################################################################
################################################################################
################################################################################

# Move things around...
mkdir sorted
mkdir filtered

ls | grep pe.sorted.filtered.bam | awk -F "." '{print $1}' | sort | uniq | while read -r LINE; do
    mv ./"$LINE".pe.sorted.bam ./sorted/
    mv ./"$LINE".pe.sorted.filtered.bam ./filtered/
    rm "$LINE".pe.bam
done ;


################################################################################
################################################################################
################################################################################

# Figure out how many records exist in each BAM and reomve ones that drop out.
ls | grep pe.sorted.filtered.bam | awk -F "." '{print $1}' | sort | uniq | while read -r LINE; do
    paste <(echo "$LINE".pe.sorted.filtered.bam) <(samtools view -h "$LINE".pe.sorted.filtered.bam | grep -v "@" | wc -l) >> records.txt &
done 

less -S records.txt

# remove empty bam files
rm PM_17_DCJ_5_1.pe.sorted.filtered.bam


################################################################################
################################################################################
################################################################################
# Call genotypes with gstacks

mkdir /scratch/seth/lamprey/genos
ls | grep "sorted.filtered.bam" > bams.list
gstacks -I ./ -S .pe.sorted.filtered.bam --threads 24 -M /scratch/seth/lamprey/lamprey_popmap.txt -O /scratch/seth/lamprey/genos &

# output vcf using populations, then convert to proper vcf using custom script from seth.


################################################################################
################################################################################
################################################################################
#Call genos with freebayes
mkdir freebayes_calls

ls | grep "sorted.filtered.bam" | awk -F "." '{print $1}' | sort | uniq | while read -r LINE; do
    echo "samtools index "$LINE".pe.sorted.filtered.bam" >> Index_Commands.txt ;
done ;
parallel < Index_Commands.txt &

cat /scratch/seth/Genomes/lamprey/pmarinus_genomic.fa.fai | awk '{print $1":"0"-"$2}' | while read -r REGION; do
	echo "freebayes -L lamprey.list --region "$REGION" -f /scratch/seth/Genomes/lamprey/pmarinus_genomic.fa --binomial-obs-priors-off --no-population-priors --prob-contamination 0.01 --min-coverage 1000 --min-base-quality 20 -G 2 --strict-vcf --genotype-qualities | bcftools view -i 'INFO/AC > 1 && INFO/AF < 1 && INFO/AF > 0 && INFO/AB > 0 && INFO/AB < 1' - > ./freebayes_calls/"$REGION".freebayes.raw.vcf" >> freebayes.commands;
done

parallel -j 16 < 

ls | grep vcf | while read -r VCF; do
	cat "$VCF" | grep -v "#" >> temp.vcf ;
done &

cat `ls | grep vcf | head -n 1` | grep "#" > temp.vcf.header

cat temp.vcf.header temp.vcf > lamprey_freebayes.baseline.vcf &