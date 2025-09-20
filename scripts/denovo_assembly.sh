#!/bin/bash

SECONDS=0
# Set working directory
cd /mnt/f/hypnoidus_abbreviatus

#set -e # Exit on error

# Step 1: Quality Check on Raw Reads
mkdir -p raw_reads/fastq_report
mkdir -p raw_reads/multiqc_report

#fastqc raw_reads/*.fastq.gz -o raw_reads/fastq_report
multiqc -d -s raw_reads/fastq_report/*_fastqc.zip -o raw_reads/multiqc_report
echo "finished running quality check on raw reads!"


mkdir -p trimmed_reads
mkdir -p unique/fastq_report
mkdir -p unique/multiqc_report

# Step 2: Trim adapters and poor quality reads
ls raw_reads/*_R1.fastq.gz | parallel -j 3 "
    base_name=\$(basename {} _R1.fastq.gz)
    R1_read=raw_reads/\${base_name}_R1.fastq.gz
    R2_read=raw_reads/\${base_name}_R2.fastq.gz
    trim_R1=trimmed_reads/\${base_name}_R1.fastq.gz
    trim_R2=trimmed_reads/\${base_name}_R2.fastq.gz
    unique_R1=unique/\${base_name}_R1.fastq.gz
    unique_R2=unique/\${base_name}_R2.fastq.gz
    report_html=trimmed_reads/\${base_name}.html
    report_json=trimmed_reads/\${base_name}.json

    echo 'Processing: ' \$base_name

    # Trim adapters, poor quality reads & low complexity regions
    fastp --in1 \$R1_read --in2 \$R2_read \
          --out1 \$trim_R1 --out2 \$trim_R2 \
          --detect_adapter_for_pe \
          --low_complexity_filter \
          --length_required 30 \
          --qualified_quality_phred 20 \
          --n_base_limit 5 \
          --thread 8 \
          --html \$report_html --json \$report_json

    # Remove duplicated sequences
    seqkit rmdup -s -o \$unique_R1 \$trim_R1
    seqkit rmdup -s -o \$unique_R2 \$trim_R2
"

echo "Fastp succesful trimming adapter, poor quality reads & low complexity regions!!"
echo "Seqkit successfull removing duplicated sequences!!!"

rm -r trimmed_reads #Optional to create room for space in the hard drive

# Quality Check on Trimmed Reads
fastqc unique/*.fastq.gz -o unique/fastq_report
multiqc -d -s unique/fastq_report/*_fastqc.zip -o unique/multiqc_report

echo "Finished quality check on trimmed reads!!"


# Step 3: De novo Assembly
# Assemble Transcripts with Trinity

mkdir -p trinity_assembly

cd unique

# Run Trinity
Trinity --seqType fq --max_memory 120G \
        --left Habb_WWA_6_R1.fastq.gz,Habb_WWA_7_R1.fastq.gz,Habb_WWA_8_R1.fastq.gz \
        --right Habb_WWA_6_R2.fastq.gz,Habb_WWA_7_R2.fastq.gz,Habb_WWA_8_R2.fastq.gz \
        --CPU 16 \
        --output /mnt/f/hypnoidus_abbreviatus/trinity_assembly

cd ..

# Access Assembly Quality
TrinityStats.pl trinity_assembly/Trinity.fasta > stats_trinity.log
echo "Trinity finished running!"

# Assess the raw assembly completeness
busco -i trinity_assembly/Trinity.fasta -o busco_out_raw -m transcriptome -l insecta_odb10 --cpu 8


# Step 4: Filtering (Length, Coverage, Clustering)

# Rename the identifier (TRINITY_...) to Habb_XXXXX
# Generate a mapping between original IDs and new Habb_XXXXX IDs
grep "^>" trinity_assembly/Trinity.fasta | awk '{split($0,a," "); print a[1]}' | sed 's/^>//' \
| awk '{printf("Habb_%05d\t%s\n", NR, $1)}' > id_map.txt

# Rename IDs
awk '
  NR==FNR {map[$2]=$1; next}
  /^>/ {
    split($0, a, " ");
    sub(/^>/, "", a[1]);
    if (a[1] in map) {
      a[1] = map[a[1]];
      printf(">%s", a[1]);
      for (i=2; i<=NF; i++) printf(" %s", a[i]);
      print ""
    } else {
      print $0
    }
    next
  }
  {print}
' id_map.txt trinity_assembly/Trinity.fasta > Habb_renamed_contigs.fasta


# Cluster Redundant Contigs
cd-hit-est -i Habb_renamed_contigs.fasta -o Habb_clustered_contigs.fasta -c 0.90 -n 8

# Length Filtering: Filter Short Contigs >=250
seqkit seq -m 250 Habb_clustered_contigs.fasta > Habb_filtered_contigs.fasta



# Filter for coverage >=5
# Use Salmon to quantify expression levels, filter for coverage, and map contigs:

# Activate Conda environment 'salmon_env
# Note: Due to compatibility issue, you might want to create another conda environment to install tools that are conflicting. Once done you can follow this step:
# For Example: My Salmon is having compatibility issue with my MAIN conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate salmon_env

#### Index the transcriptome
salmon index -t Habb_filtered_contigs.fasta -i transcripts_index -k 25


# Create Directories
mkdir -p quants/WWA_6
mkdir -p quants/WWA_7
mkdir -p quants/WWA_8


# Read 1
salmon quant \
    -i transcripts_index \
    -l A \
    -1 unique/Habb_WWA_6_R1.fastq.gz \
    -2 unique/Habb_WWA_6_R2.fastq.gz \
    -p 8 \
    --validateMappings \
    -o quants/WWA_6

# Read 2
salmon quant \
    -i transcripts_index \
    -l A \
    -1 unique/Habb_WWA_7_R1.fastq.gz \
    -2 unique/Habb_WWA_7_R2.fastq.gz \
    -p 8 \
    --validateMappings \
    -o quants/WWA_7

# Read 3
salmon quant \
    -i transcripts_index \
    -l A \
    -1 unique/Habb_WWA_8_R1.fastq.gz \
    -2 unique/Habb_WWA_8_R2.fastq.gz \
    -p 8 \
    --validateMappings \
    -o quants/WWA_8

conda deactivate
echo "Salmon finished running!"



# Filter for significance transcripts (≥5) across all the samples
# Prepare files
paste <(tail -n +2 quants/WWA_6/quant.sf | awk '{print $1}') \
      <(tail -n +2 quants/WWA_6/aux_info/ambig_info.tsv | awk '{print $1}') > quants/WWA_6/combined.tsv

paste <(tail -n +2 quants/WWA_7/quant.sf | awk '{print $1}') \
      <(tail -n +2 quants/WWA_7/aux_info/ambig_info.tsv | awk '{print $1}') > quants/WWA_7/combined.tsv

paste <(tail -n +2 quants/WWA_8/quant.sf | awk '{print $1}') \
      <(tail -n +2 quants/WWA_8/aux_info/ambig_info.tsv | awk '{print $1}') > quants/WWA_8/combined.tsv

# Combine and filter (>=5)
join -1 1 -2 1 <(sort quants/WWA_6/combined.tsv) <(sort quants/WWA_7/combined.tsv) | \
join -1 1 -2 1 - <(sort quants/WWA_8/combined.tsv) | \
awk '{sum = $2 + $3 + $4; if (sum >= 5) print $1}' > significant_transcripts.txt

echo "Filtering complete. $(wc -l < significant_transcripts.txt) transcripts retained with UniqueCount >=5"

# Filter the FASTA for significant transcripts
seqkit grep -f significant_transcripts.txt Habb_filtered_contigs.fasta > Habb_significant_contigs.fasta

echo " Filtering Successfully!!!"



# Step 5: Assess Assembly Quality
mkdir -p quast_output

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate quast_env

quast.py Habb_significant_contigs.fasta -o quast_output

conda deactivate
echo "Assembly Quality Check Done!!"



# Step 6: Remove Contaminants, Vectors, and Adapters

mkdir -p hypnoidus
mkdir -p blast_results

# Part A

# Break the transcripts into chunks
seqkit split2 -p 8 -O hypnoidus Habb_significant_contigs.fasta
# Run BLAST using Parallel to save computation time
parallel --eta -j 8 "blastn -db /mnt/f/blastdb/nt \
    -query {} \
    -evalue 1e-5 \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames' \
    -max_target_seqs 5 \
    -num_threads 2 \
    -out blast_results/{/.}.txt" ::: hypnoidus/Habb_significant_contigs.part_*.fasta

# Combine the blast results into a single fine
cat blast_results/Habb_significant_contigs.part_*.txt > blast_results/Habb_blast.txt



# Pause for manual contamination filtering
echo "BLAST results generated. Clean contaminants using Jupyter Notebook."
echo "Press [Enter] to continue once contaminants are removed."

read
echo "Resuming the script..."

# Remove contaminated sequences from original FASTA file
seqkit grep -v -f blast_results/contaminants_ids.txt Habb_significant_contigs.fasta > Habb_high_transcripts.fasta # Remove contaminants from transcriptome
echo "Bacteria|Fungi|Viruses|Vectors contaminants removed from the transcripts!"


# Part B

# Filter for adapters/vectors using FCS-adapter from NCBI
mkdir outputdir
./run_fcsadaptor.sh --fasta-input Habb_high_transcripts.fasta --output-dir ./outputdir --euk


# Remove adapters and vectors from the transcriptome using the FCS-adapter output
cat Habb_high_transcripts.fasta | python3 ./fcs.py clean genome \
    --action-report ./outputdir/fcs_adaptor_report.txt \
    --output Habb_cleaned_transcripts.fasta \
    --contam-fasta-out contam.fasta

echo "Adapters removed!!!"


# Step 7: Validate and format

# Check for duplicates and broken formatting
grep "^>" Habb_cleaned_transcripts.fasta | sort | uniq -d
grep -B1 "^[^>ACGTN]" Habb_cleaned_transcripts.fasta

echo "Validation completed!!"

# Step 8: Prepare transcriptome for TSA submission

# Remove length & path from headers
sed -E 's/ len=[0-9]+ path=\[.*\]//g' Habb_cleaned_transcripts.fasta > Habb_final_transcripts.fasta

# Rename using TSA standard
cp Habb_final_transcripts.fasta Habb_transcriptome.fsa

# Download submission template and run tbl2asn
wget https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/attic/documentation/SubmissionTemplate.sbt

# Generate report
tbl2asn -i Habb_transcriptome.fsa \
        -t SubmissionTemplate.sbt \
        -V b \
        -M n \
        -J \
        -a s \
        -c f \
        -Z discrepancy.report

# Compress for NCBI submission
gzip -c Habb_transcriptome.fsa > Habb_transcriptome.fsa.gz

echo "Transcripts ready for NCBI submission. Check 'discrepancy.report' for any final adjustments."


cp Habb_transcriptome.fsa Habb_transcriptome.fasta


# Step 9: Assess the filtered assembly completenes
#busco -i Habb_transcriptome.fasta -o busco_out_filtered -m transcriptome -l /mnt/f/hypnoidus_abbreviatus/busco_downloads/lineages/insecta_odb10 --cpu 8
#echo "BUSCO completeness after filtering"


# Step 10: Obtain quantification metrics for downstream analysis
#Activate Conda environment 'salmon_env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate salmon_env


# Index the transcriptome
mkdir -p salmon_final && cd salmon_final
salmon index -t /mnt/f/hypnoidus_abbreviatus/Habb_transcriptome.fasta -i transcripts_index -k 25

# Create Directories
mkdir -p quants/WWA_6
mkdir -p quants/WWA_7
mkdir -p quants/WWA_8


# Read 1

salmon quant \
    -i transcripts_index \
    -l A \
    -1 /mnt/f/hypnoidus_abbreviatus/unique/Habb_WWA_6_R1.fastq.gz \
    -2 /mnt/f/hypnoidus_abbreviatus/unique/Habb_WWA_6_R2.fastq.gz \
    -p 8 \
    --validateMappings \
    -o quants/WWA_6

# Read 2
salmon quant \
    -i transcripts_index \
    -l A \
    -1 /mnt/f/hypnoidus_abbreviatus/unique/Habb_WWA_7_R1.fastq.gz \
    -2 /mnt/f/hypnoidus_abbreviatus/unique/Habb_WWA_7_R2.fastq.gz \
    -p 8 \
    --validateMappings \
    -o quants/WWA_7

# Read 3
salmon quant \
    -i transcripts_index \
    -l A \
    -1 /mnt/f/hypnoidus_abbreviatus/unique/Habb_WWA_8_R1.fastq.gz \
    -2 /mnt/f/hypnoidus_abbreviatus/unique/Habb_WWA_8_R2.fastq.gz \
    -p 8 \
    --validateMappings \
    -o quants/WWA_8

cd ..
conda deactivate

echo "Salmon finished running!"



#Step 11: Pairwise BLAST for top 3 annotated species
mkdir -p pairwise

cd pairwise

# Download the transcripts for the top species
# Photinus pyralis
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/802/855/GCF_008802855.1_Ppyr1.3/GCF_008802855.1_Ppyr1.3_rna.fna.gz
gunzip GCF_008802855.1_Ppyr1.3_rna.fna.gz
mv GCF_008802855.1_Ppyr1.3_rna.fna Photinus_pyralis.fna

# Agrilus planipennis
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/699/045/GCF_000699045.2_Apla_2.0/GCF_000699045.2_Apla_2.0_rna.fna.gz
gunzip GCF_000699045.2_Apla_2.0_rna.fna.gz
mv GCF_000699045.2_Apla_2.0_rna.fna Agrilus_planipennis.fna

# Onthophagus taurus
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/711/975/GCF_036711975.1_IU_Otau_3.0/GCF_036711975.1_IU_Otau_3.0_rna.fna.gz
gunzip GCF_036711975.1_IU_Otau_3.0_rna.fna.gz
mv GCF_036711975.1_IU_Otau_3.0_rna.fna Onthophagus_taurus.fna

# Tribolium castaneum
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/031/307/605/GCF_031307605.1_icTriCast1.1/GCF_031307605.1_icTriCast1.1_rna.fna.gz
gunzip GCF_031307605.1_icTriCast1.1_rna.fna.gz
mv GCF_031307605.1_icTriCast1.1_rna.fna Tribolium_castaneum.fna

cd ..


# Build BLAST databases
### Create a BLAST Database for ALL sequences downloaded
echo "Creating custom database for each sequence part..."

makeblast_dir="/mnt/f/hypnoidus_abbreviatus/pairwise/"

makeblastdb -in "${makeblast_dir}Photinus_pyralis.fna" -dbtype nucl -out Photinus_pyralis_db &
makeblastdb -in "${makeblast_dir}Agrilus_planipennis.fna" -dbtype nucl -out Agrilus_planipennis_db &
makeblastdb -in "${makeblast_dir}Onthophagus_taurus.fna" -dbtype nucl -out Onthophagus_taurus_db &
makeblastdb -in "${makeblast_dir}Tribolium_castaneum.fna" -dbtype nucl -out Tribolium_castaneum_db &

wait
echo "Finished creating custom database for each sequence part..."


# Run tblastx
tblastx_dir="/mnt/f/hypnoidus_abbreviatus/"

tblastx -query "${tblastx_dir}Habb_transcriptome.fasta" -db Photinus_pyralis_db -evalue 1e-5 -outfmt 6 -num_threads 4 -out "${makeblast_dir}Ppyri_blast.txt" &
tblastx -query "${tblastx_dir}Habb_transcriptome.fasta" -db Agrilus_planipennis_db -evalue 1e-5 -outfmt 6 -num_threads 4 -out "${makeblast_dir}Apla_blast.txt" &
tblastx -query "${tblastx_dir}Habb_transcriptome.fasta" -db Onthophagus_taurus_db -evalue 1e-5 -outfmt 6 -num_threads 4 -out "${makeblast_dir}Otau_blast.txt" &
tblastx -query "${tblastx_dir}Habb_transcriptome.fasta" -db Tribolium_castaneum_db -evalue 1e-5 -outfmt 6 -num_threads 4 -out "${makeblast_dir}Tcas_blast.txt" &

wait
echo "BLAST search (tblastx) completed successfully!"
echo "Transcriptome Successfully Completed !!!



# Functional Enrichment

# Create working directories
mkdir -p Transdecoder
mkdir -p eggnog_results

# Predict CDS/Prot using TransDecoder
mkdir -p Transdecoder # Create a working directory
TransDecoder.LongOrfs -t Habb_transcriptome.fasta --output_dir Transdecoder
TransDecoder.Predict -t Habb_transcriptome.fasta --output_dir Transdecoder

# Predict the functional enrichment using the .pep output from TransDecoder

# Step 1: Download PFAM + MMseqs2 additional features + hmmer taxID
# Note: Execute each commandline in the terminal and follow the prompt for effective download into eggnog_data directory

download_eggnog_data.py -P -M --data_dir eggnog_data
download_eggnog_data.py -H -d 33208 --data_dir eggnog_data # Download hmmer taxid database for Metazoa (33208)

# Step 2: Run EggNOG_Mapper
emapper.py -i Transdecoder/Habb_transcriptome.fasta.transdecoder.pep \
  --itype proteins \
  --cpu 16 \
  --data_dir /mnt/f/hypnoidus_abbreviatus/eggnog_data \
  --pfam_realign denovo \
  --output_dir eggnog_results \
  -o Habb_annot

echo "EGGNOG-Mapper for Hypnoidus abbreviatus completed"
echo "BLAST search (tblastx) completed successfully!"


echo "Transcriptome Successfully Completed !!!















duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
