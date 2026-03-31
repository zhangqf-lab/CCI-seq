#!bin/bash

# ==============================================================================
# Description: This script performs deduplication, cell barcode whitelisting, 
#              UMI extraction, adapter trimming, and final format conversion.
# ==============================================================================

# 1. Deduplication: Remove physical duplicates based on sequence content using seqkit.
function Deduplicate(){
     zcat $1 | seqkit rmdup -s -o $2
    }

# 2. Whitelist: Identify valid Cell Barcodes (CB) from the raw pool.
function Whitelist(){
    infile=$1
    outfile=$2
    umi_tools whitelist --stdin ${infile} \
                         --stdout ${outfile} \
                         --bc-pattern='(?P<cell_1>.{16})(?P<umi_1>.{12})' \ # Pattern: 16bp Cell Barcode + 12bp UMI.
                         --extract-method=regex \
                         --ed-above-threshold=correct \
                         --method=umis \
                         --plot-prefix=${log_path}/threshold \
                         --log ${log_path}/whitelist.out \
                         --error ${log_path}/whitelist.err
}

# 3. Extract: Move CB and UMI information from the sequence to the Read ID.
function Extract(){
    infile=$1
    outfile=$2
    whitelist=$3
    
    umi_tools extract --bc-pattern='(?P<cell_1>.{16})(?P<umi_1>.{12})' \
                       --extract-method=regex \
                       --stdin ${infile} \
                       --stdout ${outfile} \
                       --filter-cell-barcode \
                       --whitelist=${whitelist}
}

# 4. Cutadapt: Trim 5' and 3' adapters and filter by specific read length.
function Cutadapt(){
    cutadapt -g $3 \
              -a $4 \
              -e 0.1 \
              -O 5 \
              --discard-untrimmed \
              -n 2 \
              $1 \
              -o $2 \
              -m 63 \  
              -M 63 \  
              # -m 71 \ #for hgmm_mix data. An additional round of barcodes indicating species was added.
              # -M 71 \ #for hgmm_mix data
              -j 1

# 5. Processed: Convert processed FASTQ to a tab-delimited text format ({BC}\t{UMI}\t{conbinatorial index}).
function Processed(){ 
    zcat $1 | awk -F ' |_' 'FNR%4==1{print $2"\t"$3}' >> $2/CellID_1.fq
    zcat $1 | awk 'FNR%4==2{print $1}' >> $2/CellID_2.fq
    paste $2/CellID_1.fq $2/CellID_2.fq  >> $2/processed.txt
    rm $2/CellID_1.fq $2/CellID_2.fq
}


# Define directory paths
input_path=./data/   
output_path=./output/
export log_path=./log/

# Step A: Remove PCR/physical duplicates
Deduplicate ${input_path}/CellID.fq.gz ${input_path}/CellID_deduplicate.fq.gz

# Step B: Generate the Cell Barcode whitelist
Whitelist ${input_path}/CellID_deduplicate.fq.gz ${output_path}/whitelist.txt

# Step C: Extract tags and filter reads
Extract ${input_path}/CellID_deduplicate.fq.gz ${input_path}/CellID_extracted.fq.gz ${output_path}/whitelist.txt

# Step D: Trim adapters (e.g., PCR primers or internal spacers)
adapter1=CGGCCTTAAAGC
adapter2=AGATCGGAAGAG
Cutadapt ${input_path}/CellID_extracted.fq.gz ${input_path}/CellID_trimmed.fq.gz ${adapter1} ${adapter2}

# Step E: Final conversion to text format for downstream analysis
Processed ${input_path}/CellID_trimmed.fq.gz ${output_path}
