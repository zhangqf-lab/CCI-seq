#!bin/bash

function Deduplicate(){
     zcat $1 | seqkit rmdup -s -o $2
    }

function Whitelist(){
    infile=$1
    outfile=$2
    umi_tools whitelist --stdin ${infile} \
                         --stdout ${outfile} \
                         --bc-pattern='(?P<cell_1>.{16})(?P<umi_1>.{12})' \
                         --extract-method=regex \
                         --ed-above-threshold=correct \
                         --method=umis \
                         --plot-prefix=${log_path}/threshold \
                         --log ${log_path}/whitelist.out \
                         --error ${log_path}/whitelist.err
}


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

function Cutadapt(){
    cutadapt -g $3 \
              -a $4 \
              -e 0.1 \
              -O 5 \
              --discard-untrimmed \
              -n 2 \
              $1 \
              -o $2 \
              # -m 63 \  # for other dataset.
              # -M 63 \  # for other dataset.
              -m 71 \ #for hgmm_mix data. An additional round of barcodes indicating species was added.
              -M 71 \ #for hgmm_mix data
              -j 1
}


function Processed(){ 
    zcat $1 | awk -F ' |_' 'FNR%4==1{print $2"\t"$3}' >> $2/CellID_1.fq
    zcat $1 | awk 'FNR%4==2{print $1}' >> $2/CellID_2.fq
    paste $2/CellID_1.fq $2/CellID_2.fq  >> $2/processed.txt
    rm $2/CellID_1.fq $2/CellID_2.fq
}


input_path=./data/
output_path=./output/
export log_path=./log/

Deduplicate ${input_path}/CellID.fq.gz ${input_path}/CellID_deduplicate.fq.gz
Whitelist ${input_path}/CellID_deduplicate.fq.gz ${output_path}/whitelist.txt
Extract ${input_path}/CellID_deduplicate.fq.gz ${input_path}/CellID_extracted.fq.gz ${output_path}/whitelist.txt

adapter1=CGGCCTTAAAGC
adapter2=AGATCGGAAGAG
Cutadapt ${input_path}/CellID_extracted.fq.gz ${input_path}/CellID_trimmed.fq.gz ${adapter1} ${adapter2}

Processed ${input_path}/CellID_trimmed.fq.gz ${output_path}
