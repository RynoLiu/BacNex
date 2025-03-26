#!/usr/bin/env bash

function preprocess_help(){
    echo "
Usage: BacNex.sh preprocess [options]

Parameters:
    -i,    --input           The raw fastq directory.
    -o,    --outdir          The output directory.
    -t,    --threads         The number of threads.
    -db,   --database        The built database.
    -kd,   --kneaddata       Kneaddata options.
    -hum,  --humann          Humann options.
    -mpa,  --metaphlan       Metaphlan options.
    -u,    --unmatched       Contains unmatched sequence.
    -h,    --help            Display the help message.
"
}


# bacnex args array
BACNEX_ARGS=(-i --input -o --outdir -t --threads -db --database -kd --kneaddata -hum --humann -mpa --metaphlan -u --unmatched -h --help)

if [[ ${#} == 0  ]]; then
    preprocess_help
fi

# default parameters
threads=1
hum_options="--remove-temp-output"
unmatched=false

while [ ${#} -gt 0 ]; do
    error_message="Error: a value is needed for '$1'";
    case $1 in
        -i | --input ) seq_dir=${2:?$error_message}; shift 2;;
        -o | --outdir ) outdir=${2:?$error_message}; shift 2;;
        -t | --threads ) threads=${2:?$error_message}; shift 2;;
        -db | --database ) db_dir=${2:?$error_message}; shift 2;;
        -kd | --kneaddata )
            shift
            KD_options=""
            while [[ $# -gt 0 ]]; do
                if [[ " ${BACNEX_ARGS[@]} " =~ " $1 " ]]; then
                    break
                fi
                KD_options="$KD_options $1"
                shift
            done
            KD_options=$(echo "$KD_options" | sed 's/^ //') ;;
        -hum | --humann )
            shift
            hum_options=""
            while [[ $# -gt 0 ]]; do
                if [[ " ${BACNEX_ARGS[@]} " =~ " $1 " ]]; then
                    break
                fi
                hum_options="$hum_options $1"
                shift
            done
            hum_options=$(echo "$hum_options" | sed 's/^ //') ;;
        -mpa | --metaphlan )
            shift
            mpa_options=""
            while [[ $# -gt 0 ]]; do
                if [[ " ${BACNEX_ARGS[@]} " =~ " $1 " ]]; then
                    break
                fi
                mpa_options="$mpa_options $1"
                shift
            done
            mpa_options=$(echo "$mpa_options" | sed 's/^ //') ;;
        -u | --unmatched ) unmatched=true; shift 1;;
        -h | --help) preprocess_help; exit 0;;
        * ) echo "Unknown parameter passed: $1"; preprocess_help; exit 1 ;;
    esac
done

# default mpa_options
if [ -z $mpa_options ]; then
  mpa_options="--bowtie2db ${db_dir}/metaphlan"
else
  mpa_options="--bowtie2db ${db_dir}/metaphlan ${mpa_options}"
fi

# ====== kneaddata main ======
# KD_options string to array
IFS=" " read -r -a KD_options_array <<< "$KD_options"

exec > >(tee -a "${outdir}/preprocess.log") 2>&1
echo "$(date '+%Y-%m-%d %H:%M:%S') - Preprocess starts..."
echo "$(date '+%Y-%m-%d %H:%M:%S') - kneaddata processing..."
KD_out=$outdir/"kneaddata_out"
mkdir $KD_out
R1files=$(ls $seq_dir | awk '/_1/ {print $NF}' | sort)
R2files=$(ls $seq_dir | awk '/_2/ {print $NF}' | sort)
R1vec=()
for x in $R1files; do   R1vec+=($x); done
R2vec=()
for x in $R2files; do   R2vec+=($x); done
Len=${#R1vec[@]}

# for loop
for (( i=0; i<${Len}; i++ ))
do
  R1file=${R1vec[i]}
  R2file=${R2vec[i]}
  sample=${R2file%%_2*}
  outPathRun=$KD_out/$sample
  if [ -d "$outPathRun" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - kneaddata $sample is done."
  else
    kneaddata -i1 $seq_dir/$R1file -i2 $seq_dir/$R2file -o $outPathRun -t $threads --output-prefix $sample -db $db_dir/"kneaddata" "${KD_options_array[@]}"
    kneaddata_read_count_table --input $outPathRun --output $outPathRun/"kneaddata_table.tsv"
    cd $outPathRun
    #rm $outPathRun/*trimmed* $outPathRun/*removed* $outPathRun/*hg37dec*
    gzip $sample'_paired_1.fastq'
    gzip $sample'_paired_2.fastq'
    gzip $sample'_unmatched_1.fastq'
    gzip $sample'_unmatched_2.fastq'
  fi
done


# =====HUMAnN 3 main=====
echo "$(date '+%Y-%m-%d %H:%M:%S') - humann3 processing..."

# KD_options string to array
IFS=" " read -r -a hum_options_array <<< "$hum_options"

# database config
humann_config --update database_folders nucleotide "${db_dir}/humann/chocophlan"
humann_config --update database_folders protein "${db_dir}/humann/uniref"
humann_config --update database_folders utility_mapping "${db_dir}/humann/utility_mapping"

hum_out=$outdir/"humann_out"
tmp_dir=$hum_out/"tmp_files"
mkdir $hum_out
mkdir $tmp_dir
samples=$(ls -l $KD_out |awk '/^d/ {print $NF}')
for sample in $samples
do
  outPathRun=$hum_out/$sample
  # check existed samples
  if [ -d "$outPathRun" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - humann $sample is done."
  else
    # unmatced option
    if [ $unmatched == "on" ]; then
      zcat $KD_out/$sample/"${sample}_paired_1.fastq.gz" $KD_out/$sample/"${sample}_paired_2.fastq.gz" \
      $KD_out/$sample/"${sample}_unmatched_1.fastq.gz" $KD_out/$sample/"${sample}_unmatched_2.fastq.gz" > $tmp_dir/"${sample}_merge.fastq"
    else
      zcat $KD_out/$sample/"${sample}_paired_1.fastq.gz" $KD_out/$sample/"${sample}_paired_2.fastq.gz" > $tmp_dir/"${sample}_merge.fastq"
    fi

    # main
    inseqFile=$tmp_dir/$sample'_merge.fastq'
    humann --input $inseqFile --output $outPathRun --threads $threads --output-format "tsv" "${hum_options_array[@]}" --metaphlan-options "${mpa_options}"
    cd $outPathRun
    mv $sample'_merge_genefamilies.tsv' $sample'_genefamilies.tsv'
    mv $sample'_merge_pathabundance.tsv' $sample'_pathabundance.tsv'
    mv $sample'_merge_pathcoverage.tsv' $sample'_pathcoverage.tsv'
    humann_renorm_table --input $sample'_genefamilies.tsv' --output $sample'_genefamilies_relab.tsv' --units relab
    humann_regroup_table --input $sample'_genefamilies_relab.tsv' -g uniref90_ko --output ko_table.tsv
    humann_regroup_table --input $sample'_genefamilies_relab.tsv' -g uniref90_go --output go_table.tsv
    humann_regroup_table --input $sample'_genefamilies_relab.tsv' -g uniref90_eggnog --output cog_table.tsv
    rm $inseqFile
    echo "$(date '+%Y-%m-%d %H:%M:%S') - humann3 done!"
  fi
done


# merge hum output
merge_dir=$hum_out/"final"
mkdir $merge_dir

## pathway
echo "$(date '+%Y-%m-%d %H:%M:%S') - Merge pathway abundance table..."
for sample in $samples; do
    cp $hum_out/$sample/*_pathabundance.tsv $tmp_dir/$sample"_pw.tsv"
done
humann_join_tables --input $tmp_dir --output $merge_dir/"pathway.tsv"
rm -rf $tmp_dir/*

## COG
echo "$(date '+%Y-%m-%d %H:%M:%S') - Merge COG abundance table..."
for sample in $samples; do
    cp $hum_out/$sample/"cog_table.tsv" $tmp_dir/$sample"_tmp.tsv"
done
humann_join_tables --input $tmp_dir --output $merge_dir/"cog.tsv"
rm -rf $tmp_dir/*

## KO
echo "$(date '+%Y-%m-%d %H:%M:%S') - Merge KO abundance table..."
for sample in $samples; do
    cp $hum_out/$sample/"ko_table.tsv" $tmp_dir/$sample"_tmp.tsv"
done
humann_join_tables --input $tmp_dir --output $merge_dir/"ko.tsv"
rm -rf $tmp_dir/*

## GO
echo "$(date '+%Y-%m-%d %H:%M:%S') - Merge GO abundance table..."
for sample in $samples; do
    cp $hum_out/$sample/"go_table.tsv" $tmp_dir/$sample"_tmp.tsv"
done
humann_join_tables --input $tmp_dir --output $merge_dir/"go.tsv"
rm -r $tmp_dir

echo "$(date '+%Y-%m-%d %H:%M:%S') - Preprocess all done!"
