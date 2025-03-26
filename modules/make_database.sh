#!/usr/bin/env bash

function make_database_help(){
    echo "
Usage: BacNex.sh make_database [options]

Parameters:
    -o,    --outdir          The output file directory.
    -h,    --help            Display the help message.
"
}

if [[ ${#} == 0  ]]; then
    make_database_help
fi

# parameters
while [ ${#} -gt 0 ]; do
    error_message="Error: a value is needed for '$1'";
    case $1 in
        -o | --outdir ) outdir=${2:?$error_message}; shift 2;;
        -h | --help) make_database_help; exit 0;;
        * ) echo "Unknown parameter passed: $1"; make_database_help; exit 1 ;;
    esac
done


# ===== prepare =====
echo "$(date '+%Y-%m-%d %H:%M:%S') - Make database starts..."
if [ ! -d $outdir ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: the output directory does not exist."
    exit 1
fi
db_dir=$outdir/"bacnex_db"
mkdir -p $db_dir

# ===== kneaddata =====
mkdir -p $db_dir/kneaddata
kneaddata_database --download human_genome bowtie2 $db_dir/kneaddata

# ===== humann3 =====
mkdir -p $db_dir/humann
humann_databases --download chocophlan full $db_dir/humann
humann_databases --download uniref uniref90_diamond $db_dir/humann
humann_databases --download utility_mapping full $db_dir/humann

# ===== metaphlan =====
mkdir -p $db_dir/metaphlan
metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202403 --bowtie2db $db_dir/metaphlan # latest version

echo "$(date '+%Y-%m-%d %H:%M:%S') - All done."
