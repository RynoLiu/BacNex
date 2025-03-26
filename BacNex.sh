#!/usr/bin/env bash

VERSOIN="beta 2.1"

main_usage(){
echo "
Usage: BacNex.sh [modules]

Modules:

    preprocess               Preprocess Module the input data.
    make_table               Generate BacNex standard input files.
    make_database            Download BacNex preprocess database.
    app                      Run BacNex pipeline application.

    -h,  --help              Display this help message.
    -v,  --version           Display the version of BacNex.
"
}

SCRIPT_DIR=$(dirname "$(realpath "$0")")

if [ "$1" = "preprocess" ]; then
	time ${SCRIPT_DIR}/modules/preprocess.sh ${@:2}
elif [ "$1" = "make_table" ]; then
    time ${SCRIPT_DIR}/modules/make_table.py ${@:2}
elif [ "$1" = "app" ]; then
    time ${SCRIPT_DIR}/modules/app.sh ${@:2}
elif [ "$1" = "make_database" ]; then
    time ${SCRIPT_DIR}/modules/make_database.sh ${@:2}
elif [ "$1" = "-v" ] || [ "$1" = "--version" ]; then
    echo "BacNex version: ${VERSOIN}"
elif [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	main_usage
fi
