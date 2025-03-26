#!/usr/bin/env bash

function app_help(){
    echo "
Usage: BacNex.sh ko_search [options]

Parameters:
    -w,    --work-dir        The folder from make_table module.
    -t,    --threads         The number of threads to use.
    -b,    --browser         The browser to use. (google-chrome, firefox)
    -h,    --help            Display the help message.
"
}

if [[ ${#} == 0  ]]; then
    app_help
fi

# default parameters
browser="firefox"
working_dir=$(pwd)

# parameters
while [ ${#} -gt 0 ]; do
    error_message="Error: a value is needed for '$1'";
    case $1 in
        -w | --work-dir ) working_dir=${2:?$error_message}; shift 2;;
        -t | --threads ) threads=${2:?$error_message}; shift 2;;
        -b | --browser ) browser=${2:?$error_message}; shift 2;;
        -h | --help) app_help; exit 0;;
        * ) echo "Unknown parameter passed: $1"; app_help; exit 1 ;;
    esac
done


# ===== app main =====
script_dir=$(dirname "$(realpath "$0")")
# env var
export WORKING_DIR="${working_dir}"
export SCRIPT_DIR="${script_dir}"
export THREADS="${threads}"
app_dir=$(realpath "${script_dir}/app.R")

# browser
open_browser=""
if [ $browser == "google-chrome" ]; then
    open_browser="/usr/bin/google-chrome"
elif [ $browser == "firefox" ]; then
    open_browser="/usr/bin/firefox"
else
    echo "Error: browser not supported"; exit 1
fi

# open app
Rscript -e "shiny::runApp('${app_dir}', launch.browser=function(url) { utils::browseURL(url, browser='$open_browser') })"
