#!/bin/bash
# Run each of the benchmarks and calculate an average

# Create a new tempfile for a copy of the output
tmpfile=$(mktemp /tmp/sombrero_XXXXXXXXX)

# Clean finish
function finish {
  rm -rf "$tmpfile"
}
trap finish EXIT

# Print usage and quit
function usage {
    echo "./sombrero.sh { -n <num-cores> | -H <hostfile> }  [ -w ] [ -s small | medium | large | very_large ]" >&2
    exit 1
}

# safe run
trap "exit 1" TERM
export TOP_PID=$$
function safe_run {
    $@
    if [ $? -ne 0 ]; then
        kill -s TERM $TOP_PID
    fi
}

#### main
# check arguments 
shopt -s extglob
while getopts 'n:H:s:whl:p:' opt ; do
    case "$opt" in
        n)
            case "$OPTARG" in
                (+([0-9])) nodes="-n $OPTARG " ;;
                "") usage ;;
            esac ;;
        H)
            case "$OPTARG" in
                "") usage ;;
                *) nodes="-hostfile $OPTARG " ;;
            esac ;;
        s)
            case "$OPTARG" in
                "small" | "medium" | "large" | "very_large") size="-s $OPTARG " ;;
                *) usage ;;
            esac ;;
        l)
            case "$OPTARG" in
                "") usage ;;
                *) size="-l $OPTARG " ;;
            esac ;;
        p)
            case "$OPTARG" in
                "") usage ;;
                *) partition="-p $OPTARG " ;;
            esac ;;
        w) weak="-w " ;;
        h|?) usage ;;
    esac
done

# From https://stackoverflow.com/a/246128/3113564
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# Run each of the benchmarks, copy output to stdout and the tempfile
touch $tmpfile
for i in {1..6}
do 
    safe_run mpirun $nodes \
    $SCRIPT_DIR/sombrero/sombrero$i $weak $size $partition -v result |\
    tee -a $tmpfile
done

# Calculate the sums of the timings and the flop counts
time=$(grep RESULT $tmpfile | grep -v "Gflops/seconds" | awk '{s+=$7} END {print s}') 
Gflops=$(grep RESULT $tmpfile | grep -v "Gflops/seconds" | awk '{s+=$4} END {print s}') 
echo [RESULT] SUM $Gflops Gflops in $time seconds &&
echo [RESULT] SUM $(echo $Gflops/$time | bc ) Gflops/seconds



