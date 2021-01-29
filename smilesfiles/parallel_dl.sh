#!/bin/bash

#read curl lines into array 'lines'
#IFS=$'\n' read -d '' -r -a lines < ZINC-downloader-2D-smi.curl
IFS=$'\n' read -d '' -r -a lines < short.curl

#count the number of lines:
#numline=`wc ZINC-downloader-2D-smi.curl | awk '{print $1}' -`
numline=`wc short.curl | awk '{print $1}' -`

#how many jobs to run in parallel:
N=4

for i in $(seq 1 $numline); do
    (
        # .. do your stuff here
        echo "starting task $i.."
	echo "task is ${lines[$i]}"
	`${lines[$i]}`
    ) &

    # allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
        wait -n
    fi

done

# no more jobs to be started but wait for pending jobs
# (all need to be finished)
wait

echo "all done"
