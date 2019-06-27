#!/bin/bash

# arg 1 : st_time : Start time of video slice (in mm:ss format)
# arg 2 : duration : Duration of video slice (in mm:ss format)
# arg 3 : filename : video filename, no quotes
# arg 4 : rootname : Root of filenames for naming extracted frames 
#                    (eg.run01_%04d.jpg)
# arg 5 (optional) : fps : frame frequency of video, default is 30 fps

st_time=$1
duration=$2
filename=$3
rootname=$4
fps=${5:-30}

ffmpeg -ss $st_time -t $duration -i $filename -r $fps $rootname
