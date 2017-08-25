#!/bin/bash

nframes=`grep 'frame' WATER.MOLINDEX | wc -l`
echo $nframes

csplit -s -n 5 WATER.MOLINDEX '/frame/' '{*}'
mkdir split_frames
mv xx* split_frames/

#for i in {00001..00005}; do
for i in {00001..25000}; do
  sort -n -k 1 split_frames/xx$i > sort_frame.$i

  tail -n +2 sort_frame.$i | awk '{print $1}' > tail.$i
  sort -n ../1-512.dat tail.$i | uniq -u | awk '{print $1, 1}' > 0dh_tmp.$i

  uniq -u -w 12 sort_frame.$i > 1DH_tmp.$i
  cat 1DH_tmp.$i 0dh_tmp.$i > 01dh_tmp.$i
  sort -n -k 1 01dh_tmp.$i > 01DH_tmp.$i

  uniq -D -w 12 sort_frame.$i > 2dh_tmp.$i
  head -1 sort_frame.$i > 2DH_tmp.$i
  cat 2dh_tmp.$i >> 2DH_tmp.$i

  odh=`wc -l 01DH_tmp.${i} | awk '{print $1-1}'`
  tdh=`wc -l 2DH_tmp.${i} | awk '{print $1-1}'`

#  echo "$odh" "$tdh"
  awk -v var=$odh 'NR==1 {$3=var};1' 01DH_tmp.${i} > 01DH_tmp && mv 01DH_tmp 01DH_tmp.${i}
  awk -v var=$tdh 'NR==1 {$3=var};1' 2DH_tmp.${i} > 2DH_tmp && mv 2DH_tmp 2DH_tmp.${i}

#  cat sort_frame.$i >> SORT_WATER.MOLINDEX
  cat 01DH_tmp.${i} >> 01DH_SORT_WATER.MOLINDEX
  cat 2DH_tmp.${i} >> 2DH_SORT_WATER.MOLINDEX 
done

