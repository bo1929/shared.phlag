#!/bin/bash
nw_topology $1 | nw_distance -n -m m - | awk 'NR>1 {
    row=$1;
    for (col=2; col<=NF; col++) {
      if ($col == 3) {
        print row, header[col-1];
      }
    }
  } NR==1 { for (i=1; i<=NF; i++) header[i]=$i; }' \
    | awk '$1>$2 {print $0}' \
    | python sort_triplets_ultrametricity.py $1 -
