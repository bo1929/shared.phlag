#!/bin/bash
awk -F'[:,]' '
NR==1 { 
    val=$NF; gsub(/^[ \t]+/, "", val) 
} 
NR>1 { 
    for(i=1; i<=NF; i++) a[i, NR-1] = $i; 
    if(NF>max_f) max_f=NF; 
    rows++ 
} 
END { 
    for(i=1; i<=max_f; i++) { 
        for(j=1; j<=rows; j++) { 
            printf "%s,", a[i,j] 
        } 
        print val 
    } 
  }' <(tail -n+3 $1)
