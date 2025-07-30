#!/bin/bash

NPBOS_DIR=$1
MASS=$2
Bn=$3
eps=$4
kappa=$5

cd $NPBOS_DIR
bash 2outputs.sh $MASS $Bn $eps $kappa

sleep 0.5

# ファイルを１行ずつ読み込む
FILE="./out1.dat"
if [ ! -f "$FILE" ]; then
    echo "File not found!"
    exit 0
fi
COL=$(wc -l < "$FILE")

while read LINE ; do
    COL=$(($COL - 1))
    if [ ${COL} -lt 9 ] ; then
        LINE_SUB=$(echo "$LINE" | cut -c7-16)
        if [ "$LINE_SUB" = "2  +  ( 1)" ] ; then
            f_2=$(echo "$LINE" | cut -c23-28)
        elif [ "$LINE_SUB" = "4  +  ( 1)" ] ; then
            f_4=$(echo "$LINE" | cut -c23-28)
        elif [ "$LINE_SUB" = "6  +  ( 1)" ] ; then
            f_6=$(echo "$LINE" | cut -c23-28)
        elif [ "$LINE_SUB" = "0  +  ( 2)" ] ; then
            s_0=$(echo "$LINE" | cut -c23-28)
        fi
    fi
done < "$FILE"

echo "${f_2} ${f_4} ${f_6} ${s_0}"