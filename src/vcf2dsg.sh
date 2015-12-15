#!/bin/bash

src=/dev/stdin			# source vcf
dst=/dev/stdout			# target dgz
ske=0

## helper
function help()
{
    echo "tranform the VCF into dosage data."
    echo "Usage: $0 -s <VCF> -d <DGZ> -hpf"

    echo "Options:"

    t="    "
    echo "$t-s source *.vcf or *.vcf.gz file to read from."
    echo "$t-- by default it is the standard input."
    echo
    
    echo "$t-d target to write the dosage dataset."
    echo "$t-- by default the standard output is chosen."
    echo
    
    echo "$t-e skip existing output file."
    echo

    echo "$t-h show this help."
    echo
}

## get parameters
function opts()
{
    while getopts ":hpfs:d:" opt; do
	case $opt in
            h)help;exit 0;;
            p)dry=1;;
            b)skp=1;;
            s)src="$OPTARG";;
            d)dst="$OPTARG";;
            \?)
		echo "Invalid option: -$OPTARG" >&2
		help;exit 1;;
	esac
    done
}

## argument check and report
function args()
{
    echo "SRC=$src"
    echo "DST=$dst"
    echo "SKE=$ske"

    ## check the srouce file
    if [ ! -e "$src" ]; then
	echo "E: the source file $src does not exist." >&2
	exit 1
    fi

    ## check the target file
    if [ -e "$dst" -a $ske -ne 0 ]; then
	echo "W: skip existing target file $dst." >&2
	exit 2
    fi
}

function main()
{
    tmp=$(mktemp -dp ./)
    ## copy the input if it came from STDIN
    if [ $src = /dev/stdin ]; then
	cat $src > $tmp/src
	src=$tmp/src
    fi

    ## extract VCF header
    bcftools view -h $src | head -n-1 | gzip > $tmp/hdr.txt.gz

    ## extract subject IDs
    bcftools query -l $src | gzip > $tmp/sid.txt.gz

    ## query map and GT(s) from source VCF
    fmt='%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n'
    bcftools query -f $fmt $src | gzip > $tmp/qry
    qry=$tmp/qry

    ## cut out genomic map
    zcat $qry | cut -f-5 | gzip > $tmp/map.txt.gz

    ## cut out GTs, tranform them into dosage values
    ptn=""				# sed tranformation pattern
    ptn="$ptn;s:[.|/]\{3\}:3:g"		# missing values -> 4
    ptn="$ptn;s:[0|/]\{3\}:0:g"		# AA -> 0
    ptn="$ptn;s:[1-9|/]\{3\}:2:g"	# Aa|aA -> 1
    ptn="$ptn;s:[^\t]\{3\}:1:g"		# aa -> 2
    zcat $qry | cut -f6- | sed $ptn | gzip > $tmp/mtx.txt.gz

    ## cut out genomic variant IDs
    zcat $tmp/map.txt.gz | cut -f3 | gzip > $tmp/vid.txt.gz

    ## pack everything
    dst_dir=$(dirname $dst)
    cd $tmp
    tar -zcf../$ *.txt.gz "$dst"
    echo "xt: success"
}

opts $@
args
main
