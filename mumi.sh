#!/bin/bash 
#
# Compute the MUMi similarity value between two given complete genome sequences. If a genome sequence is
# contained within a directory, all chromosomes of the genome sequences are merged before the genomes are
# compared with each other.
#
# INPUT
#   seq1  GenBank file or directory containing GenBank files for the same genome
#   seq2  GenBank file or directory containing GenBank files for the same genome
#   -p    optional prefix used for creation of temporary file names (default: "MUMI")
#   -t    optional directory for storage of temporary files (default: "/tmp")
#
# syntax: mumi [-p prefix] [-t tmp_dir] seq1 seq2

# process command line options
prefix="MUMi"
tmp_dir="/tmp"
while getopts 'p:t:' option
do
  case ${option} in
    p) prefix=${OPTARG};;
    t) tmp_dir=`echo "${OPTARG}" | sed -e 's/\/*$//'`;;
    ?) echo "Usage: mumi [-p prefix] [-t tmp_dir] seq1 seq2" >&2
       exit 1;;
  esac
done
let "numoptions = ${OPTIND}-1"
shift ${numoptions}

# process command line arguments
if [ $# -ne 2 ]
then
  echo "$0: invalid number of arguments" >&2
  echo "Usage: mumi [-p prefix] [-t tmp_dir] seq1 seq2" >&2
  exit 1
fi
if [ ! -f $1 -a ! -d $1 ]
then
  echo "$0: illigal argument" >&2
  echo "Usage: mumi [-p prefix] [-t tmp_dir] seq1 seq2" >&2
  exit 1
fi
if [ ! -f $2 -a ! -d $2 ]
then
  echo "$0: illigal argument" >&2
  echo "Usage: mumi [-p prefix] [-t tmp_dir] seq1 seq2" >&2
  exit 1
fi
seq1_file=`echo "$1" | sed -e 's/\/*$//'`                         # remove final slash from directory name
seq2_file=`echo "$2" | sed -e 's/\/*$//'`
seq1_name=`echo "${seq1_file}" | awk -F'/' '{print $NF}'`         # extract name as last part of file or directory name
seq2_name=`echo "${seq2_file}" | awk -F'/' '{print $NF}'`
seq1_fasta="${tmp_dir}/${prefix}_${seq1_name}.fasta"              # construct temporary FASTA file names
seq2_fasta="${tmp_dir}/${prefix}_${seq2_name}.fasta"
mumfile="${tmp_dir}/${prefix}_${seq1_name}_${seq2_name}.mummer"   # construct temporary file name to output results

# convert GenBank files to (concatenated) files in FASTA format
# echo "converting GenBank files to (concatenated) FASTA files ..."
echo ">${seq1_name}" > ${seq1_fasta}
if [ -d ${seq1_file} ]
then
  for seqfile in `grep -H '^DEFINITION' ${seq1_file}/*.gbk | grep -v 'plasmid' | sort | cut -d':' -f1`
  do
    seqret -sequence ${seqfile} -sformat gb -osf fasta -stdout -auto 2> /dev/null | tail -n +2 >> ${seq1_fasta}
    while [ $? -ne 0 ]
    do
      seqret -sequence ${seqfile} -sformat gb -osf fasta -stdout -auto 2> /dev/null | tail -n +2 >> ${seq1_fasta}
    done
  done
else
  seqret -sequence ${seq1_file} -sformat gb -osf fasta -stdout -auto 2> /dev/null | tail -n +2 >> ${seq1_fasta}
  while [ $? -ne 0 ]
  do
    seqret -sequence ${seq1_file} -sformat gb -osf fasta -stdout -auto 2> /dev/null | tail -n +2 >> ${seq1_fasta}
  done
fi

echo ">${seq2_name}" > ${seq2_fasta}
if [ -d ${seq2_file} ]
then
  for seqfile in `grep -H '^DEFINITION' ${seq2_file}/*.gbk | grep -v 'plasmid' | sort | cut -d':' -f1`
  do
    seqret -sequence ${seqfile} -sformat gb -osf fasta -stdout -auto 2> /dev/null | tail -n +2 >> ${seq2_fasta}
    while [ $? -ne 0 ]
    do
      seqret -sequence ${seqfile} -sformat gb -osf fasta -stdout -auto 2> /dev/null | tail -n +2 >> ${seq2_fasta}
    done
  done
else
  seqret -sequence ${seq2_file} -sformat gb -osf fasta -stdout -auto 2> /dev/null | tail -n +2 >> ${seq2_fasta}
  while [ $? -ne 0 ]
  do
    seqret -sequence ${seq2_file} -sformat gb -osf fasta -stdout -auto 2> /dev/null | tail -n +2 >> ${seq2_fasta}
  done
fi

# process sequences by mummer
mummer -mum -b -c -l 19 ${seq1_fasta} ${seq2_fasta} > ${mumfile} 2> /dev/null

# get sequence length
seq1_len=`tail -n +2 ${seq1_fasta} | tr -d '\n\r ' | wc -c`
seq2_len=`tail -n +2 ${seq2_fasta} | tr -d '\n\r ' | wc -c`

# process mummer output
# echo "processing mummer output ..."
awk -v seq1_len=${seq1_len} -v seq2_len=${seq2_len} -v seq1_name=${seq1_name} -v seq2_name=${seq2_name} '
# forward or reverse hit for second sequence
/^>/ { if ($0 ~ /Reverse/) reverse=1; next }

# mark positions covered by MUMs
{
  len+=$3
  for(i=$1;i<$1+$3;++i) seq1[i-1]=1
  if (reverse==1)
    for(i=$2-$3+1;i<=$2;++i) seq2[i-1]=1
  else
    for(i=$2;i<$2+$3;++i) seq2[i-1]=1
}

# determine MUM-index
END {
  # compute MUM-coverages of both genomes
  for(i=0;i<seq1_len;++i) seq1_cov+=seq1[i]
  for(i=0;i<seq2_len;++i) seq2_cov+=seq2[i]

  # compute different versions of MUMi similarity value
  sim1=seq1_cov/seq1_len
  sim2=seq2_cov/seq2_len
  sim3=(seq1_cov + seq2_cov)/(seq1_len + seq2_len)
  sim4=0.5*(sim1 + sim2)

  # output results
  printf("%s\t%s\t%d\t%d\t%d\t%d\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n",seq1_name,seq2_name,seq1_len,seq2_len,seq1_cov,seq2_cov,sim1,sim2,sim3,sim4)
}
' ${mumfile}

# remove temporary files
rm -f ${seq1_fasta}
rm -f ${seq2_fasta}
rm -f ${mumfile}
