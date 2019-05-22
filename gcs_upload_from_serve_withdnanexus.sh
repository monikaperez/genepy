#!bin/bash 
# listing the files
dx ls "immediate/bam/*.WholeGenome*" >list.txt
# getting the urls
for i in $(cat list.txt):
do dx make_download_url --duration "1w" "immediate/bam/$i" >> urls.txt; done
# for 8 of them at a time, downloading them directly to the bucket in parrallel 
nohup cat urls.txt | xargs -P 8 -n 1 wget -qO- | gsutil cp - 'gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/' > nohup.out 2>&1 &
