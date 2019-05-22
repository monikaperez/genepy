#!bin/bash 
# listing the files
echo "be carefull, in this DB two files can have the same name!! you won't see it until it is too late"
dx ls "immediate/bam/*.WholeGenome*" >list.txt
# getting the urls
for i in $(cat list.txt);
do dx make_download_url --duration "1w" "immediate/bam/$i" >> urls.txt; done
# for 8 of them at a time, downloading them directly to the bucket in parrallel 
