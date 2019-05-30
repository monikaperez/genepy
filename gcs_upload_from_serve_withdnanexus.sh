#!bin/bash 
# listing the files
echo "be carefull, in this DB two files can have the same name!! you won't see it until it is too late"
dx ls "immediate/bam/*.WholeGenome*" >list.txt
# getting the urls
for i in $(cat list.txt); do 
	dx make_download_url --duration "1w" "immediate/bam/$i" >> urls.txt; 
done;
# for 8 of them at a time, downloading them directly to the bucket in parrallel 


for i in $(cat urls.txt);  do 
	curl -L -s $i | gsutil cp - "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/"$(echo $i | cut -d \/ -f 7) &  
done;


for i in $(cat urls.txt); do
	echo $i;
	sizeA=$(dx describe --json "immediate/bam/"$(echo $i | cut -d \/ -f 7) | jq '.size');
	sizeB=$(gsutil du "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/"$(echo $i | cut -d \/ -f 7) | cut -d \  -f 1);
	if (($sizeA!=$sizeB)); then
		echo "$i $sizeB" > filters.txt;
	fi;
	echo "$sizeB, $sizeA";
done;

IFS=$'\n'
for i in $(cat filters.txt); do
	file=$(echo ${i} | awk '{ print $1}');
	size=$(echo ${i} | awk '{ print $2 }');
	echo "downloading $file of size: $size";
	name=$(echo $file | cut -d \/ -f 7);
	curl -C $size -L -s $file | gsutil cp - "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/retry_"$name;
	gsutil compose "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/"$name "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/retry_"$name "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/"$name
	gsutil rm "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/retry_"$name
done;




## VERSION WITHOUT BREAKING LINES
for i in $(cat urls.txt); do echo $i; sizeA=$(dx describe --json "immediate/bam/"$(echo $i | cut -d \/ -f 7) | jq '.size'); sizeB=$(gsutil du "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/"$(echo $i | cut -d \/ -f 7) | cut -d \  -f 1); if (($sizeA!=$sizeB)); then echo "$i, $sizeB" > filters.txt; fi; echo "$sizeB, $sizeA"; done;