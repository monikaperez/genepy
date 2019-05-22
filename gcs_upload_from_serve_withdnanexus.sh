#!bin/bash 
# listing the files
echo "be carefull, in this DB two files can have the same name!! you won't see it until it is too late"
dx ls "immediate/bam/*.WholeGenome*" >list.txt
# getting the urls
for i in $(cat list.txt);
do dx make_download_url --duration "1w" "immediate/bam/$i" >> urls.txt; done
# for 8 of them at a time, downloading them directly to the bucket in parrallel 
for i in $(cat urls.txt); 
do if ((counter%8 != 0)); 
	then curl -L -s $i | gsutil cp - "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/"$(echo $i | cut -d \/ -f 7) & 
	else 
		curl -L -s $i | gsutil cp - "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/"$(echo $i | cut -d \/ -f 7); 
	fi; 
	((counter++)); 
done;
