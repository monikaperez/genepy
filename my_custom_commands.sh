#!/bin/bash
# trial prints the input
function print_my_input() {
  echo 'Your input: ' $1
}

# download from youtube
function update_my_playlists() {
  echo 'Code it first you fool '
}

# convert to gif
function convert2gif() {
  echo 'Converting';
  mkdir frames;
  ffmpeg -i $1 -r 5 'frames/frame-%03d.jpg';
  cd frames;
  convert -delay 20 -loop 0 *.jpg $1.gif;
  mv $1.gif ..;
  cd ..;
  rm frames;
}

function retrieve_remote_file_size(){
  $ch = curl_init($url);

  curl_setopt($ch, CURLOPT_RETURNTRANSFER, TRUE);
  curl_setopt($ch, CURLOPT_HEADER, TRUE);
  curl_setopt($ch, CURLOPT_NOBODY, TRUE);

  $data = curl_exec($ch);
  $size = curl_getinfo($ch, CURLINFO_CONTENT_LENGTH_DOWNLOAD);

  curl_close($ch);
  return $size;
}

function create_dx_urls_from(){
  dx ls $path$file >list.txt
  # getting the urls
  for i in $(cat list.txt); do 
    dx make_download_url --duration "1w" "$path$i" >> urls.txt; 
  done;
}

function gcpconnect(){
  #fswatch -o . | while read f; do rsync -zrq --max-size=200m  ./ $1:current; done & \
  ssh -L 8157:127.0.0.1:8890 \-R 52698:localhost:52698 $1 
}

function upload_urls_to_gcp(){ 
  echo "if you don't have paralllel run the following: \
  (wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - http://pi.dk/3) | bash"
  cat urls.txt | parallel -j 10 'curl -L -s {} | gsutil cp - \
  $bucket$folder$(echo {} | cut -d \/ -f 7)'
}


function get_unuploadedfiles_from_dx_to_gcp(){
  for i in $(cat $urls); do
    echo $i;
    sizeA=$(dx describe --json $dxfolder$(echo $i | cut -d \/ -f 7) | jq '.size');
    sizeB=$(gsutil du $bucket$gcpfolder$(echo $i | cut -d \/ -f 7) | cut -d \  -f 1);
    if (($sizeA!=$sizeB)); then
      echo "$i $sizeB" > filters.txt;
    fi;
    echo "$sizeB, $sizeA";
  done;
}


function continue_upload(){
  IFS=$'\n'
  for i in $(cat $file); do
    file=$(echo ${i} | awk '{ print $1}');
    size=$(echo ${i} | awk '{ print $2 }');
    echo "downloading $file of size: $size";
    name=$(echo $file | cut -d \/ -f 7);
    curl -C $size -L -s $file | gsutil cp - "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/retry_"$name;
    gsutil compose "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/"$name "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/retry_"$name "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/"$name
    gsutil rm "gs://fc-secure-fd4c24cf-6bfd-410a-9bca-e02642da12f8/immediate/bam_wg/retry_"$name
  done;
}

function rename_bunch(){for file in $folder; do mv "$file" "${file/mapped./}"; done}


