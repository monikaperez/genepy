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


function launchandco(){
  gcloud compute instances start broadproject2 --zone us-east1-b &&\
  gcloud compute instances list &&\
  gcloud compute config-ssh &&\
  gcpconnect broadproject2.us-east1-b.jkobject
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
  echo $1,$2,$3,$4,$5,$6
  cat $1 | parallel -j $5 "/usr/local/opt/curl-openssl/bin/curl -L "$6" {} | gsutil cp - "$2$3'$(echo {} | cut -d \/ -f '$4')'
}


function gitpush(){ 
  git add . && git commit -m $1 && git push
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


####### AWESOME COMMANDS

function extract {
 if [ -z "$1" ]; then
    # display usage if no parameters given
    echo "Usage: extract <path/file_name>.<zip|rar|bz2|gz|tar|tbz2|tgz|Z|7z|xz|ex|tar.bz2|tar.gz|tar.xz>"
    echo "       extract <path/file_name_1.ext> [path/file_name_2.ext] [path/file_name_3.ext]"
    return 1
 else
    for n in $@
    do
      if [ -f "$n" ] ; then
          case "${n%,}" in
            *.tar.bz2|*.tar.gz|*.tar.xz|*.tbz2|*.tgz|*.txz|*.tar) 
                         tar xvf "$n"       ;;
            *.lzma)      unlzma ./"$n"      ;;
            *.bz2)       bunzip2 ./"$n"     ;;
            *.rar)       unrar x -ad ./"$n" ;;
            *.gz)        gunzip ./"$n"      ;;
            *.zip)       unzip ./"$n"       ;;
            *.z)         uncompress ./"$n"  ;;
            *.7z|*.arj|*.cab|*.chm|*.deb|*.dmg|*.iso|*.lzh|*.msi|*.rpm|*.udf|*.wim|*.xar)
                         7z x ./"$n"        ;;
            *.xz)        unxz ./"$n"        ;;
            *.exe)       cabextract ./"$n"  ;;
            *)
                         echo "extract: '$n' - unknown archive method"
                         return 1
                         ;;
          esac
      else
          echo "'$n' - file does not exist"
          return 1
      fi
    done
fi
}

myinfo () {
  printf "CPU: "
  cat /proc/cpuinfo | grep "model name" | head -1 | awk '{ for (i = 4; i <= NF; i++) printf "%s ", $i }'
  printf "\n"

  cat /etc/issue | awk '{ printf "OS: %s %s %s %s | " , $1 , $2 , $3 , $4 }'
  uname -a | awk '{ printf "Kernel: %s " , $3 }'
  uname -m | awk '{ printf "%s | " , $1 }'
  kded4 --version | grep "KDE Development Platform" | awk '{ printf "KDE: %s", $4 }'
  printf "\n"
  uptime | awk '{ printf "Uptime: %s %s %s", $3, $4, $5 }' | sed 's/,//g'
  printf "\n"
  cputemp | head -1 | awk '{ printf "%s %s %s\n", $1, $2, $3 }'
  cputemp | tail -1 | awk '{ printf "%s %s %s\n", $1, $2, $3 }'
  #cputemp | awk '{ printf "%s %s", $1 $2 }'
}

kp () {
  ps aux | grep $1 > /dev/null
  mypid=$(pidof $1)
  if [ "$mypid" != "" ]; then
    kill -9 $(pidof $1)
    if [[ "$?" == "0" ]]; then
      echo "PID $mypid ($1) killed."
    fi
  else
    echo "None killed."
  fi
  return;
}

ssd () {
  echo "Device         Total  Used  Free  Pct MntPoint"
  df -h | grep "/dev/sd"
  df -h | grep "/mnt/"
}

