#!/bin/bash
# trial prints the input
function print_my_input() {
  echo 'Your input: ' $1
}

# download from youtube 
# based on a playlist ID (see in the url list\=ID) $1 
# and the start number of the playlist $2 
function update_my_playlists() {
  echo 'Code it first you fool '
  mkdir $1;
  cd $1;
  start=${2:-0}
  youtube-dl --no-overwrite -i --playlist-start $start --write-thumbnail --yes-playlist --output \
  '%(title)s.%(ext)s' --restrict-filenames --extract-audio --audio-quality 0 \
  --audio-format "mp3" --no-check-certificate \
  --user-agent "default" "https://www.youtube.com/playlist\?list\="$1 &&\
  for i in *.jpg; do eyeD3 --add-image $i:FRONT_COVER ${i%.jpg}.mp3; done &&\
  rm *.jpg
}


cpp-run() {
    echo "Compiling file..."
    g++ -o "$1" "$1.cpp"
    echo "Compiled! Enter input :D"
    ./"$1"
}
# cpp-run filename

c-run() {
    echo "Compiling file..."
    gcc -o "$1" "$1.c"
    echo "Compiled! Enter input :D"
    ./"$1"
}

# convert all the set of images in the current folder
# into a gif names $1
function convert2gif() {
  echo 'Converting';
  mkdir frames;
  ffmpeg -i $1 -r 24 'frames/frame-%03d.jpg';
  cd frames;
  convert -delay 24 -loop 0 *.jpg $1.gif;
  mv $1.gif ..;
  cd ..;
  rm -rf frames;
}


function launchandco(){
  # will launch and connect with jkconnect to a gcp instance 
  instance=${1:-broadproject2}
  zone=${2:-us-east1-b}
  name=${3:-jkobject}
  gcloud compute instances start $instance --zone $zone &&\
  gcloud compute instances list &&\
  gcloud compute config-ssh &&\
  jkconnect $instance'.'$zone'.'$name
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

function removeFirstLineOf(){
  echo "$(tail -n +2 $file)" > $file
}

function create_dx_urls_from(){
  # given a dx folder url $1, will create a list of downloadable file links from the files $2 within it
  dx ls $1$2 >list.txt
  # getting the urls
  for i in $(cat list.txt); do 
    dx make_download_url --duration "1w" "$1$i" >> urls.txt; 
  done;
}

function jkconnect(){
  #connects to the server by also connecting rsublime and more
  #fswatch -o . | while read f; do rsync -zrq --max-size=200m  ./ $1:current; done & \
  ssh -L 8157:127.0.0.1:8890 \
  -R 52698:localhost:52698 $1 
}

function upload_urls_to_gcp(){ 
  # args:
  #   1 filepath with list of server filepath to curl
  #   2 gs://bucket/
  #   3 bucket/filepaht/
  #   4 position of the name in the serverfilepath (1/2/3/4)
  #   5 threads number
  #   6 additional curl parameters
  echo "if you don't have paralllel run the following: \
  (wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - http://pi.dk/3) | bash \
  you also need to have a brewed curl in /usr/local/opt/curl-openssl/bin/curl "
  echo $1,$2,$3,$4,$5,$6
  cat $1 | parallel -j $5 "/usr/local/opt/curl-openssl/bin/curl -L "$6" {} | gsutil cp - "$2$3'$(echo {} | cut -d \/ -f '$4')'
}


function gitpush(){ 
  # add.commitpush with 1 being commit message
  git add . && git commit -m $1 && git push
}

function get_unuploadedfiles_from_dx_to_gcp(){
  # gets the flies in fx not in gcp from a list of filepaths in dx
  # Args:
  #  dxfolder url to dxfolder
  #  gcpfolder names
  #  urls filepath of filepath urls
  # WIP ###############################
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
  # will continue a stopped gcp upload with the compose command from a set of curl files
  # Args:
  #   1 filepath to list of filepaths in gcp
  #   2 position in the filepaths where the filename is 1/2/3/4
  #   3 gcp bucket path to download into
  IFS=$'\n'
  cut=${2:-7}
  for i in $(cat $1); do
    file=$(echo ${i} | awk '{ print $1}');
    size=$(echo ${i} | awk '{ print $2}');
    echo "downloading $file of size: $size";
    name=$(echo $file | cut -d \/ -f $cut);
    curl -C $size -L -s $file | gsutil cp - $3'retry_'$name;
    gsutil compose $3$name $3"retry_"$name $3$name
    gsutil rm $3"retry_"$name
  done;
}


function rename_bunch(){for file in $1; do mv "$file" "${file/$2/}"; done}
# rename all files in $1 with the addition $2

####### AWESOME COMMANDS

function extract {
 # extract almlost anything
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
  # displays all my info
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
  # kill the process with name $1
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

ssd () { # gets my ssd storage left (only in linux)
  echo "Device         Total  Used  Free  Pct MntPoint"
  df -h | grep "/dev/sd"
  df -h | grep "/mnt/"
}

remove_port_rsub () {
  username=${1:-jeremie}
  ps -u username
  # kill # of first process named sshd 
}  
