#!/bin/sh 

debug=""
while ( test $# -gt 0 ) ; do
 case $1 in 
   -dbg) debug=gdb
         core=$2
         shift
	 ;;
      *)    echo $1
   echo $1 | grep -q -e "-DCMD=" && cmd=`echo $1 | sed -e 's/\"//g'`
   echo $1 | grep -q -e "-DENV=" && test_env=`echo $1 | sed -e 's/\"//g'`
	 
 esac
 shift	 

done
 

#read parts

part=""
for part in $parts; do 
   #echo $part
   echo $part | grep -q -e "-DCMD=" && cmd=`echo $part | sed -e 's/\"//g'`
   echo $part | grep -q -e "-DENV=" && test_env=`echo $part | sed -e 's/\"//g'`
   
done

  #echo "command= $cmd"
  #echo "environment= $test_env"
  
for var in `echo $test_env | sed -e 's/-DENV=//' | sed -e 's/#/ /g'`; do
  val=`echo $var | awk -F@ '{print $2}'`
  var=`echo $var | awk -F@ '{print $1}'`
  export $var=$val 
done

env | grep G4
echo "====>>> cmd=$cmd"
cmdargs=`echo $cmd | sed -e 's/#/ /' | awk '{print $2}' | sed -e 's/#/ /g'`
cmd=`echo $cmd | sed -e's/#/ /' | awk '{print $1}' | sed -e 's/-DCMD=//'`

cmddir=`dirname $cmd`
cmd=`echo $cmd | awk -F/ '{ print $NF }'`
echo dir=$cmddir

if test -d $cmddir ; then 
   cd $cmddir || exit
   echo "Run in directory  $cmddir"
   echo "Command used:     $cmd"
   echo "Arguments:        $cmdargs"
   if test -z "$debug" ; then 
      ./$cmd $cmdargs
   else
      echo starting bash, exit when done.
      /bin/bash   
   fi
else
   echo "Directory for test $cmddir not found"
fi   
