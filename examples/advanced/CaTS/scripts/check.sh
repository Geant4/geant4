#!/bin/bash -l
msg="=== $BASH_SOURCE :"

info(){ cat << EOI
$msg checking the environment 
   OPTICKS_HOME   : $OPTICKS_HOME
   CMAKE_PREFIX_PATH : 
$(echo $CMAKE_PREFIX_PATH | tr ":" "\n")
The CMAKE_PREFIX_PATH is expected to contain about seven prefix directories including::
   OPTICKS_PREFIX            : $OPTICKS_PREFIX
   OPTICKS_PREFIX/externals  : $OPTICKS_PREFIX/externals
Note that environment setup must be done in login scripts : .bashrc .bash_profile .opticks_config etc.. 
Doing environment setup just in the current session will not work as scripts often invoke login scripts. 
EOI
}

rc=0
if [ -z "$OPTICKS_HOME" -o -z "$OPTICKS_PREFIX" -o -z "$CMAKE_PREFIX_PATH" ]; then 
   echo $msg missing required envvars : your need to source .opticks_config 
   rc=1 
fi 
if [ "$CMAKE_PREFIX_PATH" == ${CMAKE_PREFIX_PATH/$OPTICKS_PREFIX:/} ]; then 
   echo $msg CMAKE_PREFIX_PATH does not contain the expected prefix : $OPTICKS_PREFIX : you need to invoke opticks-setup from .opticks_config
   rc=2
fi
if [ "$CMAKE_PREFIX_PATH" == ${CMAKE_PREFIX_PATH/$OPTICKS_PREFIX\/externals\:/} ]; then 
   echo $msg CMAKE_PREFIX_PATH does not contain the expected prefix : $OPTICKS_PREFIX/externals : you need to invoke opticks-setup from .opticks_config
   rc=3
fi 

if [ $rc -ne 0 ]; then
   info
   echo $msg environment check FAILED : rc $rc 
else
   [ -n "$VERBOSE" ] && info 
   echo $msg environment check PASSED : rc $rc
fi 
exit $rc
