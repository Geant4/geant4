#!/bin/sh
identify_release() {
   local libc=`basename /lib/libc-*.so .so`
   local release=slc5
   case $libc in
       libc-2.3*)
           release=slc4
           ;;
       libc-2.5*)
           ;;
       *)
           echo "Unable to identify correct binary "\
"package for this system (libc version ${libc}): using $release" 1>&2
   esac
   echo "$release"
}
