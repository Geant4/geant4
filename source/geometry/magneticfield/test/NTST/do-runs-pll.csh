#
#  Run the test program with different input files
#   Check output for killed tracks
#
set runid=$1; shift
set n=$1; shift

echo "$0 : Script parameter: runid= $runid,  runno=$n"

# Check the environment
if( $?G4BIN ) then
   if( -d $G4BIN ) then 
      set exe=$G4BIN/$G4SYSTEM/testNTST
   else 
      echo "Directory G4BIN=$G4BIN does not exist. "
   endif 
else
   ## Report issue
   if( ! $?G4BIN ) then    
      echo "Environment variable G4BIN is not set."
   endif

   ## Try the second
   if( $?G4WORKDIR ) then 
      if( -d $G4WORKDIR/bin ) then 
	 set exe=$G4WORKDIR/bin/$G4SYSTEM/testNTST
      else
         echo "Directory G4WORKDIR/bin ( $G4WORKDIR/bin ) does not exist. "
         echo "Directory for executable is not known, or cannot be found. "
         echo "Fatal ERROR.  Exiting."
         exit 1
      endif
   else
      ## Report issue
      echo "Environment variable G4WORKDIR is not set."
      exit 2
   endif
endif
if( ! -x $exe ) then
   echo "Fatal ERROR> Cannot locate executable $exe "
   exit 3
endif

echo "Machine name, load and users"
hostname
w | grep -v days
# uptime

foreach inptype ( 2x 2 )
  uptime
  foreach var ( a b c )
    set inp=${inptype}${var}
    set suffix=run${inp}-${runid}.n$n     # g460c0dev.pcg2.DetCpr.ChF.n$n
    echo "Running with input ${inp} "
 
    time $exe run$inp.mac > & oer.$suffix   &
  end
  echo "Waiting for the results of all runs labelled >" $inptype  "< ....."
  wait
  echo "Jobs done. Analysing."
  forech var ( a b c )
    set inp=${inptype}${var}
    set suffix=run${inp}-${runid}.n$n     # g460c0dev.pcg2.DetCpr.ChF.n$n
    tail -18 oer.$suffix | \
                      egrep -v 'delet|[Rr]eport|Proce|trajectories|getFieldStats|/vis' | \
                      grep -v 'Number of Vertices =   23.1 +-  0.515' | \
 		      grep -v 'User='  | \
 		      sort -u
    set kildtrax=`grep -c  killing oer.$suffix`
    echo "Killed $kildtrax number of tracks"

    grep 'track has' oer.$suffix | \
      awk ' /MeV/ { en += $4; } \!/MeV/ { print; } END { print "Total energy killed= ", en, " MeV"; } '
 
    if( `grep -c G4Prop oer.$suffix` ) then
       echo "Iterations with many trials: "
       ./count.maxIt-PiF  oer.$suffix
    else
       echo "Iterations of trials: no statistics collected"
    endif

    echo " "   ## Separate the runs
  end
end
w | grep -v days
exit

#  gprof $exe > gprof.$suffix
#  mv gmon.out gmon.$suffix
#  gzz gmon.$suffix
