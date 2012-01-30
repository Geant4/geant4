#
#  Run the test program with different input files
#   Check output for killed tracks
#
set exe=$G4BIN/$G4SYSTEM/testNTST
set runid=$1; shift
set n=$1; shift

echo "$0 : Script parameter: runid= $runid,  runno=$n"

echo "Machine name, load and users"
hostname
w | grep -v days
# uptime

foreach inp ( 2xa 2xb 2xc 2a 2b 2c )
   uptime
   set suffix=run${inp}-${runid}.n$n     # g460c0dev.pcg2.DetCpr.ChF.n$n
   echo "Running with input ${inp} "

   time $exe run$inp.mac > & oer.$suffix
   #    ****
   tail -18 oer.$suffix | egrep -v 'delet|[Rr]eport|Proce|trajectories|getFieldStats|/vis' | \
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
w | grep -v days
exit

#  gprof $exe > gprof.$suffix
#  mv gmon.out gmon.$suffix
#  gzz gmon.$suffix
