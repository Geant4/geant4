# runlimit.sh is sourced from a sh script 
# and calls ulimit as the bash/sh implementation.
# (tcsh etc call limit) and is builtin to the shell (so man bash).
#
# prevent a coredump, 
# signal looping tests to stop
# prevent memory leaks from crashing pcgeant4

#ulimit -t 600
# after cut_per_mat...    
ulimit -t 6000
# check without!
#ulimit -v 500000
# may be it's usefull to restrict output file? Don't think so...
# even dengerous: lost output (only 5 Mb :)
ulimit -f 10000
# BTW for core preventing 
ulimit -c 0 # by default (?)

echo "STT:RunLimit  $G4SYSTEM  $G4LARGE_N $TMPDIR"
ulimit -a

# export TMPDIR=/afs/cern.ch/user/s/stesting/stt/tmpdir
