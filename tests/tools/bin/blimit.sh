# blimit.sh is sourced from a sh script 
# and calls ulimit as the bash/sh implementation.
# (tcsh etc call limit) and is builtin to the shell (so man bash).
#
# prevent a coredump, 
# signal looping tests to stop
# prevent memory leaks from crashing pcgeant4

ulimit -t 600     # the longest test is 508 on hpplus.
ulimit -v 500000  # is this the way to stop memory leaks ?
ulimit -f 200000  # largest library/executable is ...

echo "STT:BatchLimit  $G4SYSTEM  $G4LARGE_N $TMPDIR"
ulimit -a

# export TMPDIR=/afs/cern.ch/user/s/stesting/stt/tmpdir
