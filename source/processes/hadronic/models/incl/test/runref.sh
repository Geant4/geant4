#!/bin/sh

# Program for running INCL/ABLA FORTRAN version from command line. It
# generates cu42.in run configuration file based on command line
# arguments. The command line arguments define target type, projectile
# type/energy, the number of events and output file.

#=====================================================================
writeRunFile() # Prepare setup file cu42.in
#=====================================================================
# Arg_1 = Target mass number
# Arg_2 = Target charge number
# Arg_3 = Bullet type
# Arg_4 = Bullet energy
# Arg_5 = Number of events
# Arg_6 = Log file (ASCII)
# Arg_7 = Output of the calculation (HBOOK)
{
rm -f cu42.in
cat > cu42.in <<EOF
'$6'	* name of output file (a few numbers for identification of the calculation)
979678188 1		* first random, choice evapo (KHS=1,GEM=2)
$5 0  0	* number of events, print_event_i, seeds given 0/1=No/Yes (IF yes, 20 values given after)
$3 $4		* type projectile (1:p, 2:n, 6:d), total kinetic energy (MeV)
$1 $2	   	* A Z of the target
45. 1. -2  8. 0	* Nuclear potential (MeV), t=t0*fact, NOSURF, XFOISA, NPAULSTR
'./data/'	* path to get the .tab et .dat (look at "open" in abla_v3p.f and setup_AB.f)	
'./$7'		* ntuple
0  'C12_avat.hbk'	* sortie ntuple avatars 0/1=No/Yes (be careful, yes for a few events only!)     
EOF
}

#=====================================================================
runProgram() # Run INCL/ABLA program
#=====================================================================
# Arg_1 = program name
{
$1
rm -f cu42.in
}

if [ $# -ne 7 ]; then
	echo "Usage: "
	echo -n $0
	echo " <cascade|full> <A> <Z> <projectile> <energy> <events> <file (no suffix)>"
	echo "<cascade|full>   cascade = INCL only, full = INCL/ABLA"
	echo "<projectile>     1 = proton, 2 = neutron, 3 = pi+, 4 = pi-, 5 = pi0"
	echo "                 6 = H2, 7 = H3, 8 = He3, 9 = He4"
	echo "<energy>         Range 100 - 3000 in MeV"
	echo "<file>           Filename without suffix (generates file.hbook and file.out)"
	exit 127
fi

# Main:
cascade_evaporation_program="./cugnon42_khs_gem"
cascade_only_program="./cugnon42_noevapfis"
evaporation_only_program="./evap"

if [ $1 == "cascade" ]; then
    program=$cascade_only_program
fi
if [ $1 == "full" ]; then
    program=$cascade_evaporation_program
fi
if [ $1 == "abla" ]; then
    program=$evaporation_only_program
fi

massnumber=$2
chargenumber=$3
projectiletype=$4
projectileenergy=$5
events=$6 # 10k for production run
logfile=$7".out"
output=$7".hbook"

# writeRunFile $massnumber $chargenumber $targettype $targetenergy $events $logfile $output
# runProgram $program

writeRunFile $massnumber $chargenumber $projectiletype $projectileenergy $events $logfile $output
runProgram $program

h2root $output
