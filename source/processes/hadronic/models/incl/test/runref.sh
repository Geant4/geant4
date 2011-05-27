#!/bin/bash

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
if [ $3 -eq 1 ]; then
    projA=1
    projZ=1
fi
if [ $3 -eq 2 ]; then
    projA=1
    projZ=0
fi
if [ $3 -eq 3 ]; then
    projA=-1
    projZ=1
fi
if [ $3 -eq 4 ]; then
    projA=-1
    projZ=0
fi
if [ $3 -eq 5 ]; then
    projA=-1
    projZ=-1
fi
if [ $3 -eq 6 ]; then
    projA=2
    projZ=1
fi
if [ $3 -eq 7 ]; then
    projA=3
    projZ=1
fi
if [ $3 -eq 8 ]; then
    projA=3
    projZ=2
fi
if [ $3 -eq 9 ]; then
    projA=4
    projZ=2
fi
if [ $3 -eq -12 ]; then
    projA=12
    projZ=6
fi

if [ $1 -le 18 ]; then
    evapo=4
else
    evapo=1
fi

rm -f cu43.in
cat > cu43.in <<EOF
'$6'	* fichier de sortie
'INCL43 special HI for GEANT4'    * Comment
38035 $evapo  	* first random, choice evapo (KHS=1,GEM=2,Dresner=4)
$5  0  0	* Beam projectiles, print_event_i, seeds given 0/1=N/Y
$projA $projZ $4	* Projectile A and Z (A=-1 for pions), TOTAL kinetic energy (MeV)
$1 $2		* A Z of the target
45. 1. -2  10. 0	* Nuclear Pot. (MeV), t=t0*fact, NOSURF, XFOISA, NPAULSTR
0		* Coulomb barrier IN (0=NO,1=LAHET,2=INCL)
0  0  0.        * Not used, clusters (0=NO,1=Coul bar R0,2=coul bar R0+RMS+FNECK),FNECK  
0     		* Forced absorption (0/1=NO/YES)
'./data/'	* Path for .tab et .dat
'$7'		* ntuple physics
0  'pb1000_cu43_abla_loce_t5_avat.hbk'	* ntuple avatars 0/1=NO/Yes
4  2.51  100.  58.17	* Z, rms_R, rms_P, Binding_E
4 2.36 100. 58.17 10.	* Z, rms_R, rms_P, Binding_E, Amin for breakup
6 2.44 100. 92.17 10.	* Z, rms_R, rms_P, Binding_E, Amin for breakup 
6 2.44 100. 92.17 10.	* Z, rms_R, rms_P, Binding_E, Amin for breakup 
EOF
}

#=====================================================================
runProgram() # Run INCL/ABLA program
#=====================================================================
# Arg_1 = program name
{
$1
#rm -f cu43.in
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
cascade_evaporation_program="./cugnon43_khsv3p"
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
