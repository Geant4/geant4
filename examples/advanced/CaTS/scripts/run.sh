#!/bin/bash -l 
msg="=== $BASH_SOURCE :"

./check.sh 
[ $? -ne 0 ] && echo $msg check failed && exit 1

extra_logging(){
    # switch on extra logging for these classes
    export G4Opticks=INFO
    export OpMgr=INFO
    export OpPropagator=INFO
    export OpEngine=INFO
    export OScene=INFO
    export OEvent=INFO
    export OSensorLib=INFO
}
#extra_logging

genstep_skips(){
    # from okc:OpticksGentep.h
    local OpticksGenstep_G4Cerenkov_1042=1
    local OpticksGenstep_G4Scintillation_1042=2
    local OpticksGenstep_DsG4Cerenkov_r3971=3
    local OpticksGenstep_DsG4Scintillation_r3971=4
    local OpticksGenstep_TORCH=5
    # uncomment one of the below to skip cerenkov OR scintillation gensteps at G4Opticks collection
    #export OPTICKS_SKIP_GENCODE=${OpticksGenstep_G4Cerenkov_1042},${OpticksGenstep_DsG4Cerenkov_r3971}
    #export OPTICKS_SKIP_GENCODE=${OpticksGenstep_G4Scintillation_1042},${OpticksGenstep_DsG4Scintillation_r3971}
}
genstep_skips


# defines the embedded commandline picking between production and devlopment(with lots of .npy saves)
#export OPTICKS_EMBEDDED_COMMANDLINE="dev"
export OPTICKS_EMBEDDED_COMMANDLINE="pro"

# appends to the embedded commandline
#export OPTICKS_EMBEDDED_COMMANDLINE_EXTRA="--rngmax 10 --trivial"   # trivial changes the generate.cu kernel 
export OPTICKS_EMBEDDED_COMMANDLINE_EXTRA="--rngmax 100"


#SY=50
#SY=100
SY=5000
#SY=10000
#SY=50000

cmd="CaTS $PWD/gdml/G4Opticks_$SY.gdml macros/muon_noIO.mac"
#cmd="gdb -ex r --args $cmd"   # uncomment to run under debugger

echo $cmd
eval $cmd 
