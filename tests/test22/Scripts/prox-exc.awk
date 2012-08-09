#!/usr/bin/awk -f 
#
#  NOTE:  The corresponding script in Perl prox-exc.pl, is derived from this, updated and faster.
#           John Apostolakis,  29 June 2012
#
#  Process the output of a run in which checking for Energy and momentum conservation is activated.
#  Summarise each issue reported in one of two ways : 
#      either from default check for 'catastrophic' violation of Energy conservation
#      or for problem with 'developer' level checks of E/p, charge and baryon number conservation.
#
#  NOTE: Since 9.6-beta both relative and absolute E/p checks must fail for an event to report a failure
#        (  Version before 9.5-ref-04 or so reported an issues if only one failed. )
#
#  First   version:   John Apostolakis,  31 May  2012   
#  Latest revision:   John Apostolakis,  15 June 2012  ( ca )
#
#  INPUT Format:
#rint process, model, primary, primaryEn, nucleus, limitR, energyR, momentumR, limitA, energyA, momentumA, chargeB, baryonB; 
#
#   1        2              3      4      5            6      7       8        9   10       11  12        13    14
#  Process: ElectroNuclear  ,      Model: CHIPSElectroNuclear
#  Primary: e-            (11),    E=     39.5593,     target nucleus (82,207)
#  fail     relative,     limit    0.01,  values       E      /       p     (MeV)  =  0.390098  /   0.385144
#  fail     absolute,     limit    (MeV)  1,           values E       /         p  MeV)      =  5.2327     /   15.2347
#  pass     charge/baryon number   balance 0 / 0 

/Process/ && /Model/ { process= $2; model= $5; } 
/Primary/   { primary=$2; primaryEn=$5; nucleus=$8; } ##            (11), E= 39.5593, target nucleus (82,207)
/relative/  { limitR=$4; energyR= $11; momentumR=$13; } 
/absolute/  { limitA=$5; energyA= $12; momentumA=$14; } 
/charge/  { chargeB=$5; baryonB=$7;  } 
/charge/  { report_it(); } 


## Problem is found in Excpetion
/G4Exception/ && /START/ { InException= 1; } 

## /Primary/ && InException { particle= $2; pdg= $11; En=$5;  nucleus=$8; } 
/initial/ && /final/ && InException {   deltaE= $5; unit= $6;  
         energyA= deltaE; 
         energyR= deltaE / primaryEn;  
         momentumA = momentumR = 0.0; 
         chargeB=0; baryonB=0; 
}  

##  E(initial - final) = 51787.2 MeV.
/Process/ && /Model/ && InException { process= $4; model= $6; 
         limitR= 0;  limitA= 5000; 
} 
##  Process / Model: ElectroNuclear / CHIPSElectroNuclear



/G4Exception/ &&  /END/ { 
	    ## print "*** From Exception *** "; 
	    report_it(); 
            InException= 0; 
} 

## E(initial - final) = 15001.8 MeV.

## -------- WWWW ------- G4Exception-START -------- WWWW -------
## *** G4Exception : had012
##       issued by : G4HadronicProcess:CheckResult()
## Warning: Bad energy non-conservation detected, will re-sample the interaction
##  Process / Model: ElectroNuclear / CHIPSElectroNuclear
##  Primary: e- (11), E= 44754.6, target nucleus (82,208)
##  E(initial - final) = 15001.8 MeV.
## 
## *** This is just a warning message. ***
## -------- WWWW -------- G4Exception-END --------- WWWW -------

 /Event/   { if( verbose ) { print; } } 
 /\/gun/ || /\gun\/energy/  { if( verbose ) { print; } } 
 /Material/  { if( verbose ) { print; } } 

BEGIN      { verbose=1; } 

func report_it() { 
    #print process, model, primary, primaryEn, nucleus, limitR, energyR, momentumR,  limitA, energyA, momentumA, chargeB, baryonB; 
    #print prox  mod  pr  prEn  nucl  limR enR,  pR,  limA  enA   momA  chargeB baryonB
    printf "%-18s %-20s %6s  %8.2f  %9s, %8s  %9.4f %9.4f %8s  %9.2f %9.2f  %3s  %3s\n",
            process, model, primary, primaryEn, nucleus, limitR, energyR, momentumR,  limitA, energyA, momentumA, chargeB, baryonB; 
    #rint prox mod           pr     prEn     nucl     limR  energyR, pR       lA energyA limA, chargeB, baryonB; 
}
# ElectroNuclear CHIPSElectroNuclear  30.5987, (82,208) 0.01, 0.352011 0.346255 1, 10.5912 10.5935 0 0
BEGIN{  printf "%-18s %-20s %6s  %8s  %9s, %8s  %9s %9s %8s  %9s %9s  %3s  %3s\n",
            "Process", "Model", "primary", "primaryEn", "nucleus", "limitR", "energyR", "momentumR",  "limitA", "energyA", "momentumA", "chargeB", "baryonB"; 
     }

    #rint process, model, primary, primaryEn, nucleus, limitR, energyR, momentumR, limitA, energyA, momentumA, chargeB, baryonB; 

#   1        2              3      4      5            6      7       8        9   10       11  12        13    14
#  Process: ElectroNuclear  ,      Model: CHIPSElectroNuclear
#  Primary: e-            (11),    E=     39.5593,     target nucleus (82,207)
#  fail     relative,     limit    0.01,  values       E      /       p     (MeV)  =  0.390098  /   0.385144
#  fail     absolute,     limit    (MeV)  1,           values E       /         p  MeV)      =  5.2327     /   15.2347
#  pass     charge/baryon number   balance 0 / 0 

## Original lines:
#  Process: ElectroNuclear , Model: CHIPSElectroNuclear
# Primary: e- (11), E= 39.5593, target nucleus (82,207)
#   fail relative, limit 0.01, values E / p (MeV) = 0.390098 / 0.385144
#   fail absolute, limit (MeV) 1, values E / p (MeV) = 15.2327 / 15.2347
#   pass charge/baryon number balance 0 / 0 
