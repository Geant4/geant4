#!/bin/awk -f
##
##  Create summary report for processed output 
##     gzcat out-proc.en100GeV-it1k.n2-part3.gz | grep -v Elec  | ./report.awk > out-GammaNucl..
##     gzcat out-proc.en100GeV-it1k.n2-part3.gz | grep -v Photo | ./report.awk > out-ElecNucl...
##
##   John Apostolakis,  12 June 2012
BEGIN { allEnergy=0; }
/gun/                  { 
    if( newConfig || ((nEv>0) || (nBad>0)) ) { 
       report(); 
       newConfig=0; 
       resetAll(); 
    }
    # print; 
}  
/gun\/energy/        { primaryEn= $2;   unit = $3; newConfig=1; } 
/gun\/particle/      { primaryType= $2;            newConfig=1; } 
/[sS]electMaterial/  { material= $2;       }
/Event/              { nEv++; newConfig=0; }
! ( /gun/ || /Event/ || /Material/ || /Process/ ) { 
    if( /gamma/ || /e[+-]/ ) {
	nBad++;
	enDeficit=$10; pDeficit=$11; 
	newConfig=1; 
	# print  "Line of data", nBad, enDeficit, newConfig;

        if( enDeficit > 0.0 ) {
            nBad_E++; 
	    sumE  += enDeficit;
	    sumE2 += enDeficit*enDeficit;
	    abs_enDeficit= fabs(enDeficit);
	    if( abs_enDeficit>eDMaxAbs){   ## Must compare with Absolute !
	       eDAbs3rd=eDAbs2nd;
	       eDAbs2nd=eDMaxVal;
	       eDMaxVal=enDeficit;   
	       eDMaxAbs=fabs(enDeficit); 
	       maxline=$0;
	       # print "New Max> ", enDeficit, maxline; 
	    } 
	    else if( abs_enDeficit>fabs(eDAbs2nd)){   ## Must compare with Absolute !
	       eDAbs3rd= eDAbs2nd;
	       eDAbs2nd= enDeficit;
	    } 
	    else if( abs_enDeficit>fabs(eDAbs3rd)){   ## Must compare with Absolute !
	       eDAbs3rd= enDeficit;
	    }
	    if( enDeficit < eDMin){ eDMin2nd= eDMin;  eDMin = enDeficit;  } 
	    else {
	       if ( enDeficit<eDMin2nd){ eDMin2nd= enDeficit;  } 
	    }
	    ##  Now keep truly largest
	    if( enDeficit>eDMax){ eDMax2nd= eDMax; eDMax = enDeficit; } 
	    else {
	       if ( enDeficit>eDMax2nd){ eDMax2nd= enDeficit;  } 
	    }
	    if( fabs(enDeficit) > 0.5 * primaryEn * 1000 ) {   ## enDeficit is in MeV , primary in GeV
	       nOverHalf++;
	    } 
	} 

	if( pDeficit != 0.0 ){
            nBad_p++; 
	    sumP  += pDeficit;
	    sumP2 += pDeficit*pDeficit;
	    abs_pDeficit= fabs(pDeficit);
	    if( abs_pDeficit>fabs(pDMaxVal)){   ## Must compare with Absolute !
	       pDMaxVal=pDeficit;   
	    } 
	    if( pDeficit < pDMin){ pDMin = pDeficit;  } 
	    if( pDeficit > pDMax){ pDMax = pDeficit; } 
        }
    }else{
        print "Unknown type of line", $0;
    } 
}  
## /gun/ && newConfig     { report(); newConfig=0; } 
 END                   { report(); } 

func report(){ 
     #rint "Max= ", enDeficit, " Line> ", maxline; 
     # if( nBad  ==0 ) {
        # print "No bad events found. " 
     # }
     if( nBad > 0 ) {
        # print "Energy"; 
        if( nBad_E != 0.0) {
	    averE  = sumE  / nBad_E;
	    averE2 = sumE2 / nBad_E;
        }else{
	    averE  = 0.0 ; 
	    averE2 = 0.0
        }
        stdE2 = averE2 - averE * averE;
        sdvE  = stdE2  > 0 ? sqrt( stdE2 ) : 0.0; 
        rmsE  = averE2 > 0 ? sqrt( averE2 ) : 0.0; 
        #printf "%8s  %8s  %10.5f  %6d  %6d  %4d  %8.5g  %8.5g  %8.5g  %8.2f  %8.2f\n", material, primaryType, primaryEn, nBad, nOverHalf, nEv, eDMax, eD2nd, averE, sdvE, rmsE;   ## " Line> ", maxline; 
        printf "%8s  %8s  %10.5f  %6d  %6d  %4d  ", material, primaryType, primaryEn, nBad_E, nOverHalf, nEv; 
        if( eDMax > 1000.00 ) { 
            printf "%8.1f  ", eDMaxVal;
            if( allEnergy ){ 
               printf "%8.0f  %8.0f  %8.2f  %8.2f  %8.1f", eDMaxVal, eDAbs2nd, averE, sdvE, rmsE;   
               printf "%9.3f  %9.3f  ", eDMin,    eDMax; 
            }
        } else { 
            printf "%8.5g  ", eDMaxVal;
            if( allEnergy ){ 
               printf "%8.5g  %8.5g  %8.5f  %8.5f", eDAbs2nd, averE, sdvE, rmsE;   
               printf "%9.5g  %9.5g  ", eDMin,    eDMax; 
            }
        }
        # print "Momentum"; 
        if( nBad_p != 0.0) {
	    averP  = sumP  / nBad_p;
	    averP2 = sumP2 / nBad_p;
        }else{
	    averP  = averP2 =0.0; 
        }
        stdP2 = averP2 - averP * averP;
        sdvP  = stdP2  > 0 ? sqrt( stdP2 ) : -sqrt( -stdP2 ); // 0.0; 
        rmsP  = averP2 > 0 ? sqrt( averP2 ) : 0.0; 
        printf "%6d %8.5g  %8.5g  %8.5g  %8.5f", nBad_p, pDMaxVal, averP, sdvP, rmsP;   
        printf "%9.5g  %9.5g  ", pDMin,    pDMax; 
        printf "\n";
        ## print " Line> ", maxline; 
     }else if(nEv > 0){
        printf "%8s  %8s  %10.5f  %6d  %6d  %4d  %8.5g  %8.5g  %8.5g  %8.2f %8.2f   %8.2f  %8.2f\n", 
          material, primaryType, primaryEn, nBad, nOverHalf, nEv, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; 
     }
}

func resetAll(){
     nEv=0;
     nBad=0 ;
     nBad_E=0 
     eDMaxAbs= 0.0;
     eDMaxVal= 0.0;
     eDMax= -100000000.0;
     eDMin= 10.0e+25;
     eDAbs2nd=0.0; 
     eDAbs3rd=0.0; 
     sumE= sumE2 = 0.0; 
     maxline="";
     nOverHalf=0

     nBad_p=0 
     pDMaxVal= 0.0;
     sumP= sumP2 = 0.0; 
     pDMax= 0.0;
     pDMin= 0.0;
}

func fabs( x ) { return ((x >= 0) ? x : -x); } 

BEGIN{  resetAll();  
        #printf "# %8s %8s %10s  %6s  %7s  %4s  %8s  %8s  %8s  %8s\n", 
        #         "Material", "Primary", "E_primary", "nBad", "n>E_p/2", "nEv", "eDMax", "eD2nd", "<E>", "stddev(E) rms(E)";   ## " Line> ", maxline; 
        printf  "# Material Primary   E_primary    nBad   n>0.5E_p nEv   dE Max     ";# "nBadP   Max|deltaP|    <dP>    stddev(P),   rms(P)   p_Min,  pMax  "
        if( allEnergy ) printf " dE 2nd    <dE>    stddev(E)   rms(E)    dE min     dE max   "; 
        print   "nBadP   Max|deltaP|    <dP>    stddev(P),   rms(P)   p_Min,  pMax  "
        print   "#                      (GeV)                             (MeV)              (MeV)         (MeV)    (MeV)     (MeV)      (MeV)      (MeV)"; 
     }

#
# Process            Model                primary  primaryEn    nucleus,   limitR    energyR momentumR   limitA    energyA momentumA  chargeB  baryonB
# /gun/direction 1 0 0
# /gun/particle gamma
# /gun/energy 100 GeV
# ElectroNuclear     CHIPSElectroNuclear      e-  44754.60   (82,208),        0     0.3352    0.0000     5000   15001.80      0.00    0    0
