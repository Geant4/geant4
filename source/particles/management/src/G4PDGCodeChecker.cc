// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PDGCodeChecker.cc,v 1.2 1999-08-19 08:18:31 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, based on object model of
//      17 Aug 1999 H.Kurashige
// **********************************************************************

#include <fstream.h>
#include <iomanip.h>

#include "G4PDGCodeChecker.hh"

/////////////
G4PDGCodeChecker::G4PDGCodeChecker()
{
  code = 0;
  verboseLevel = 3;
}

/////////////
G4int  G4PDGCodeChecker::CheckPDGCode( G4int    PDGcode, 
				       G4String particleType)
{
  code = PDGcode;
  theParticleType = particleType;

  // clear QuarkContents
  G4int flavor;
  for (flavor=0; flavor<NumberOfQuarkFlavor; flavor++){
    theQuarkContent[flavor] =0;
    theAntiQuarkContent[flavor] =0;
  }

  // get each digit number
  GetDigits(code);

  // check code
  if (theParticleType =="quarks") {
    return CheckForQuarks();

  } else if  (theParticleType =="diquarks") {
    return CheckForDiQuarks();

  } else if (theParticleType =="gluons") {
    // gluons 
    //   do not care about
    return code;

  } else if (theParticleType == "meson") {
    return CheckForMesons();

  } else if (theParticleType == "baryon"){
    return CheckForBaryons();

  }
  // No check
  return code;
}
 
/////////////
G4int G4PDGCodeChecker::CheckForBaryons()
{
  G4int   tempPDGcode = code;

  if ((quark1==0)||(quark2==0)||(quark3==0)){ 
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cout << " meson has three quark ";
      G4cout << " PDG code=" << code <<endl;
    }
#endif
    return 0;
  }
 
 //exceptions
  if (abs(tempPDGcode)%10000 == 3122) { 
    // Lambda
    quark2=2;  quark3 = 1; spin = 1;
  } else if (abs(tempPDGcode)%10000 == 4122) { 
    // Lambda_c
    quark2=2;  quark3 = 1; spin = 1;
  } else if (abs(tempPDGcode)%10000 == 4132) { 
    // Xi_c0
    quark2=3;  quark3 = 1; spin = 1;
  } else if (abs(tempPDGcode)%10000 == 4232) { 
    // Xi_c+
    quark2=3;  quark3 = 2; spin = 1;
  } else if (abs(tempPDGcode)%10000 == 2122) { 
    // Delta+ (spin 1/2) 
    quark2=2;  quark3 = 1; spin = 1;
  } else if (abs(tempPDGcode)%10000 == 1212) { 
    // Delta0 (spin 1/2) 
    quark1=2;  quark2 = 1; spin = 1;
  } else if (abs(tempPDGcode)%10000 == 2126) { 
    // Delta+ (spin 5/2) 
    quark2=2;  quark3 = 1; spin = 5;
  } else if (abs(tempPDGcode)%10000 == 1216) { 
    // Delta0 (spin 5/2) 
    quark1=2;  quark2 = 1; spin = 5;
  } else if (abs(tempPDGcode)%10000 == 2128) { 
    // Delta+ (spin 7/2) 
    quark2=2;  quark3 = 1; spin = 7;
  } else if (abs(tempPDGcode)%10000 == 1218) { 
    // Delta0 (spin 7/2) 
    quark1=2;  quark2 = 1; spin = 7;
  } else if (abs(tempPDGcode)%10000 == 2124) { 
    // N*+ (spin 3/2) 
    quark2=2;  quark3 = 1; spin = 3;
  } else if (abs(tempPDGcode)%10000 == 1214) { 
    // N*0 (spin 3/2) 
    quark1=2;  quark2 = 1; spin = 3;
  } 

    // check quark flavor
  if ((quark1<quark2)||(quark2<quark3)||(quark1<quark3)) { 
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cout << " illegal code for baryon ";
      G4cout << " PDG code=" << code <<endl;
    }
#endif
    return 0;
  }
  if (quark1> NumberOfQuarkFlavor) {
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cout << " ??? unknown quark ";
      G4cout << " PDG code=" << code <<endl;
    }
#endif
    return 0;
  }
  

  // Fill Quark contents
  if (tempPDGcode >0) {
    theQuarkContent[quark1-1] ++;
    theQuarkContent[quark2-1] ++;
    theQuarkContent[quark3-1] ++;
  } else {
    theAntiQuarkContent[quark1-1] ++;
    theAntiQuarkContent[quark2-1] ++;
    theAntiQuarkContent[quark3-1] ++;
  }

  return code;
}
 
/////////////
G4int G4PDGCodeChecker::CheckForMesons()
{
  G4int   tempPDGcode = code;

  //   -- exceptions --
  if (tempPDGcode == 310) spin = 0;        //K0s
  if (tempPDGcode == 130) {     //K0l
    spin = 0;        
    quark2 = 3;
    quark3 = 1;
  }
  
  // 
  if ((quark1 !=0)||(quark2==0)||(quark3==0)){ 
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cout << " meson has only quark and anti-quark pair";
      G4cout << " PDG code=" << code <<endl;
    }
#endif
    return 0;
  } 
  if (quark2<quark3) { 
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cout << " illegal code for meson ";
      G4cout << " PDG code=" << code <<endl;
    }
#endif
    return 0;
  }

  // check quark flavor
  if (quark2> NumberOfQuarkFlavor){
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cout << " ??? unknown quark ";
      G4cout << " PDG code=" << code <<endl;
    }
#endif
    return 0;
  }

 
  // check heavier quark type
  if (quark2 & 1) {
    // down type qurak
    if (tempPDGcode >0) {
      theQuarkContent[quark3-1] =1;
      theAntiQuarkContent[quark2-1] =1;
    } else {
      theQuarkContent[quark2-1] =1;
      theAntiQuarkContent[quark3-1] =1;
    }
  } else {
    // up type quark
    if (tempPDGcode >0) {
      theQuarkContent[quark2-1] =1;
      theAntiQuarkContent[quark3-1] =1;
    } else {
      theQuarkContent[quark3-1] =1;
      theAntiQuarkContent[quark2-1] =1;
    }
  }
  return code;
}

  

/////////////
G4int G4PDGCodeChecker::CheckForDiQuarks()
{
  if ((quark1 ==0) || (quark2 ==0) || (quark3 !=0)) {
    // quark3 should be 0
    //  --- code is wrong 
    return 0;

  } else if (quark1 < quark2) {
    //  --- code is wrong 
    return 0;

  } else if (quark2>NumberOfQuarkFlavor){
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cout << " ??? unknown quark ";
      G4cout << " PDG code=" << code <<endl;
    }
#endif
    return 0;

  }

  // Fill Quark Contents
  if (code>0){
    theQuarkContent[quark1-1] +=1;
    theQuarkContent[quark2-1] +=1;
  } else {
    theAntiQuarkContent[quark1-1] +=1;
    theAntiQuarkContent[quark2-1] +=1;
  }

  return code;
}
 
/////////////
G4int G4PDGCodeChecker::CheckForQuarks()
{
  if ( abs(quark1)>NumberOfQuarkFlavor ) {
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cout << " ??? unknown quark ";
      G4cout << " PDG code=" << code <<endl;
    }
#endif
    //  --- code is wrong 
    return 0;

  } 

  quark1 = abs(code);

  // Fill Quark Contents
  if (code>0){
    theQuarkContent[quark1-1] =1;
  } else {
    theAntiQuarkContent[quark1-1] =1;
  }
  return code;
}

/////////////
G4bool G4PDGCodeChecker::CheckCharge(G4double thePDGCharge) const
{
  // check charge
  G4double totalCharge = 0.0;
  for (G4int flavor= 0; flavor<NumberOfQuarkFlavor-1; flavor+=2){
    totalCharge += (-1./3.)*eplus*theQuarkContent[flavor];
    totalCharge += 1./3.*eplus*theAntiQuarkContent[flavor];
    totalCharge += 2./3.*eplus*theQuarkContent[flavor+1];
    totalCharge += (-2./3.)*eplus*theAntiQuarkContent[flavor+1];
  }

  if (abs(totalCharge-thePDGCharge)>0.1*eplus) { 
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cout << " illegal electric charge " << thePDGCharge/eplus;
      G4cout << " PDG code=" << code <<endl;
    }
#endif
    return false;
  }
  return true;
}

/////////////
void G4PDGCodeChecker::GetDigits(G4int PDGcode)
{
  G4int temp = abs(PDGcode);
  
  higherSpin = temp/10000000;
  temp -= G4int(higherSpin*10000000);

  exotic = temp/1000000;
  temp -= G4int(exotic*1000000);

  radial = temp/100000;
  temp -= G4int(radial*100000);

  multiplet = temp/10000;
  temp -= G4int(multiplet*10000);

  quark1 = temp/1000;
  temp -= G4int(quark1*1000);

  quark2 = temp/100;
  temp -= G4int(quark2*100);

  quark3 = temp/10;
  temp -= G4int(quark3*10);

  spin= temp;
  if ((spin ==0) && ( higherSpin !=0 )) {
    spin =  higherSpin-1;
  } else {
    spin -= 1;
  }
}
