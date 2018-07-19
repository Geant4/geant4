//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4PDGCodeChecker.cc 105720 2017-08-16 12:38:10Z gcosmo $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      17 Aug 1999 H.Kurashige
// **********************************************************************

#include <fstream>
#include <iomanip>

#include "G4PDGCodeChecker.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

/////////////
G4PDGCodeChecker::G4PDGCodeChecker()
  :code(0),theParticleType(""),
   higherSpin(0),
   exotic(0),radial(0),multiplet(0),
   quark1(0),quark2(0),quark3(0),spin(0)
{
  verboseLevel = 1;
  // clear QuarkContents
  G4int flavor;
  for (flavor=0; flavor<NumberOfQuarkFlavor; flavor++){
    theQuarkContent[flavor] =0;
    theAntiQuarkContent[flavor] =0;
  }
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

  // check code for nuclei
  if ((theParticleType == "nucleus")||(theParticleType == "anti_nucleus")) {
    return CheckForNuclei();
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
    if (verboseLevel>0) {
      G4cout << " G4PDGCodeChecker::CheckPDGCode : ";
      G4cout << " meson has three quark ";
      G4cout << " PDG code=" << code <<G4endl;
    }
#endif
    return 0;
  }
 
 //exceptions
  if (std::abs(tempPDGcode)%10000 == 3122) { 
    // Lambda
    quark2=2;  quark3 = 1; spin = 1;
  } else if (std::abs(tempPDGcode)%10000 == 3124) { 
    // Lambda*
    quark2=2;  quark3 = 1; spin = 3;
  } else if (std::abs(tempPDGcode)%10000 == 3126) { 
    // Lambda*
    quark2=2;  quark3 = 1; spin = 5;
  } else if (std::abs(tempPDGcode)%10000 == 3128) { 
    // Lambda*
    quark2=2;  quark3 = 1; spin = 7;
  } else if (std::abs(tempPDGcode)%10000 == 4122) { 
    // Lambda_c
    quark2=2;  quark3 = 1; spin = 1;
  } else if (std::abs(tempPDGcode)%10000 == 5122) { 
    // Lambda_b
    quark2=2;  quark3 = 1; spin = 1;
  } else if (std::abs(tempPDGcode)%10000 == 4132) { 
    // Xi_c0
    quark2=3;  quark3 = 1; spin = 1;
  } else if (std::abs(tempPDGcode)%10000 == 4232) { 
    // Xi_c+
    quark2=3;  quark3 = 2; spin = 1;
  } else if (std::abs(tempPDGcode)%10000 == 5132) { 
    // Xi_b0
    quark2=3;  quark3 = 1; spin = 1;
  } else if (std::abs(tempPDGcode)%10000 == 5232) { 
    // Xi_b+
    quark2=3;  quark3 = 2; spin = 1;
  } else if (std::abs(tempPDGcode)%10000 == 2122) { 
    // Delta+ (spin 1/2) 
    quark2=2;  quark3 = 1; spin = 1;
  } else if (std::abs(tempPDGcode)%10000 == 1212) { 
    // Delta0 (spin 1/2) 
    quark1=2;  quark2 = 1; spin = 1;
  } else if (std::abs(tempPDGcode)%10000 == 2126) { 
    // Delta+ (spin 5/2) 
    quark2=2;  quark3 = 1; spin = 5;
  } else if (std::abs(tempPDGcode)%10000 == 1216) { 
    // Delta0 (spin 5/2) 
    quark1=2;  quark2 = 1; spin = 5;
  } else if (std::abs(tempPDGcode)%10000 == 2128) { 
    // Delta+ (spin 7/2) 
    quark2=2;  quark3 = 1; spin = 7;
  } else if (std::abs(tempPDGcode)%10000 == 1218) { 
    // Delta0 (spin 7/2) 
    quark1=2;  quark2 = 1; spin = 7;
  } else if (std::abs(tempPDGcode)%10000 == 2124) { 
    // N*+ (spin 3/2) 
    quark2=2;  quark3 = 1; spin = 3;
  } else if (std::abs(tempPDGcode)%10000 == 1214) { 
    // N*0 (spin 3/2) 
    quark1=2;  quark2 = 1; spin = 3;
  } 

    // check quark flavor
  if ((quark1<quark2)||(quark2<quark3)||(quark1<quark3)) { 
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << " G4PDGCodeChecker::CheckPDGCode : ";
      G4cout << " illegal code for baryon ";
      G4cout << " PDG code=" << code <<G4endl;
    }
#endif
    return 0;
  }
  if (quark1> NumberOfQuarkFlavor) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << " G4PDGCodeChecker::CheckPDGCode : ";
      G4cout << " ??? unknown quark ";
      G4cout << " PDG code=" << code <<G4endl;
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
    if (verboseLevel>0) {
      G4cout << " G4PDGCodeChecker::CheckPDGCode : ";
      G4cout << " meson has only quark and anti-quark pair";
      G4cout << " PDG code=" << code <<G4endl;
    }
#endif
    return 0;
  } 
  if (quark2<quark3) { 
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << " G4PDGCodeChecker::CheckPDGCode : ";
      G4cout << " illegal code for meson ";
      G4cout << " PDG code=" << code <<G4endl;
    }
#endif
    return 0;
  }

  // check quark flavor
  if (quark2> NumberOfQuarkFlavor){
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << " G4PDGCodeChecker::CheckPDGCode : ";
      G4cout << " ??? unknown quark ";
      G4cout << " PDG code=" << code <<G4endl;
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
    if (verboseLevel>0) {
      G4cout << " G4PDGCodeChecker::CheckPDGCode : ";
      G4cout << " ??? unknown quark ";
      G4cout << " PDG code=" << code <<G4endl;
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
  quark1 = std::abs(code);

  if ( std::abs(quark1)>NumberOfQuarkFlavor ) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << " G4PDGCodeChecker::CheckPDGCode : ";
      G4cout << " ??? unknown quark ";
      G4cout << " PDG code=" << code <<G4endl;
    }
#endif
    //  --- code is wrong 
    return 0;

  } 

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

  if (std::fabs(totalCharge-thePDGCharge)>0.1*eplus) { 
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << " G4PDGCodeChecker::CheckCharge  : ";
      G4cout << " illegal electric charge " << thePDGCharge/eplus;
      G4cout << " PDG code=" << code <<G4endl;
    }
#endif
    return false;
  }
  return true;
}

/////////////
G4int G4PDGCodeChecker::CheckForNuclei()
{
  G4int pcode = std::abs(code);
  if (pcode < 1000000000) {
    // non-nuclei   
    return 0;
  }

  pcode -= 1000000000;
  G4int LL = pcode/10000000;
  pcode -= 10000000*LL;
  G4int Z = pcode/10000;
  pcode -= 10000*Z;
  G4int A = pcode/10;
  
  // Allow neutron balls
  // if (A < 2 || Z > A-LL || LL>A || Z<=0 ) {
  if (A < 2 || Z > A-LL || LL>A ) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << " G4PDGCodeChecker::CheckPDGCode : ";
      G4cout << " ???  Illegal PDG encoding for nucleus ";
      G4cout << " PDG code=" << code <<G4endl;
    }
#endif
    return 0;
  }

  G4int n_up   = 2*Z +   (A-Z-LL) + LL;
  G4int n_down =   Z + 2*(A-Z-LL) + LL;
  G4int n_s    =   LL;

  // Fill Quark contents
  if (code>0) {
    theQuarkContent[0] = n_up;
    theQuarkContent[1] = n_down;
    theQuarkContent[2] = n_s;
   } else {
    // anti_nucleus
    theAntiQuarkContent[0] = n_up;
    theAntiQuarkContent[1] = n_down;
    theAntiQuarkContent[2] = n_s;
  }
  return code;
}
 
/////////////
void G4PDGCodeChecker::GetDigits(G4int PDGcode)
{
  G4int temp = std::abs(PDGcode);
  
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
