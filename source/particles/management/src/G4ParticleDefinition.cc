// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleDefinition.cc,v 1.2 1999-04-13 08:00:29 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:
//	CERN, IT Division, ASD Group
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
//      ---------------- G4ParticleDefinition -----------------
//      first implementation by Makoto Asai, 29 January 1996
//      revised by G.Cosmo, 29 February 1996
//      revised by H.Kurashige, 19 April 1996
//      Code uses operators (+=, *=, ++, -> etc.) correctly, P. Urban, 26/6/96
//      revised by H.Kurashige, 4 July 1996
//      revised by H.Kurashige, 16 Feb 1997
//      revised by H.Kurashige, 10 Nov 1997
//      remove new/delete G4ProcessManager   by H.Kurashige  06 June 1998 
//      added  Resonance flag and ApplyCuts flag  H.Kurashige 27  June 1998
//      modify FillQuarkContents() for quarks/diquarks H.Kurashige 30 June 1998
//      modify encoding rule H.Kurashige 23 Oct. 98
//      commented out G4cout/G4cerr in the constructor 10 Nov.,98 H.Kurashige
//       modify FillQuarkContents() for deltas         25 Nov.,98 H.Kurashige
// --------------------------------------------------------------


#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"

G4ParticleDefinition::G4ParticleDefinition(
		     const G4String&     aName,  
		     G4double            mass,
		     G4double            width,
                     G4double            charge,   
		     G4int               iSpin,
                     G4int               iParity,    
		     G4int               iConjugation,
                     G4int               iIsospin,   
		     G4int               iIsospin3, 
		     G4int               gParity,
		     const G4String&     pType,
                     G4int               lepton,      
		     G4int               baryon,
		     G4int               encoding,
		     G4bool              stable,
		     G4double            lifetime,
		     G4DecayTable        *decaytable,
		     G4bool              shortlived)
		 : theParticleName(aName), 
		   thePDGMass(mass),
		   thePDGWidth(width),
		   thePDGCharge(charge),
		   thePDGiSpin(iSpin),
		   thePDGSpin(iSpin*0.5),
		   thePDGiParity(iParity), 
		   thePDGiConjugation(iConjugation),
		   thePDGiGParity(gParity),
		   thePDGiIsospin(iIsospin),
		   thePDGiIsospin3(iIsospin3),
		   thePDGIsospin(iIsospin*0.5),
		   thePDGIsospin3(iIsospin3*0.5),
		   theParticleType(pType), 
		   theLeptonNumber(lepton),
		   theBaryonNumber(baryon),
		   thePDGEncoding(encoding),
		   theAntiPDGEncoding(-1*encoding),
		   thePDGStable(stable), 
		   thePDGLifeTime(lifetime), 
                   theDecayTable(decaytable),
		   theProcessManager(0),
		   fShortLivedFlag(shortlived),
		   fApplyCutsFlag(false),
                   verboseLevel(1)
{
   // check name and register this particle into ParticleTable
   theParticleTable = G4ParticleTable::GetParticleTable();
   G4ParticleDefinition *ptr = theParticleTable->Insert(this);
   //#ifdef G4VERBOSE
   //if (ptr != this) {
   //  // Fail to register to ParticleTable
   //  if (ptr != 0) {
   //    // particle with same name already exists in ParticleTable 
   //    if (verboseLevel>0) {
   //      G4cerr << "Particle name of " << aName << " has already defined " <<endl;
   //      G4cerr << "This particle is not registered in ParticleTable " << endl;
   //    } else {
   //	   if (verboseLevel>0) {
   //  	     G4cerr << "Particle name of " << aName; 
   //	     G4cerr << " is not registered in ParticleTable" <<endl; 
   //	     if (verboseLevel>1) this->DumpTable();
   //	   }
   //    }
   //  }
   //}
   //#endif

   // check quark contents
#ifdef G4VERBOSE
   if (this->FillQuarkContents() != thePDGEncoding) {
     if (verboseLevel>0) {
   //    G4cerr << "Particle " << aName << " has a strange PDGEncoding " <<endl;
       cerr << "Particle " << aName << " has a strange PDGEncoding " <<endl;
   //    if (verboseLevel>1) this->DumpTable();
     }
   }
#endif
}

G4ParticleDefinition::G4ParticleDefinition(const G4ParticleDefinition &right)
{
  G4Exception("You call Copy Constructor of G4ParticleDefinition ");
  //G4cerr << "You call Copy Constructor of G4ParticleDefinition " << endl;
  //#ifdef G4VERBOSE
  //  if (verboseLevel>0) right.DumpTable();
  //#endif
}

G4ParticleDefinition::G4ParticleDefinition()
{
  G4Exception("You call Default Constructor of G4ParticleDefinition ");
  //G4cerr << "You call Default Constructor of G4ParticleDefinition " << endl;
}


G4ParticleDefinition::~G4ParticleDefinition() 
{
  if (theDecayTable!= 0) delete theDecayTable;
}


const G4ParticleDefinition & G4ParticleDefinition::operator=(const G4ParticleDefinition &right)
{
  if (this != &right)  {
  } return right;
}

G4int G4ParticleDefinition::operator==(const G4ParticleDefinition &right) const
{
  return (this->theParticleName == right.theParticleName);
}

G4int G4ParticleDefinition::operator!=(const G4ParticleDefinition &right) const
{
  return (this->theParticleName != right.theParticleName);
}



G4int G4ParticleDefinition::FillQuarkContents()
      //  calculate quark and anti-quark contents
      //  return value is PDG encoding for this particle.
      //  It means error if the return value is differnt from
      //  this->thePDGEncoding.
{
  G4int tempPDGcode = thePDGEncoding;

  for (G4int flavor=0; flavor<NumberOfQuarkFlavor; flavor++){
    theQuarkContent[flavor] =0;
    theAntiQuarkContent[flavor] =0;
  }

  G4int temp = abs(tempPDGcode);
  
  G4int higherSpin = temp/10000000;
  temp -= G4int(higherSpin*10000000);
  G4int exotic = temp/1000000;
  temp -= G4int(exotic*1000000);
  G4int radial = temp/100000;
  temp -= G4int(radial*100000);
  G4int multiplet = temp/10000;
  temp -= G4int(multiplet*10000);
  G4int quark1 = temp/1000;
  temp -= G4int(quark1*1000);
  G4int quark2 = temp/100;
  temp -= G4int(quark2*100);
  G4int quark3 = temp/10;
  temp -= G4int(quark3*10);
  G4int spin= temp;
  if ((spin ==0) && ( higherSpin !=0 )) {
    spin =  higherSpin-1;
  } else {
    spin -= 1;
  }
  
  if (theParticleType =="quarks") {
    if ((quark1 !=0) || (quark2 !=0) || (quark3 !=0)) { 
      //#ifdef G4VERBOSE
      //if (verboseLevel>0) {
      //  G4cerr << " ??? unknown quark ";
      //  G4cerr << " PDG code=" << thePDGEncoding <<endl;
      //}
      //#endif
	//  --- thePDGEncoding is wrong 
	tempPDGcode = 0;
    } else {
      quark1 = abs(tempPDGcode);
      if (quark1>NumberOfQuarkFlavor){
	//#ifdef G4VERBOSE
	//if (verboseLevel>0) {
	//  G4cerr << " ??? unknown quark ";
	//  G4cerr << " PDG code=" << thePDGEncoding <<endl;
	//}
	//#endif
	//  --- thePDGEncoding is wrong 
	tempPDGcode = 0;
      } else {
	if (tempPDGcode>0){
	  theQuarkContent[quark1-1] =1;
	} else {
	  theAntiQuarkContent[quark1-1] =1;
	}
      }
    }

  } else if (theParticleType =="diquarks") {
    if ((quark1 ==0) || (quark2 ==0) || (quark3 !=0)) {
      // quark3 should be 0
      //  --- thePDGEncoding is wrong 
      tempPDGcode = 0;
    } else if (quark1 < quark2) {
      //  --- thePDGEncoding is wrong 
      tempPDGcode = 0;
    } else if (quark2>NumberOfQuarkFlavor){
      //#ifdef G4VERBOSE
      //if (verboseLevel>0) {
      //  G4cerr << " ??? unknown quark ";
      //  G4cerr << " PDG code=" << thePDGEncoding <<endl;
      //}
      //#endif
      tempPDGcode = 0;
    } else {
      if (tempPDGcode>0){
	theQuarkContent[quark1-1] +=1;
	theQuarkContent[quark2-1] +=1;
      } else {
	theAntiQuarkContent[quark1-1] +=1;
	theAntiQuarkContent[quark2-1] +=1;
      }
    }

  } else if (theParticleType =="gluons") {
    // gluons 
    //   do not care about

  } else if (theParticleType == "meson") {
     // check meson or not

    //   -- exceptions --
    if (tempPDGcode == 310) spin = 0;        //K0s
    if (tempPDGcode == 130) {     //K0l
      spin = 0;        
      quark2 = 3;
      quark3 = 1;
    }

    if (quark1 !=0) {
      //#ifdef G4VERBOSE
      //if (verboseLevel>0) {
      //  G4cerr << " meson has only quark and anti-quark pair";
      //  G4cerr << " PDG code=" << thePDGEncoding <<endl;
      //}
      //#endif
      tempPDGcode = 0;
    } 
     if ((quark2==0)||(quark3==0)){ 
       //#ifdef G4VERBOSE
       //if (verboseLevel>0) {
       // G4cerr << " meson has quark and anti-quark pair";
       // G4cerr << " PDG code=" << thePDGEncoding <<endl;
       //}
       //#endif
      tempPDGcode = 0;
     }
    // check spin 
    if ( spin != thePDGiSpin) {
      //#ifdef G4VERBOSE
      //if (verboseLevel>0) {
      //  G4cerr << " illegal SPIN (" << thePDGiSpin << "/2)";
      //  G4cerr << " PDG code=" << thePDGEncoding <<endl;
      //}
      //#endif
      tempPDGcode = 0;
    }
    if (quark2<quark3) { 
      //#ifdef G4VERBOSE
      //if (verboseLevel>0) {
      //  G4cerr << " illegal code for meson ";
      //  G4cerr << " PDG code=" << thePDGEncoding <<endl;
      //}
      //#endif
      tempPDGcode = 0;
    }
    // check quark flavor
    if (quark2> NumberOfQuarkFlavor){
      //#ifdef G4VERBOSE
      //if (verboseLevel>0) {
      //  G4cerr << " ??? unknown quark ";
      //  G4cerr << " PDG code=" << thePDGEncoding <<endl;
      //}
      //#endif
      tempPDGcode = 0;
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
    // check charge
    G4double totalCharge = 0.0;
    for (G4int flavor= 0; flavor<NumberOfQuarkFlavor-1; flavor+=2){
      totalCharge += (-1./3.)*eplus*theQuarkContent[flavor];
      totalCharge += 1./3.*eplus*theAntiQuarkContent[flavor];
      totalCharge += 2./3.*eplus*theQuarkContent[flavor+1];
      totalCharge += (-2./3.)*eplus*theAntiQuarkContent[flavor+1];
    }
    if (abs(totalCharge-thePDGCharge)>0.1*eplus) { 
      //#ifdef G4VERBOSE
      //if (verboseLevel>0) {
      //  G4cerr << " illegal charge for meson " << thePDGCharge/eplus;
      //  G4cerr << " PDG code=" << thePDGEncoding <<endl;
      //}
      //#endif
      tempPDGcode = 0;
    }
  } else if (theParticleType == "baryon"){
    // check meson or not
    if ((quark1==0)||(quark2==0)||(quark3==0)){ 
      //#ifdef G4VERBOSE
      //if (verboseLevel>0) {
      //  G4cerr << " meson has three quark ";
      //  G4cerr << " PDG code=" << thePDGEncoding <<endl;
      //}
      //#endif
      tempPDGcode = 0;
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

    // check spin 
    if (spin != thePDGiSpin) {
      //#ifdef G4VERBOSE
      //if (verboseLevel>0) {
      //  G4cerr << " illegal SPIN (" << thePDGiSpin << "/2";
      //  G4cerr << " PDG code=" << thePDGEncoding <<endl;
      //}
      //#endif
      tempPDGcode = 0;
    }
    // check quark flavor
    if ((quark1<quark2)||(quark2<quark3)||(quark1<quark3)) { 
      //#ifdef G4VERBOSE
      //if (verboseLevel>0) {
      // G4cerr << " illegal code for baryon ";
      //  G4cerr << " PDG code=" << thePDGEncoding <<endl;
      //}
      //#endif
      tempPDGcode = 0;
    }
    if (quark1> NumberOfQuarkFlavor) {
      //#ifdef G4VERBOSE
      //if (verboseLevel>0) {
      //  G4cerr << " ??? unknown quark ";
      // G4cerr << " PDG code=" << thePDGEncoding <<endl;
      //}
      //#endif
      tempPDGcode = 0;
    }
    if (tempPDGcode >0) {
      theQuarkContent[quark1-1] ++;
      theQuarkContent[quark2-1] ++;
      theQuarkContent[quark3-1] ++;
    } else {
      theAntiQuarkContent[quark1-1] ++;
      theAntiQuarkContent[quark2-1] ++;
      theAntiQuarkContent[quark3-1] ++;
    }
    // check charge 
    G4double totalCharge = 0.0;
    for (G4int flavor= 0; flavor<NumberOfQuarkFlavor-1; flavor+=2){
      totalCharge += (-1./3.)*eplus*theQuarkContent[flavor];
      totalCharge += 1./3.*eplus*theAntiQuarkContent[flavor];
      totalCharge += 2./3.*eplus*theQuarkContent[flavor+1];
      totalCharge += (-2./3.)*eplus*theAntiQuarkContent[flavor+1];
    }
    if (abs(totalCharge-thePDGCharge)>0.1*eplus) { 
      //#ifdef G4VERBOSE
      //if (verboseLevel>0) {
      //  G4cerr << " illegal charge for baryon " << thePDGCharge/eplus;
      //  G4cerr << " PDG code=" << thePDGEncoding <<endl;
      //}
      //#endif
      tempPDGcode = 0;
    }
  } else {
  }
  return tempPDGcode;
}

void G4ParticleDefinition::DumpTable() const
{
  G4cout << endl;
  G4cout << "--- G4ParticleDefinition ---" << endl;
  G4cout << " Particle Name : " << theParticleName << endl;
  G4cout << " PDG particle code : " << thePDGEncoding;
  G4cout << " [PDG anti-particle code: " << this->GetAntiPDGEncoding() << "]"<< endl;
  G4cout << " Mass [GeV/c2] : " << thePDGMass/GeV ;
  G4cout << "     Width : " << thePDGWidth/GeV << endl;
  G4cout << " Lifetime [nsec] : " << thePDGLifeTime/ns << endl;
  G4cout << " Charge [e]: " << thePDGCharge/eplus << endl;
  G4cout << " Spin : " << thePDGiSpin << "/2" << endl;
  G4cout << " Parity : " << thePDGiParity << endl;
  G4cout << " Charge conjugation : " << thePDGiConjugation << endl;
  G4cout << " Isospin : (I,Iz): (" << thePDGiIsospin <<"/2";
  G4cout << " , " << thePDGiIsospin3 << "/2 ) " << endl;
  G4cout << " GParity : " << thePDGiGParity << endl;
  G4cout << " Quark contents     (d,u,s,c,b,t) : " << theQuarkContent[0];
  G4cout << ", " << theQuarkContent[1];
  G4cout << ", " << theQuarkContent[2];
  G4cout << ", " << theQuarkContent[3];
  G4cout << ", " << theQuarkContent[4];
  G4cout << ", " << theQuarkContent[5] << endl;
  G4cout << " AntiQuark contents               : " << theAntiQuarkContent[0];
  G4cout << ", " << theAntiQuarkContent[1];
  G4cout << ", " << theAntiQuarkContent[2];
  G4cout << ", " << theAntiQuarkContent[3];
  G4cout << ", " << theAntiQuarkContent[4];
  G4cout << ", " << theAntiQuarkContent[5] << endl;
  G4cout << " Lepton number : " << theLeptonNumber;
  G4cout << " Baryon number : " << theBaryonNumber << endl;
  G4cout << " Particle type : " << theParticleType << endl;

  if ( fShortLivedFlag ){
    G4cout << " ShortLived : ON" << endl;
  }

  if ( thePDGStable ){
    G4cout << " Stable : stable" << endl;
  } else {
    if( theDecayTable != 0 ){
      theDecayTable->DumpInfo();
    } else {
      G4cout << "Decay Table is not defined !!" <<endl;
    }
  }

  if ( fApplyCutsFlag ){
    G4cout << " ApplyCuts : ON" << endl;
  } else {
    G4cout << " ApplyCuts : OFF" << endl;
  }
}











