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
// $Id$
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
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
//      modify FillQuarkContents() for deltas      25 Nov.,98 H.Kurashige
//
//      modify FillQuarkContents() to use G4PDGCodeChecker 17 Aug. 99 H.Kurashige
// --------------------------------------------------------------


#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4DecayTable.hh"
#include "G4PDGCodeChecker.hh"
#include "G4StateManager.hh"

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
		     G4bool              shortlived,
                     const G4String&     subType,
                     G4int               anti_encoding,
		     G4double            magneticMoment)

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
		   thePDGMagneticMoment(magneticMoment),
		   theLeptonNumber(lepton),
		   theBaryonNumber(baryon),
		   theParticleType(pType), 
		   theParticleSubType(subType), 
		   thePDGEncoding(encoding),
		   theAntiPDGEncoding(-1*encoding),
		   fShortLivedFlag(shortlived),
		   thePDGStable(stable), 
		   thePDGLifeTime(lifetime), 
                   theDecayTable(decaytable),
		   theProcessManager(0),
                   theAtomicNumber(0),
                   theAtomicMass(0),
                   verboseLevel(1),
  		   fApplyCutsFlag(false)
{
  static G4String nucleus("nucleus");
  theParticleTable = G4ParticleTable::GetParticleTable();
   
   //set verboseLevel equal to ParticleTable 
   verboseLevel = theParticleTable->GetVerboseLevel();

   if (anti_encoding !=0) theAntiPDGEncoding = anti_encoding;

   // check quark contents
   if (this->FillQuarkContents() != thePDGEncoding) {
#ifdef G4VERBOSE
     if (verboseLevel>0) {
       // Using G4cerr expecting that it is available in construction of static objects 
       G4cerr << "Particle " << aName << " has a strange PDGEncoding " <<G4endl;
     }
#endif
     G4Exception( "G4ParticleDefintion::G4ParticleDefintion",
		  "PART102", JustWarning, 
		  "Strange PDGEncoding ");
   }

   // check initialization is in Pre_Init state except for ions
   G4ApplicationState currentState = G4StateManager::GetStateManager()->GetCurrentState();

   if ( !fShortLivedFlag && (theParticleType!=nucleus) && (currentState!=G4State_PreInit)){
#ifdef G4VERBOSE
     if (GetVerboseLevel()>0) {
       G4cerr << "G4ParticleDefintion (other than ions and shortlived) should be created in Pre_Init state  " 
              << aName << G4endl;
     }
#endif
     G4Exception( "G4ParticleDefintion::G4ParticleDefintion",
		  "PART101", JustWarning, 
		  "G4ParticleDefinition should be created in PreInit state");
   }

   
   if (theParticleTable->GetIonTable()->IsIon(this)) {
     SetAtomicNumber( G4int(GetPDGCharge()/eplus) );
     SetAtomicMass( GetBaryonNumber() );
   }
  
   if (theParticleTable->GetIonTable()->IsAntiIon(this)) {
     SetAtomicNumber( std::abs(G4int(GetPDGCharge()/eplus)) );
     SetAtomicMass( std::abs(GetBaryonNumber()) );
   }
   
   // check name and register this particle into ParticleTable
   theParticleTable->Insert(this);

}

G4ParticleDefinition::G4ParticleDefinition(const G4ParticleDefinition &)
{
  G4Exception("G4ParticleDefinition::G4ParticleDefinition()",
	      "PART001", FatalException,
	      "Illegal call of copy Constructor for G4ParticleDefinition ");
}

G4ParticleDefinition::G4ParticleDefinition()
{
  G4Exception("G4ParticleDefinition::G4ParticleDefinition()",
	      "PART001", FatalException,
	      "Illegal call of default Constructor for G4ParticleDefinition ");
}


G4ParticleDefinition::~G4ParticleDefinition() 
{
  if (G4ParticleTable::GetParticleTable()->GetReadiness()) {
    G4StateManager* pStateManager = G4StateManager::GetStateManager();
    G4ApplicationState currentState = pStateManager->GetCurrentState();
    if (currentState != G4State_PreInit) {
      G4String msg = "Request of deletion for ";
      msg += GetParticleName();  
      msg += " has No effects because readyToUse is true.";
      G4Exception("G4ParticleDefinition::~G4ParticleDefinition()",
		  "PART117", JustWarning, msg);
      return ;
    } else {
#ifdef G4VERBOSE
      if (verboseLevel>0){
	G4cout << GetParticleName()
	       << " will be deleted " << G4endl;
      }
#endif
    }
  }

  if (theDecayTable!= 0) delete theDecayTable;
}


const G4ParticleDefinition & G4ParticleDefinition::operator=(const G4ParticleDefinition &right)
{
  if (this != &right)  {
  } 
  return *this;
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
  G4int flavor;
  for (flavor= 0; flavor<NumberOfQuarkFlavor; flavor++){
    theQuarkContent[flavor]     = 0;
    theAntiQuarkContent[flavor] = 0;
  }

  G4PDGCodeChecker checker;
  checker.SetVerboseLevel(verboseLevel);

  G4int temp = checker.CheckPDGCode(thePDGEncoding, theParticleType);

  if ( temp != 0) {
    for (flavor= 0; flavor<NumberOfQuarkFlavor; flavor++){
      theQuarkContent[flavor]     = checker.GetQuarkContent(flavor);
      theAntiQuarkContent[flavor] = checker.GetAntiQuarkContent(flavor);
    }
    if ((theParticleType == "meson")||(theParticleType == "baryon")) {
      // check charge
      if (!checker.CheckCharge(thePDGCharge) ){
	temp = 0;
	G4Exception( "G4ParticleDefintion::G4ParticleDefintion",
		  "PART103", JustWarning, 
		  "Inconsistent charge against PDG code ");
#ifdef G4VERBOSE
	if (verboseLevel>0) {
	  G4cerr << "G4ParticleDefinition::FillQuarkContents  : "
	         << " illegal charge (" << thePDGCharge/eplus
	         << " PDG code=" << thePDGEncoding <<G4endl;
	}
#endif
      }
      // check spin 
      if (checker.GetSpin() != thePDGiSpin) {
	temp=0;
	G4Exception( "G4ParticleDefintion::G4ParticleDefintion",
		  "PART104", JustWarning, 
		  "Inconsistent spin against PDG code ");
#ifdef G4VERBOSE
	if (verboseLevel>0) {
	  G4cerr << "G4ParticleDefinition::FillQuarkContents  : "
	         << " illegal SPIN (" << thePDGiSpin << "/2"
	         << " PDG code=" << thePDGEncoding <<G4endl;
	}
#endif
      }
    }
  }
  return temp;
}

void G4ParticleDefinition::DumpTable() const
{
  G4cout << G4endl;
  G4cout << "--- G4ParticleDefinition ---" << G4endl;
  G4cout << " Particle Name : " << theParticleName << G4endl;
  G4cout << " PDG particle code : " << thePDGEncoding;
  G4cout << " [PDG anti-particle code: " << this->GetAntiPDGEncoding() << "]"<< G4endl;
  G4cout << " Mass [GeV/c2] : " << thePDGMass/GeV ;
  G4cout << "     Width : " << thePDGWidth/GeV << G4endl;
  G4cout << " Lifetime [nsec] : " << thePDGLifeTime/ns << G4endl;
  G4cout << " Charge [e]: " << thePDGCharge/eplus << G4endl;
  G4cout << " Spin : " << thePDGiSpin << "/2" << G4endl;
  G4cout << " Parity : " << thePDGiParity << G4endl;
  G4cout << " Charge conjugation : " << thePDGiConjugation << G4endl;
  G4cout << " Isospin : (I,Iz): (" << thePDGiIsospin <<"/2";
  G4cout << " , " << thePDGiIsospin3 << "/2 ) " << G4endl;
  G4cout << " GParity : " << thePDGiGParity << G4endl;
  if (thePDGMagneticMoment != 0.0) {
    G4cout << " MagneticMoment [MeV/T] : " << thePDGMagneticMoment/MeV*tesla << G4endl;
  }
  G4cout << " Quark contents     (d,u,s,c,b,t) : " << theQuarkContent[0];
  G4cout << ", " << theQuarkContent[1];
  G4cout << ", " << theQuarkContent[2];
  G4cout << ", " << theQuarkContent[3];
  G4cout << ", " << theQuarkContent[4];
  G4cout << ", " << theQuarkContent[5] << G4endl;
  G4cout << " AntiQuark contents               : " << theAntiQuarkContent[0];
  G4cout << ", " << theAntiQuarkContent[1];
  G4cout << ", " << theAntiQuarkContent[2];
  G4cout << ", " << theAntiQuarkContent[3];
  G4cout << ", " << theAntiQuarkContent[4];
  G4cout << ", " << theAntiQuarkContent[5] << G4endl;
  G4cout << " Lepton number : " << theLeptonNumber;
  G4cout << " Baryon number : " << theBaryonNumber << G4endl;
  G4cout << " Particle type : " << theParticleType ;
  G4cout << " [" << theParticleSubType << "]" << G4endl;

  if (   (theParticleTable->GetIonTable()->IsIon(this)) 
      || (theParticleTable->GetIonTable()->IsAntiIon(this)) ) {
    G4cout << " Atomic Number : " << GetAtomicNumber();
    G4cout << "  Atomic Mass : " << GetAtomicMass()  << G4endl;
  }
  if ( fShortLivedFlag ){
    G4cout << " ShortLived : ON" << G4endl;
  }

  if ( thePDGStable ){
    G4cout << " Stable : stable" << G4endl;
  } else {
    if( theDecayTable != 0 ){
      theDecayTable->DumpInfo();
    } else {
      G4cout << "Decay Table is not defined !!" <<G4endl;
    }
  }

}

void G4ParticleDefinition::SetApplyCutsFlag(G4bool flg)
{
  if(theParticleName=="gamma"
  || theParticleName=="e-"
  || theParticleName=="e+"
  || theParticleName=="proton")
  { fApplyCutsFlag = flg; }
  else
  {
    G4cerr
     << "G4ParticleDefinition::SetApplyCutsFlag() for " << theParticleName
     << G4endl;
    G4cerr
     << "becomes obsolete. Production threshold is applied only for "
     << "gamma, e- ,e+ and proton." << G4endl;
  }
}
