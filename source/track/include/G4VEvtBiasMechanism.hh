// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VEvtBiasMechanism.hh,v 1.1 1999-01-07 16:14:24 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
// 
// ------------------------------------------------------------
//   Implemented for the new scheme            17 Nov. 1998  H.Kurahige
// 

#ifndef G4VEvtBiasMechanism_h
#define G4VEvtBiasMechanism_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4Step;
class G4VParticleChange;
class G4ParticleDefinition;

class G4VEvtBiasMechanism
{
 public:
  G4VEvtBiasMechanism(const G4String& name = ""){ theEBName = name;}
  G4VEvtBiasMechanism(const G4VEvtBiasMechanism& right){ theEBName = right.theEBName;}

  virtual ~G4VEvtBiasMechanism(){};

  virtual G4VParticleChange* ApplyMath( G4VParticleChange*, const G4Step& ) =0;
  virtual G4bool IsApplicable(G4ParticleDefinition*) const = 0;
  
  G4String GetName(){ return theEBName;}
  G4int GetVerboseLevel() { return verboseLevel; }
  void  SetVerboseLevel(G4int value) { verboseLevel = value; }
 private:
  G4String theEBName;

 private:
  G4int verboseLevel;
  
};

#endif

