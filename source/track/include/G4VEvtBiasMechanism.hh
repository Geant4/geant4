//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VEvtBiasMechanism.hh,v 1.4 2001-07-11 10:08:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//   
// -------------------------------------------------------
//  Class Description   
//   This class is a base class for all "event biasing mechanism". 
//  
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
public: // with description

  // constructors
  G4VEvtBiasMechanism(const G4String& name = "");
  G4VEvtBiasMechanism(const G4VEvtBiasMechanism& right);

public: 
  //destructors
  virtual ~G4VEvtBiasMechanism();


public: // with description
  // pure virtual functions
  
  // ApplyMath method will be invoked in G4VParticleChange::UpdateStepInfo()
  // if G4VParticleChange::fUseEB is set
  virtual G4VParticleChange* ApplyMath( G4VParticleChange*, const G4Step& ) =0;
     
  // IsApplicable method returns 'true' if this Event Bias Mechanism is 
  // valid for the particle type
  virtual G4bool IsApplicable(G4ParticleDefinition*) const = 0;

  
  // name of the biasing mechanism
  const G4String& GetName(){ return theEBName;}

  // Set/Get Verbose level 
  G4int GetVerboseLevel() const { return verboseLevel; }
  void  SetVerboseLevel(G4int value) { verboseLevel = value; }

 private:
  G4String theEBName;

 private:
  G4int verboseLevel;
  
};

#endif

