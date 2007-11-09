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
// $Id: G4DNAProcess.hh,v 1.5 2007-11-09 16:20:04 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author:  Maria Grazia Pia (Maria.Grazia.Pia@ge.infn.it)
//
// History:
// -----------
// 29 Sep 2007 MGP   First (incomplete) implementation
//
// -------------------------------------------------------------------

// Class description:
// Host class for DNA physics policies
// Further documentation available in PAPER REF. TO BE ADDED

// -------------------------------------------------------------------

#ifndef G4DNAPROCESS_HH
#define G4DNAPROCESS_HH 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4FinalStateProduct.hh"
#include <vector>

template < class TCrossSection, class TFinalState>
class G4DNAProcess : public G4VDiscreteProcess {
  
public:

  // ---- MGP ---- Note: process name to be replaced with a better identifying mechanism  
  G4DNAProcess(const G4String& processName = "DNAProcess"): G4VDiscreteProcess(processName) { /* nop */; }
  
  ~G4DNAProcess() { /* nop */; }

  //  ---- MGP ---- Dummy initially: process is always applicable
  virtual G4bool IsApplicable(const G4ParticleDefinition&) { return true; } 
  
  //  ---- MGP ---- Dummy initially: no PhysicsTable (usefulness to be verified)
  virtual void BuildPhysicsTable(const G4ParticleDefinition& /* particle */) { /* nop */; }
 
  // Implemented in terms of TFinalState
  virtual G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);
 
  // For testing purpose only
  G4double DumpMeanFreePath(const G4Track& aTrack, 
			    G4double previousStepSize, 
			    G4ForceCondition* condition) 
  { return GetMeanFreePath(aTrack, previousStepSize, condition); }

protected:
  
  // Implemented in terms of TCrossSection
  virtual G4double GetMeanFreePath(const G4Track& track, 
			   G4double previousStepSize, 
			   G4ForceCondition* condition);
private:

 // Hide copy constructor and assignment operator as private 
  G4DNAProcess& operator=(const G4DNAProcess& right);
  G4DNAProcess(const G4DNAProcess& );

  // Policy classes
  TCrossSection crossSection;
  TFinalState finalState;

};

#include "G4DNAProcess.icc"

#endif



