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
// $Id: G4MuonMinusAtomicCapture.hh 66367 2012-12-18 09:18:08Z gcosmo $
//
//---------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4MuonMinusAtomicCapture
//
// 20160701 K.L. Genser - New process using G4MuonicAtom
//
// Class Description:
//
// Stopping of mu-
//
// Modifications: 
//   20160912 K.L. Genser made it rest process
//
//------------------------------------------------------------------------

#ifndef G4MuonMinusAtomicCapture_h
#define G4MuonMinusAtomicCapture_h 1
 
#include "globals.hh"
#include "G4VRestProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4ElementSelector.hh"
#include "G4HadronicInteraction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ForceCondition.hh"
#include "G4HadronicProcessType.hh"
#include "G4HadFinalState.hh"

class G4HadronicInteraction;

class G4MuonMinusAtomicCapture : public G4VRestProcess

{ 
public:
 
  explicit G4MuonMinusAtomicCapture(const G4String& name
                                    = "muMinusAtomicCaptureAtRest" );

  ~G4MuonMinusAtomicCapture();

  G4bool IsApplicable(const G4ParticleDefinition&);

  virtual void PreparePhysicsTable(const G4ParticleDefinition&);

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual G4double 
  AtRestGetPhysicalInteractionLength(const G4Track& track,
				     G4ForceCondition* condition);

  virtual G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&);

  void ProcessDescription(std::ostream& outFile) const;

  inline void SetElementSelector(G4ElementSelector* ptr);

  inline void SetEmCascade(G4HadronicInteraction* ptr);

protected:
  // set effective lifetime for at-rest process (default is forced action)
  // FIXME: This should be computed by subprocesses via cross-section analogue
  G4double GetMeanLifeTime(const G4Track& /*aTrack*/,
                           G4ForceCondition* /*condition*/) { return -1.0; }

private:

  // hide assignment operator as private 
  G4MuonMinusAtomicCapture& operator=(const G4MuonMinusAtomicCapture &right);
  G4MuonMinusAtomicCapture(const G4MuonMinusAtomicCapture& );

  G4ElementSelector* fElementSelector;

  G4HadronicInteraction* fEmCascade;

  G4ParticleChange* theTotalResult;

  G4HadFinalState* result;

  G4HadProjectile thePro;

  G4Nucleus targetNucleus;

};

#endif
 





