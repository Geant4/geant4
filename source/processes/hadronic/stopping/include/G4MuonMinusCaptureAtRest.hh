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
// $Id$
//
//  G4MuonMinusCaptureAtRest physics process
//  Larry Felawka (TRIUMF) and Art Olin (TRIUMF) April 1998
//---------------------------------------------------------------------
//
// Class Description:
//
// Process for nuclear capture of mu- at rest takes into account Fermi 
// model of muon capture in compounds, simplified EM cascade model, muon 
// decay from K-shell, and muon nucleus reaction. 

//
// Modifications: 
// 18/08/2000  V.Ivanchenko Update description, new method to simulate capture
// 12/12/2003  H.P.Wellisch Completly rewrite mu-nuclear part
// 17/05/2006  V.Ivanchenko Cleanup
// 14/11/2006  V.Ivanchenko Remove implementation of GetPhysicsInteractionLength
//
//-----------------------------------------------------------------------------

#ifndef G4MuonMinusCaptureAtRest_h
#define G4MuonMinusCaptureAtRest_h 1
 
#include "globals.hh"
#include "G4VRestProcess.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4StopElementSelector.hh"
#include "G4MuMinusCaptureCascade.hh"
#include "G4ReactionProductVector.hh"
#include "G4MuonMinus.hh"
#include "G4HadronicProcessType.hh"

class G4Fancy3DNucleus;
class G4ExcitationHandler;
class G4GHEKinematicsVector;

class G4MuonMinusCaptureAtRest : public G4VRestProcess
 
{ 
public:
 
  G4MuonMinusCaptureAtRest(const G4String& processName ="muMinusCaptureAtRest", 
			   G4ProcessType   aType = fHadronic );

  virtual ~G4MuonMinusCaptureAtRest();

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  virtual void PreparePhysicsTable(const G4ParticleDefinition&);

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&); 

  G4double GetMeanLifeTime(const G4Track&, G4ForceCondition*) 
  {return 0;};

private:

  G4ReactionProductVector* DoMuCapture();

  // hide assignment operator as private 
  G4MuonMinusCaptureAtRest& operator=(const G4MuonMinusCaptureAtRest &right);
  G4MuonMinusCaptureAtRest(const G4MuonMinusCaptureAtRest& );
   
  G4int      nCascade;
  G4double   targetZ;
  G4double   targetA;
  G4double   targetMass;

  G4StopElementSelector*   pSelector;
  G4MuMinusCaptureCascade* pEMCascade;
  G4GHEKinematicsVector*   Cascade;
  G4Fancy3DNucleus*        theN;
  G4ExcitationHandler*     theHandler;

  G4bool isInitialised;

};

#endif
 





