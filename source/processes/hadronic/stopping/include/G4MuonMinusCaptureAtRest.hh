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
// ------------------------------------------------------------
//      GEANT 4 class header file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4MuonMinusCaptureAtRest physics process ------
//                   by Larry Felawka (TRIUMF)
//                     E-mail: felawka@alph04.triumf.ca
//                   and Art Olin (TRIUMF)
//                     E-mail: olin@triumf.ca
//                            April 1998
//-----------------------------------------------------------------------------
//
// Class Description
// Process for nuclear capture of muon- at rest takes into account Fermi model of
// muon capture in compounds, simplified EM cascade model, muon decay from K-shell,
// and muon nucleus reaction 
// Class Description - End
//

//
// Modifications: 
// 18/08/2000  V.Ivanchenko Update description, new method to simulate capture
// 17/05/2006  V.Ivanchenko Cleanup
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
#include "G4Fancy3DNucleus.hh"
#include "G4ExcitationHandler.hh"

class G4MuonMinusCaptureAtRest : public G4VRestProcess
 
{ 
public:
 
  G4MuonMinusCaptureAtRest(const G4String& processName ="muMinusCaptureAtRest", 
			   G4ProcessType   aType = fHadronic );

  ~G4MuonMinusCaptureAtRest();

  G4bool IsApplicable(const G4ParticleDefinition&);

  void BuildPhysicsTable(const G4ParticleDefinition&) 
  {};

  G4double AtRestGetPhysicalInteractionLength(const G4Track&, G4ForceCondition*)
  {return 0;};

  G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&); 

  G4double GetMeanLifeTime(const G4Track&, G4ForceCondition*) 
  {return 0;};

private:

  // hide assignment operator as private 
  G4MuonMinusCaptureAtRest& operator=(const G4MuonMinusCaptureAtRest &right);
  G4MuonMinusCaptureAtRest(const G4MuonMinusCaptureAtRest& );
   
  G4ReactionProductVector * DoMuCapture();

  G4int      nCascade;
  G4double   targetZ;
  G4double   targetA;
  G4double   targetMass;

  G4StopElementSelector*   pSelector;
  G4MuMinusCaptureCascade* pEMCascade;
  G4GHEKinematicsVector*   Cascade;
  G4Fancy3DNucleus         theN;
  G4ExcitationHandler      theHandler;

};

#endif
 





