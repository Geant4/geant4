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
// ************************************************************
//-----------------------------------------------------------------------------

#ifndef G4MuonMinusCaptureAtRest_h
#define G4MuonMinusCaptureAtRest_h 1
// Class Description
// Process for nuclear capture of muon- at rest; 
// to be used in your physics list in case you need this physics.
// Class Description - End

 
#include "globals.hh"
#include "Randomize.hh" 
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
  private:
  // hide assignment operator as private 
      G4MuonMinusCaptureAtRest& operator=(const G4MuonMinusCaptureAtRest &right);
      G4MuonMinusCaptureAtRest(const G4MuonMinusCaptureAtRest& );
   
  public:
 
     G4MuonMinusCaptureAtRest(const G4String& processName ="MuonMinusCaptureAtRest");
    ~G4MuonMinusCaptureAtRest();

     G4bool IsApplicable(const G4ParticleDefinition&);
  // null physics table
     void BuildPhysicsTable(const G4ParticleDefinition&){}
     G4double AtRestGetPhysicalInteractionLength(const G4Track&,
						 G4ForceCondition*);

  //    G4double GetMeanLifeTime(const G4Track&, G4ForceCondition*);

     G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&); 

     virtual G4double GetMeanLifeTime(const G4Track&, G4ForceCondition*)
     {
       return 0;
     }
  private:

  //   void GetCaptureIsotope(const G4Track& track);
     G4double GetTargetMass(G4double, G4double);
     G4ReactionProductVector * DoMuCapture(G4double aEkin);

  private:

     G4int nCascade;
     G4double targetCharge;
     G4double targetAtomicMass;
  //     G4double tDelay;
     G4StopElementSelector*   pSelector;
     G4MuMinusCaptureCascade* pEMCascade;
     G4GHEKinematicsVector* Cascade;
     G4Fancy3DNucleus theN;
     G4ExcitationHandler theHandler;

};

#endif
 





