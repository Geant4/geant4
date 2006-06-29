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
#ifndef G4HadFinalState_hh
#define G4HadFinalState_hh

#include "G4ThreeVector.hh"
#include "G4DynamicParticle.hh"
#include "globals.hh"
#include <vector>
#include "G4HadSecondary.hh"
#include "G4LorentzRotation.hh"

enum G4HadFinalStateStatus{isAlive, stopAndKill, suspend};

class G4HadFinalState
{
  public:
   G4HadFinalState();
   G4int GetNumberOfSecondaries();
   void SetEnergyChange(G4double anEnergy);
   G4double GetEnergyChange();
   void SetMomentumChange(G4ThreeVector aV);
   void SetMomentumChange(G4double x, G4double y, G4double z) ;
   G4ThreeVector GetMomentumChange();
   void AddSecondary(G4DynamicParticle * aP);
   void AddSecondary(G4HadSecondary * aP);
   void SetStatusChange(G4HadFinalStateStatus aS);
   G4HadFinalStateStatus GetStatusChange();
   void Clear();
   G4LorentzRotation & GetTrafoToLab();
   void SetTrafoToLab(G4LorentzRotation & aT);
   void SetWeightChange(G4double aW);
   G4double GetWeightChange();
   G4HadSecondary * GetSecondary(size_t i);
   void SetLocalEnergyDeposit(G4double aE);
   G4double GetLocalEnergyDeposit();
   void SecondariesAreStale();
   void ClearSecondaries();
 
  private:
   G4ThreeVector theDirection;
   G4double theEnergy;
   std::vector<G4HadSecondary *> theSecs;
   G4HadFinalStateStatus theStat;
   G4LorentzRotation theT;
   G4double theW;
   G4double theEDep;
   G4bool hasStaleSecondaries;
};

#endif
