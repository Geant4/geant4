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
#include "G4HadFinalState.hh"
#include "G4HadronicException.hh"

   G4HadFinalState::G4HadFinalState()
   : theDirection(0,0,1), theEnergy(-1), theStat(isAlive), 
     theW(1.), theEDep(0.), hasStaleSecondaries(false){}

   G4int G4HadFinalState::GetNumberOfSecondaries() {return theSecs.size();}

   void G4HadFinalState::SetEnergyChange(G4double anEnergy) 
   {
     theEnergy=anEnergy;
     if(theEnergy<0) 
     {
       std::cout << "Final state energy was: E = "<<theEnergy<<G4endl;
       throw G4HadronicException(__FILE__, __LINE__, "G4HadFinalState: fatal - negative energy");
     }
   }

   G4double G4HadFinalState::GetEnergyChange() {return theEnergy;}

   void G4HadFinalState::SetMomentumChange(G4ThreeVector aV) {theDirection=aV;}

   void G4HadFinalState::SetMomentumChange(G4double x, G4double y, G4double z) 
   {
     theDirection = G4ThreeVector(x,y,z);
     if(std::fabs(theDirection.mag()-1)>0.001) 
     {
       G4cout <<"We have negative theDirection.mag() = "<<theDirection.mag()<<G4endl;
       throw G4HadronicException(__FILE__, __LINE__, "G4HadFinalState: fatal - negative direction.mag().");
     }
   }

   G4ThreeVector G4HadFinalState::GetMomentumChange() {return theDirection;}

   void G4HadFinalState::AddSecondary(G4DynamicParticle * aP) {theSecs.push_back(new G4HadSecondary(aP));}

   void G4HadFinalState::AddSecondary(G4HadSecondary * aP) {theSecs.push_back(aP);}

   void G4HadFinalState::SetStatusChange(G4HadFinalStateStatus aS){theStat=aS;}

   G4HadFinalStateStatus G4HadFinalState::GetStatusChange(){return theStat;}

   void G4HadFinalState::ClearSecondaries()
   {
     theSecs.clear();
   }
   
   void G4HadFinalState::Clear()
   {
     theDirection = G4ThreeVector(0,0,1);
     theEnergy = -1;
     theStat = isAlive;
     theW = 1.;
     theEDep = 0.;
     if(!hasStaleSecondaries) 
     {
       for(size_t i=0; i<theSecs.size(); i++) delete theSecs[i];
     }
     theSecs.clear();
     hasStaleSecondaries = false;
   }

   void G4HadFinalState::SecondariesAreStale() {hasStaleSecondaries = true;}

   G4LorentzRotation & G4HadFinalState::GetTrafoToLab() {return theT;}

   void G4HadFinalState::SetTrafoToLab(G4LorentzRotation & aT) {theT = aT;}

   void G4HadFinalState::SetWeightChange(G4double aW){ theW=aW;}

   G4double G4HadFinalState::GetWeightChange() {return theW;}

   G4HadSecondary * G4HadFinalState::GetSecondary(size_t i) 
   {
     if(i>theSecs.size())
     {
       throw G4HadronicException(__FILE__, __LINE__, 
            "Trying direct access to secondary beyond end of list");
     }
     return theSecs[i];
   }

   void G4HadFinalState::SetLocalEnergyDeposit(G4double aE) {theEDep=aE;}

   G4double G4HadFinalState::GetLocalEnergyDeposit() {return theEDep;}
