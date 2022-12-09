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
/// \file electromagnetic/TestEm5/include/EventAction.hh
/// \brief Definition of the EventAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4VProcess;

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction() = default;

  void BeginOfEventAction(const G4Event*) override;
  void   EndOfEventAction(const G4Event*) override;
    
  void AddEnergy      (G4double edep)   {fEnergyDeposit  += edep;};
  void AddTrakLenCharg(G4double length) {fTrakLenCharged += length;};
  void AddTrakLenNeutr(G4double length) {fTrakLenNeutral += length;};
    
  void CountStepsCharg () {++fNbStepsCharged;};
  void CountStepsNeutr (const G4VProcess*);
  
  void SetTransmitFlag (G4int flag) 
  {if (flag > fTransmitFlag) fTransmitFlag = flag;};
  void SetReflectFlag  (G4int flag) 
  {if (flag > fReflectFlag)   fReflectFlag = flag;};
                                             
        
private:
  G4double fEnergyDeposit;
  G4double fTrakLenCharged, fTrakLenNeutral;
  G4int fNbStepsCharged, fNbStepsNeutral;
  G4int fTransmitFlag, fReflectFlag;
  G4int fTypes[4];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
