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
// $Id: EventAction.hh,v 1.2 2004/06/21 10:57:10 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;
class HistoManager;
class EventMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* , HistoManager*);
   ~EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    void AddEnergy      (G4double edep)   {EnergyDeposit  += edep;};
    void AddTrakLenCharg(G4double length) {TrakLenCharged += length;};
    void AddTrakLenNeutr(G4double length) {TrakLenNeutral += length;};
    
    void CountStepsCharg ()               {nbStepsCharged++ ;};
    void CountStepsNeutr ()               {nbStepsNeutral++ ;};
    
    void SetTransmitFlag (G4int flag) 
                           {if (flag > TransmitFlag) TransmitFlag = flag;};
    void SetReflectFlag  (G4int flag) 
                           {if (flag > ReflectFlag)   ReflectFlag = flag;};
			           	  
    void SetDrawFlag(G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int val)  {printModulo = val;};
        
  private:
    RunAction*    runaction;
    HistoManager* histoManager;
    
    G4double EnergyDeposit;
    G4double TrakLenCharged, TrakLenNeutral;
    G4int    nbStepsCharged, nbStepsNeutral;
    G4int    TransmitFlag,   ReflectFlag; 
    
    G4String drawFlag;
    G4int    printModulo;
    
    EventMessenger* eventMessenger;                    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
