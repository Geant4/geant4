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
#ifndef EventAction_h
#define EventAction_h 1
 
// -------------------------------------------------------------
//
//
// -------------------------------------------------------------
//      GEANT4 
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- EventAction -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of  
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4UserEventAction.hh"

#include "G4Event.hh"
#include "globals.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class EventActionMessenger;
class G4UImanager;

class EventAction : public G4UserEventAction
{
public: // Without description

    EventAction();
   ~EventAction();

    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);

    void SetPrintModulo(G4int val)  {printModulo = val;};
    void SetDrawFlag(G4String val)  {drawFlag = val;};
    void AddselectedEvent(G4int val) {selectedEvents.push_back(val);
                                      nSelected++;};

  private:

    EventActionMessenger* eventMessenger;
    G4UImanager*          UI;
    G4int    nEvt;
    G4int    verbose;
    G4int    printModulo;
    G4int    nSelected;
    G4String drawFlag;
    std::vector<G4int> selectedEvents;

};

#endif


