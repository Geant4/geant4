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


