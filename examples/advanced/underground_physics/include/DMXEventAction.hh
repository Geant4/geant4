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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// EventAction header
// --------------------------------------------------------------

#ifndef DMXEventAction_h
#define DMXEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ios.hh"

#include "DMXScintHit.hh"
#include "DMXPmtHit.hh"

class DMXEventActionMessenger;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class DMXEventAction : public G4UserEventAction {

  public:
    DMXEventAction();
    virtual ~DMXEventAction();
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    void writeScintHitsToFile(void);
    void writePmtHitsToFile(const DMXPmtHitsCollection*);
    void drawTracks(const G4Event*);

  public:
    void SetDrawTrksFlag (G4String val)     {drawTrksFlag    = val;};
    G4String GetDrawTrksFlag()              {return drawTrksFlag;};

    void SetDrawColsFlag (G4String val)     {drawColsFlag    = val;};
    G4String GetDrawColsFlag()              {return drawColsFlag;};

    void SetDrawHitsFlag (G4int val)        {drawHitsFlag    = val;};
    void SetSavePmtFlag  (G4int val)        {savePmtFlag     = val;};
    void SetSaveHitsFlag (G4int val)        {saveHitsFlag    = val;};
    void SetPrintModulo  (G4int val)        {printModulo     = val;};

  private:
    G4int event_id;

    // hits collections
    G4int scintillatorCollID;                
    G4int pmtCollID;
    G4int S_hits;
    G4int P_hits;

    // event summary
    G4double aveTimePmtHits;
    G4double totEnergy;
    G4double totEnergyGammas;
    G4double totEnergyNeutrons;
    G4double hitEnergy;
    G4double firstLXeHitTime;
    G4double particleEnergy;
    G4String particleName;
    G4String firstParticleName;

    G4bool gamma_ev;
    G4bool neutron_ev;
    G4bool positron_ev;
    G4bool electron_ev;
    G4bool proton_ev;
    G4bool other_ev;
    G4bool start_gamma;
    G4bool start_neutron;

    // messenger
    G4String drawTrksFlag;
    G4String drawColsFlag;
    G4int drawHitsFlag;         
    G4int savePmtFlag;         
    G4int saveHitsFlag;         
    G4int printModulo;                         
    DMXEventActionMessenger*  eventMessenger;

};

#endif

