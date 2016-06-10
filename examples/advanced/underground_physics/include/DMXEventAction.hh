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
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
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

class DMXRunAction;
class DMXPrimaryGeneratorAction;
class DMXEventActionMessenger;
class DMXAnalysisManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class DMXEventAction : public G4UserEventAction {

  public:
    DMXEventAction();
    virtual ~DMXEventAction();
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    void writeScintHitsToFile();
    void writePmtHitsToFile(const DMXPmtHitsCollection*);
    void drawTracks(const G4Event*);

  public:
    void SetDrawTrksFlag (G4String val)     {drawTrksFlag    = val;};
    G4String GetDrawTrksFlag() const         {return drawTrksFlag;};

    void SetDrawColsFlag (G4String val)     {drawColsFlag    = val;};
    G4String GetDrawColsFlag() const  {return drawColsFlag;};

    void SetDrawHitsFlag (G4int val)        {drawHitsFlag    = val;};
    void SetSavePmtFlag  (G4int val)        {savePmtFlag     = val;};
    void SetSaveHitsFlag (G4int val)        {saveHitsFlag    = val;};
    void SetPrintModulo  (G4int val)        {printModulo     = val;};

  private:
    G4int event_id;

    const long* seeds;
    G4double energy_pri;

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
    G4double firstParticleE;
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

    const DMXRunAction*    runAct;  //pointer to run action
    const DMXPrimaryGeneratorAction* genAction; // pointer to particle generator
  std::ofstream *hitsfile;
  std::ofstream *pmtfile;
};

#endif

