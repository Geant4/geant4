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
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************

#ifndef LISAEventAction_h
#define LISAEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ios.hh"

class LISAPrimaryGeneratorAction;
class LISASteppingAction;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class LISAEventAction : public G4UserEventAction {

  public:
    LISAEventAction(LISAPrimaryGeneratorAction*, LISASteppingAction*);
    virtual ~LISAEventAction();
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    LISAPrimaryGeneratorAction* genAction;
    LISASteppingAction* stepAction;

  public:
  inline void SetFilename(const G4String& name)  {filename = name;};


  private:
    G4int event_id;
    G4double energy_pri;
    G4int  charge_in[2];
    G4int charge_out[2];
    G4int charge_tot[2];
    const long* seeds;
    G4String filename;

};

#endif


