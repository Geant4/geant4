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
#include <strstream>

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


