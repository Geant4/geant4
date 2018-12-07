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
/// \file optical/wls/include/WLSSteppingAction.hh
/// \brief Definition of the WLSSteppingAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSSteppingAction_h
#define WLSSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class WLSDetectorConstruction;
class WLSSteppingActionMessenger;

class G4Track;
class G4StepPoint;

class G4OpBoundaryProcess;

class WLSSteppingAction : public G4UserSteppingAction
{
  public:

    WLSSteppingAction(WLSDetectorConstruction*);
    virtual ~WLSSteppingAction();

    virtual void UserSteppingAction(const G4Step*);
 
    // Set the bounce limit, 0 for no limit
    void  SetBounceLimit(G4int);
 
    G4int GetNumberOfBounces();
    G4int GetNumberOfClad1Bounces();
    G4int GetNumberOfClad2Bounces();
    G4int GetNumberOfWLSBounces();
    // return number of successful events and reset the counter
    G4int ResetSuccessCounter();
 
  private:

    // Artificially kill the photon after it has bounced more than this number
    G4int fBounceLimit;
    // number of photons that reach the end
    G4int fCounterEnd;
    // number of photons that didn't make it to the end
    G4int fCounterMid;
    // total number of bounces that a photon been through
    G4int fCounterBounce;
    // number of bounces that a photon been through within the fibre
    G4int fCounterWLSBounce;
    // number of bounces that a photon been through from Cladding 1 to 2
    G4int fCounterClad1Bounce;
    // number of bounces that a photon been through from Cladding 2 to World
    G4int fCounterClad2Bounce;

    // initial gamma of the photon
    G4double fInitGamma;
    // initial theta of the photon
    G4double fInitTheta;

    G4OpBoundaryProcess* fOpProcess;

    // maximum number of save states
    static G4int fMaxRndmSave;
 
    WLSDetectorConstruction* fDetector;

    WLSSteppingActionMessenger* fSteppingMessenger;

    inline void ResetCounters()
    { 
      fCounterBounce = fCounterWLSBounce =
      fCounterClad1Bounce = fCounterClad2Bounce = 0;
      fInitGamma = fInitTheta = -1;
    }

    // save the random status into a sub-directory
    // Pre: subDir must be empty or ended with "/"
    inline void SaveRandomStatus(G4String subDir);

};

#endif
