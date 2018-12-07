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
/// \file optical/wls/include/WLSPrimaryGeneratorAction.hh
/// \brief Definition of the WLSPrimaryGeneratorAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSPrimaryGeneratorAction_h
#define WLSPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4AffineTransform.hh"

class G4GeneralParticleSource;

class G4Event;
class G4PhysicsTable;

class WLSDetectorConstruction;
class WLSPrimaryGeneratorMessenger;

class WLSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

    WLSPrimaryGeneratorAction(WLSDetectorConstruction*);
    virtual ~WLSPrimaryGeneratorAction();

  public:

    virtual void GeneratePrimaries(G4Event*);

    void BuildEmissionSpectrum();

    void SetOptPhotonPolar(G4double);

    void SetDecayTimeConstant(G4double);

  protected:

    G4PhysicsTable* fIntegralTable;

  private:

    void SetOptPhotonPolar();
    void SetOptPhotonTime();
 
    WLSDetectorConstruction*   fDetector;
    G4GeneralParticleSource*   fParticleGun;
    WLSPrimaryGeneratorMessenger* fGunMessenger;

    static G4bool fFirst;

    G4double fTimeConstant;

};

#endif
