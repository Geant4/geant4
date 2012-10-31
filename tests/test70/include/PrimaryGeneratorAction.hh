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
// -------------------------------------------------------------------
// $Id: PrimaryGeneratorAction.hh,v 1.1 2008/06/04 12:58:24 sincerti Exp $
// -------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction();

    void GeneratePrimaries(G4Event*);
    inline void SetIncidentEnergy(G4double e);
    inline G4double GetIncidentEnergy();

    inline G4ParticleDefinition* GetParticleDefinition();

private:

    G4double e0;
    G4ParticleGun*           particleGun;
};

inline void PrimaryGeneratorAction :: SetIncidentEnergy (G4double e)
{
    e0 = e*CLHEP::keV;
}

inline G4double PrimaryGeneratorAction :: GetIncidentEnergy()
{
    return e0 ;
}

inline  G4ParticleDefinition* PrimaryGeneratorAction::GetParticleDefinition()
{
    return particleGun->GetParticleDefinition();
}
#endif
