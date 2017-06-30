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
#include "G4ChannelingTrackData.hh"
#include "G4Channeling.hh"
#include "G4SystemOfUnits.hh"

G4ChannelingTrackData::G4ChannelingTrackData()
: G4VAuxiliaryTrackInformation(),
fChannelingProcess(0),
fDBL(G4ThreeVector(DBL_MAX,DBL_MAX,DBL_MAX)),
fMomCh(fDBL),
fPosCh(fDBL),
fNuD(1.),
fElD(1.),
fEFX(0.),
fEFY(0.){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ChannelingTrackData::~G4ChannelingTrackData(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ChannelingTrackData::Print() const {
    G4cout << "Nuclei Density Ratio: " << fNuD << G4endl;
    G4cout << "Electron Density Ratio: " << fElD << G4endl;
    G4cout << "Channeling Momentum (GeV/c): " << fMomCh/CLHEP::GeV << G4endl;
    G4cout << "Channeling Position (angstrom): " << fPosCh/CLHEP::angstrom << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
