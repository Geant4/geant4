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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelTrackerHit ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

#include "GammaRayTelTrackerHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreadLocal G4Allocator<GammaRayTelTrackerHit> *trackerHitAllocator{nullptr};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelTrackerHit::GammaRayTelTrackerHit() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelTrackerHit::~GammaRayTelTrackerHit() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelTrackerHit::GammaRayTelTrackerHit(const GammaRayTelTrackerHit &right) {
    stripNumber = right.stripNumber;
    siliconPlaneNumber = right.siliconPlaneNumber;
    isXPlane = right.isXPlane;
    siliconDepositedEnergy = right.siliconDepositedEnergy;
    position = right.position;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto GammaRayTelTrackerHit::operator=(const GammaRayTelTrackerHit &right) -> const GammaRayTelTrackerHit& {
    stripNumber = right.stripNumber;
    siliconPlaneNumber = right.siliconPlaneNumber;
    isXPlane = right.isXPlane;
    siliconDepositedEnergy = right.siliconDepositedEnergy;
    position = right.position;

    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto GammaRayTelTrackerHit::operator==(const GammaRayTelTrackerHit &right) const -> G4bool {
    return (
        (stripNumber == right.stripNumber) &&
        (isXPlane == right.isXPlane) &&
        (siliconDepositedEnergy == right.siliconDepositedEnergy) &&
        (position == right.position)
    );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelTrackerHit::Draw() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelTrackerHit::Print() {
}
