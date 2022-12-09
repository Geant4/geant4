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
//      ------------ GammaRayTelDigi ------
//           by F.Longo, R.Giannitrapani & G.Santin (21 oct 2001)
//
// ************************************************************

#include "GammaRayTelDigi.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreadLocal G4Allocator<GammaRayTelDigi> *digitAllocator{nullptr};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDigi::GammaRayTelDigi() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDigi::~GammaRayTelDigi() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDigi::GammaRayTelDigi(const GammaRayTelDigi &right) :
    planeType(right.planeType),
    planeNumber(right.planeNumber),
    stripNumber(right.stripNumber),
    digitType(right.digitType),
    energy(right.energy)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto GammaRayTelDigi::operator=(const GammaRayTelDigi &right) -> const GammaRayTelDigi& {
    planeType = right.planeType;
    planeNumber = right.planeNumber;
    stripNumber = right.stripNumber;
    digitType = right.digitType;
    energy = right.energy;
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto GammaRayTelDigi::operator==(const GammaRayTelDigi &right) const -> G4bool {
    return (
        (planeType == right.planeType) &&
        (planeNumber == right.planeNumber) &&
        (stripNumber == right.stripNumber) &&
        (digitType == right.digitType) &&
        (energy == right.energy)
    );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDigi::Draw() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDigi::Print() {
}
