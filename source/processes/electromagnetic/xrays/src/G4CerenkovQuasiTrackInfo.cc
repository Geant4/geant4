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
#include "G4CerenkovQuasiTrackInfo.hh"

G4Allocator<G4CerenkovQuasiTrackInfo>*& aCerenkovATIAllocator()
{
  G4ThreadLocalStatic G4Allocator<G4CerenkovQuasiTrackInfo>* _instance =
    nullptr;
  return _instance;
}

G4CerenkovQuasiTrackInfo::G4CerenkovQuasiTrackInfo(
  const G4QuasiOpticalData& aOpticalData,
  G4double aPreNumPhotons, G4double aPostNumPhotons)
  : G4VAuxiliaryTrackInformation()
  , fQuasiOpticalData(aOpticalData) 
  , fPreNumPhotons(aPreNumPhotons)
  , fPostNumPhotons(aPostNumPhotons)
{}

void G4CerenkovQuasiTrackInfo::Print() const
{
  G4cout << "Auxiliary track information for a Cerenkov step" << G4endl;
}

G4CerenkovQuasiTrackInfo* G4CerenkovQuasiTrackInfo::Cast(
  const G4VAuxiliaryTrackInformation* const aATI)
{
  G4CerenkovQuasiTrackInfo* CATI = nullptr;
  if(aATI != nullptr)
  {
    // No change will be done to the pointer and to the pointed data
    auto temp = const_cast<G4VAuxiliaryTrackInformation*>(aATI);
    CATI = dynamic_cast<G4CerenkovQuasiTrackInfo*>(temp);
  }
  return CATI;
}
