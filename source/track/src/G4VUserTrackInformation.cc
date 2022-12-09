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
// G4VUserTrackInformation class implementation
//
// Author: Makoto Asai, 2 June 2000
// --------------------------------------------------------------------

#include "G4VUserTrackInformation.hh"

// --------------------------------------------------------------------
G4VUserTrackInformation::G4VUserTrackInformation(const G4String& infoType)
{
  pType = new G4String(infoType);
}

// --------------------------------------------------------------------
G4VUserTrackInformation::~G4VUserTrackInformation()
{
  delete pType;
}

// --------------------------------------------------------------------
G4VUserTrackInformation::
G4VUserTrackInformation(const G4VUserTrackInformation& right)
{
  if(right.pType != nullptr)
    pType = new G4String(*(right.pType));
}

// --------------------------------------------------------------------
G4VUserTrackInformation&
G4VUserTrackInformation::operator=(const G4VUserTrackInformation& right)
{
  if(this != &right)
  {
    delete pType;

    if(right.pType != nullptr)
      pType = new G4String(*(right.pType));
    else
      pType = nullptr;
  }
  return *this;
}

// --------------------------------------------------------------------
const G4String& G4VUserTrackInformation::GetType() const
{
  static const G4String NOTYPE = "NONE";
  if(pType != nullptr)
    return *pType;
  return NOTYPE;
}
