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
// $Id: G4VUserTrackInformation.cc 68795 2013-04-05 13:24:46Z gcosmo $
//
// 
// --------------------------------------------------------------

#include "G4VUserTrackInformation.hh"

G4VUserTrackInformation::G4VUserTrackInformation()
  : pType(0)
{;}

G4VUserTrackInformation::G4VUserTrackInformation(const G4String& infoType)
{ 
  pType = new G4String(infoType) ;
}

G4VUserTrackInformation::~G4VUserTrackInformation()
{
  if (pType!=0) delete pType;
}

G4VUserTrackInformation::G4VUserTrackInformation(const  G4VUserTrackInformation& right) 
  :pType(0)
{
  if (right.pType!=0)  pType = new G4String(*(right.pType));
}


G4VUserTrackInformation& G4VUserTrackInformation::operator=(const G4VUserTrackInformation& right)
{
  if (this != &right) {
    if (pType !=0) delete pType;
    if (right.pType) pType = new G4String(*(right.pType));
    else pType=0;
  }
  return *this;
}
  

const G4String& G4VUserTrackInformation::GetType() const
{
  static const G4String NOTYPE="NONE";
  if(pType!=0) return *pType;
  else return NOTYPE;
}

