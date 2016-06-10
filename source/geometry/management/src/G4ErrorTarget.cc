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
// $Id: G4ErrorTarget.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// --------------------------------------------------------------------
//      GEANT 4 class implementation file
// --------------------------------------------------------------------

#include "G4ErrorTarget.hh"

G4ErrorTarget::G4ErrorTarget()
 : theType(G4ErrorTarget_GeomVolume) {}

G4ErrorTarget::~G4ErrorTarget() {}

G4double G4ErrorTarget::GetDistanceFromPoint( const G4ThreeVector&,
                                              const G4ThreeVector& ) const
{
  return DBL_MAX;
}

G4double G4ErrorTarget::GetDistanceFromPoint( const G4ThreeVector& ) const
{
  return DBL_MAX;
}

G4bool G4ErrorTarget::TargetReached(const G4Step*)
{
  return 0;
}
