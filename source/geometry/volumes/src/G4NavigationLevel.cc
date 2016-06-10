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
// $Id: G4NavigationLevel.cc 86527 2014-11-13 15:06:24Z gcosmo $
//
// 30.09.97 J.Apostolakis Initial version. 
//                    
// ----------------------------------------------------------------------

#include "G4NavigationLevel.hh"

G4ThreadLocal G4Allocator<G4NavigationLevel> *aNavigationLevelAllocator = 0;

G4NavigationLevel::G4NavigationLevel( G4VPhysicalVolume* pPhysVol,
                                const G4AffineTransform& afTransform,
                                      EVolume            volTp,
                                      G4int              repNo )
{
  fLevelRep = new G4NavigationLevelRep( pPhysVol, afTransform, volTp, repNo );
}

G4NavigationLevel::G4NavigationLevel( G4VPhysicalVolume* pPhysVol,
                                const G4AffineTransform& levelAbove,
                                const G4AffineTransform& relativeCurrent,
                                      EVolume            volTp,
                                      G4int              repNo )
{
  fLevelRep = new G4NavigationLevelRep( pPhysVol, 
                                        levelAbove, 
                                        relativeCurrent, 
                                        volTp, 
                                        repNo );
}

G4NavigationLevel::G4NavigationLevel()
{
  fLevelRep = new G4NavigationLevelRep();
}

G4NavigationLevel::G4NavigationLevel(const G4NavigationLevel& right)
  : fLevelRep( right.fLevelRep )
{
  fLevelRep->AddAReference(); 
}

G4NavigationLevel::~G4NavigationLevel()
{
  if( fLevelRep->RemoveAReference() )  { delete fLevelRep; }
}

G4NavigationLevel& G4NavigationLevel::operator=(const G4NavigationLevel &right)
{
  if ( &right != this )
  {
    right.fLevelRep->AddAReference(); 
    if( fLevelRep->RemoveAReference() )  { delete fLevelRep; }
    fLevelRep = right.fLevelRep;
  }
  return *this;
} 
