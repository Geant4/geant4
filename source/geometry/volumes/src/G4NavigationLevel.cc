//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4NavigationLevel.cc,v 1.3 2005/06/06 13:41:26 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 30.09.97 J.Apostolakis Initial version. 
//                    
// ----------------------------------------------------------------------

#include "G4NavigationLevel.hh"

G4Allocator<G4NavigationLevel> aNavigationLevelAllocator;

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
  if( fLevelRep->RemoveAReference() )
    delete fLevelRep; 
}

G4NavigationLevel& G4NavigationLevel::operator=(const G4NavigationLevel &right)
{ 
  if ( &right != this )
  {
    right.fLevelRep->AddAReference(); 
    if( fLevelRep->RemoveAReference() )
      delete fLevelRep; 
    fLevelRep = right.fLevelRep;
  }
  return *this;
} 
