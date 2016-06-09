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
// $Id: G4AuxiliaryNavServices.hh,v 1.3 2006/06/29 18:35:32 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 
// class G4NormalNavigation
//
// Class description:
//
// Utility for navigation in volumes containing only G4PVPlacement
// daughter volumes.

// History:
// - Created: Paul Kent, Aug 96
// --------------------------------------------------------------------
#ifndef G4AuxiliaryNavServices_hh
#define G4AuxiliaryNavServices_hh

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"
#include "G4AffineTransform.hh"

class G4AuxiliaryNavServices
{

 public:  // with description

   static G4bool CheckPointOnSurface( const G4VSolid* sampleSolid, 
                                      const G4ThreeVector& localPoint, 
                                      const G4ThreeVector* globalDirection, 
                                      const G4AffineTransform& sampleTransform,
                                      const G4bool locatedOnEdge);
     //
     // Is the track (point, direction) inside the solid 'sampleSolid' ? 
     // Returns true if we are going to enter the volume,
     // which is the case if:
     //   - the point is inside
     //   - the point is on the surface and the direction points inside
     //     or along it.
     // Else returns false.

 private:
 
   G4bool testOne();

};

#include "G4AuxiliaryNavServices.icc"

#endif
