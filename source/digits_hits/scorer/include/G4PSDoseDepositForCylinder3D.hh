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
//

#ifndef G4PSDoseDepositForCylinder3D_h
#define G4PSDoseDepositForCylinder3D_h 1

#include "G4PSDoseDeposit3D.hh"

////////////////////////////////////////////////////////////////////////////////
// (Description)
//   Primitive scorer class for scoring dose deposit in the geometry volume.
//
// Created: 2008-08-14  Tsukasa ASO
// 2010-07-22   Introduce Unit specification.
// 2011-03-24   Give Size and Segmentation for relicated volume in cylinder. 
// 
///////////////////////////////////////////////////////////////////////////////

class G4PSDoseDepositForCylinder3D : public G4PSDoseDeposit3D
{
 public: // with description
      G4PSDoseDepositForCylinder3D(G4String name, 
			G4int ni=1,G4int nj=1, G4int nk=1,
			G4int depi=2, G4int depj=1, G4int depk=0);
      G4PSDoseDepositForCylinder3D(G4String name, const G4String& unit,
			G4int ni=1,G4int nj=1, G4int nk=1,
			G4int depi=2, G4int depj=1, G4int depk=0);

      virtual ~G4PSDoseDepositForCylinder3D();

      void SetCylinderSize(G4double dr, G4double dz);
      void SetNumberOfSegments(G4int nSeg[3]);

  protected: // with description
      virtual G4double ComputeVolume(G4Step*, G4int idx);

  private:
  // Order of segmentation (Z PHI R) in CylinderMesh
   G4ThreeVector cylinderSize;  
   G4int nSegment[3];

};
#endif
