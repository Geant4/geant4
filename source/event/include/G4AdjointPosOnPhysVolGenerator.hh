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
// G4AdjointPosOnPhysVolGenerator
//
// Class description:
//
// This class is responsible for the generation of primary adjoint particles
// on the external surface of a user selected volume.
// The particles are generated uniformly on the surface with the angular
// distribution set to a cosine law relative to normal of the surface.
// It is equivalent to the flux going in from the surface if an isotropic flux
// is considered outside. 
// It uses ray tracking technique and can be applied to all kind of convex
// volumes. Using the ray tracking technique the area of the external surface
// is also computed. The area is needed to fix the weight of the primary
// adjoint particle.  
// At the time of the development of this class, generation of points on
// volume surface and computation of surface was limited in Geant4, therefore
// the general ray tracking technique was adopted. The direct method in
// G4VSolid could be now (2009) used instead.

// Author: L. Desorgher, SpaceIT GmbH - 01.06.2006
// Contract: ESA contract 21435/08/NL/AT
// Customer: ESA/ESTEC
// --------------------------------------------------------------------
#ifndef G4AdjointPosOnPhysVolGenerator_hh
#define G4AdjointPosOnPhysVolGenerator_hh 1

#include "G4VPhysicalVolume.hh"
#include "G4AffineTransform.hh"
#include "G4ThreeVector.hh"

class G4VSolid;

class G4AdjointPosOnPhysVolGenerator 
{
//---------   
  public:
//---------   

    static  G4AdjointPosOnPhysVolGenerator* GetInstance();
   
    G4VPhysicalVolume* DefinePhysicalVolume(const G4String& aName);
    void DefinePhysicalVolume1(const G4String& aName);
    G4double ComputeAreaOfExtSurface();
    G4double ComputeAreaOfExtSurface(G4int NStat);
    G4double ComputeAreaOfExtSurface(G4double epsilon);
    G4double ComputeAreaOfExtSurface(G4VSolid* aSolid);
    G4double ComputeAreaOfExtSurface(G4VSolid* aSolid,G4int NStat);
    G4double ComputeAreaOfExtSurface(G4VSolid* aSolid,G4double epsilon);
 
    void GenerateAPositionOnTheExtSurfaceOfASolid(G4VSolid* aSolid,
                                                  G4ThreeVector& p,
                                                  G4ThreeVector& direction);
    void GenerateAPositionOnTheExtSurfaceOfTheSolid(G4ThreeVector& p,
                                                    G4ThreeVector& direction);
    void GenerateAPositionOnTheExtSurfaceOfThePhysicalVolume(G4ThreeVector& p,
                                                    G4ThreeVector& direction);
    void GenerateAPositionOnTheExtSurfaceOfThePhysicalVolume(G4ThreeVector& p,
                                                    G4ThreeVector& direction,
                                                    G4double& costh_to_normal);

    inline void SetSolid(G4VSolid* aSolid)
      { theSolid=aSolid; }
    inline G4double GetAreaOfExtSurfaceOfThePhysicalVolume()
      { return AreaOfExtSurfaceOfThePhysicalVolume; }
    inline G4double GetCosThDirComparedToNormal()
      { return CosThDirComparedToNormal; }
  
//---------   
  private:   // private methods
//---------  
    G4AdjointPosOnPhysVolGenerator() = default;
   ~G4AdjointPosOnPhysVolGenerator() = default;
    G4double ComputeAreaOfExtSurfaceStartingFromSphere(G4VSolid* aSolid,
                                                       G4int NStat);
    G4double ComputeAreaOfExtSurfaceStartingFromBox(G4VSolid* aSolid,
                                                    G4int NStat);
    void GenerateAPositionOnASolidBoundary(G4VSolid* aSolid,
                                           G4ThreeVector& p,
                                           G4ThreeVector& direction);
    G4double GenerateAPositionOnASphereBoundary(G4VSolid* aSolid,
                                                G4ThreeVector& p,
                                                G4ThreeVector& direction);
    G4double GenerateAPositionOnABoxBoundary(G4VSolid* aSolid,
                                             G4ThreeVector& p,
                                             G4ThreeVector& direction);
    void ComputeTransformationFromPhysVolToWorld();

//---------   
  private:   // attributes
//---------   

   static G4ThreadLocal G4AdjointPosOnPhysVolGenerator* theInstance;
   G4VSolid* theSolid = nullptr;
   G4VPhysicalVolume* thePhysicalVolume = nullptr;

   G4bool UseSphere{true};
   G4String ModelOfSurfaceSource{"OnSolid"};
   G4AffineTransform theTransformationFromPhysVolToWorld;
   G4double AreaOfExtSurfaceOfThePhysicalVolume{0.};
   G4double CosThDirComparedToNormal{0.};
};

#endif
