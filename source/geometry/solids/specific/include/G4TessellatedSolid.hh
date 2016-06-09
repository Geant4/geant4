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
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * subject DEFCON 705 IPR conditions.                               *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4TessellatedSolid.hh,v 1.3 2006/06/29 18:47:35 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4TessellatedSolid.hh
//
// Date:                15/06/2005
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            UK Ministry of Defence : RAO CRP TD Electronic Systems
// Contract:            C/MAT/N03517
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
// 22 November 2005, F Lei, 
//  - Added GetPolyhedron()
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK
//  - Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class description:
//
//    G4TessellatedSolid is a special Geant4 solid defined by a number of 
//    facets (G4VFacet). It is important that the supplied facets shall form a
//    fully enclose space which is the solid. 
//    At the moment only two types of facet can be used for the construction of 
//    a G4TessellatedSolid, i.e. the G4TriangularFacet and G4QuadrangularFacet.
//
//    How to contruct a G4TessellatedSolid:
//  
//    First declare a tessellated solid:
//
//      G4TessellatedSolid* solidTarget = new G4TessellatedSolid("Solid_name");
//
//    Define the facets which form the solid
// 
//      G4double targetSiz = 10*cm ;
//      G4TriangularFacet *facet1 = new
//      G4TriangularFacet (G4ThreeVector(-targetSize,-targetSize,        0.0),
//                         G4ThreeVector(+targetSize,-targetSize,        0.0),
//                         G4ThreeVector(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      G4TriangularFacet *facet2 = new
//      G4TriangularFacet (G4ThreeVector(+targetSize,-targetSize,        0.0),
//                         G4ThreeVector(+targetSize,+targetSize,        0.0),
//                         G4ThreeVector(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      G4TriangularFacet *facet3 = new
//      G4TriangularFacet (G4ThreeVector(+targetSize,+targetSize,        0.0),
//                         G4ThreeVector(-targetSize,+targetSize,        0.0),
//                         G4ThreeVector(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      G4TriangularFacet *facet4 = new
//      G4TriangularFacet (G4ThreeVector(-targetSize,+targetSize,        0.0),
//                         G4ThreeVector(-targetSize,-targetSize,        0.0),
//                         G4ThreeVector(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      G4QuadrangularFacet *facet5 = new
//      G4QuadrangularFacet (G4ThreeVector(-targetSize,-targetSize,      0.0),
//                           G4ThreeVector(-targetSize,+targetSize,      0.0),
//                           G4ThreeVector(+targetSize,+targetSize,      0.0),
//                           G4ThreeVector(+targetSize,-targetSize,      0.0),
//                           ABSOLUTE);
//
//    Then add the facets to the solid:    
//
//      solidTarget->AddFacet((G4VFacet*) facet1);
//      solidTarget->AddFacet((G4VFacet*) facet2);
//      solidTarget->AddFacet((G4VFacet*) facet3);
//      solidTarget->AddFacet((G4VFacet*) facet4);
//      solidTarget->AddFacet((G4VFacet*) facet5);
//
//    Finally declare the solid is complete:
//
//      solidTarget->SetSolidClosed(true);
//
///////////////////////////////////////////////////////////////////////////////
#ifndef G4TessellatedSolid_hh
#define G4TessellatedSolid_hh 1

#include "G4VSolid.hh"
#include "G4VFacet.hh"
#include "G4VGraphicsScene.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4VisExtent.hh"
#include "globals.hh"

#include <iostream>
#include <vector>
#include <map>

class G4TessellatedSolid : public G4VSolid
{
  public:  // with description

    G4TessellatedSolid ();
    G4TessellatedSolid (const G4String &name);
    ~G4TessellatedSolid ();
    
    G4TessellatedSolid (const G4TessellatedSolid &s);
    const G4TessellatedSolid &operator= (const G4TessellatedSolid &s);
    const G4TessellatedSolid &operator+= (const G4TessellatedSolid &right);
    
    G4bool AddFacet (G4VFacet *aFacet);
    G4VFacet *GetFacet (size_t i) const;
    size_t GetNumberOfFacets () const;
    
//  G4double GetCubicVolume ();
//
//  void ComputeDimensions (G4VPVParameterisation* p, const G4int n,
//                          const G4VPhysicalVolume* pRep) const;
    
    EInside Inside (const G4ThreeVector &p) const;
    G4ThreeVector SurfaceNormal (const G4ThreeVector &p) const;
    G4double DistanceToIn(const G4ThreeVector &p, const G4ThreeVector &v) const;
    G4double DistanceToIn(const G4ThreeVector &p) const;
    G4double DistanceToOut(const G4ThreeVector &p, const G4ThreeVector &v,
                           const G4bool calcNorm=false,
                           G4bool *validNorm=0, G4ThreeVector *n=0) const;
    G4double DistanceToOut (const G4ThreeVector &p) const;
    G4GeometryType GetEntityType () const;
    
    void SetSolidClosed (const G4bool t);
    G4bool GetSolidClosed () const;
        
    G4bool CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform, G4double& pMin,
                           G4double& pMax) const;

    std::ostream &StreamInfo(std::ostream &os) const;

  // Functions for visualization
 
    void          DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4VisExtent   GetExtent () const;
    G4double      GetMinXExtent () const;
    G4double      GetMaxXExtent () const;
    G4double      GetMinYExtent () const;
    G4double      GetMaxYExtent () const;
    G4double      GetMinZExtent () const;
    G4double      GetMaxZExtent () const;
    G4Polyhedron* CreatePolyhedron () const;
    G4Polyhedron* GetPolyhedron      () const;
    G4NURBS*      CreateNURBS () const;
 
  protected:  // with description
 
    void DeleteObjects ();
    void CopyObjects (const G4TessellatedSolid &s);
    G4ThreeVectorList*
    CreateRotatedVertices(const G4AffineTransform& pTransform) const;
      // Create the List of transformed vertices in the format required
      // for G4VSolid:: ClipCrossSection and ClipBetweenSections.

  private: 

    mutable G4Polyhedron* fpPolyhedron;

    std::vector<G4VFacet *>  facets;
    G4GeometryType           geometryType;
    G4double                 cubicVolume;
    G4ThreeVectorList        vertexList;
    G4double                 xMinExtent;
    G4double                 xMaxExtent;
    G4double                 yMinExtent;
    G4double                 yMaxExtent;
    G4double                 zMinExtent;
    G4double                 zMaxExtent;
    G4bool                   solidClosed;
    
    G4double                 dirTolerance;
};

#endif
