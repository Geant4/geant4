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
// $Id: G4ExtrudedSolid.hh 104316 2017-05-24 13:04:23Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4ExtrudedSolid
//
// Class description:
//
// G4ExtrudedSolid is a solid which represents the extrusion of an arbitrary
// polygon with fixed outline in the defined Z sections.
// The z-sides of the solid are the scaled versions of the same polygon.
// The solid is implemented as a specification of G4TessellatedSolid.
//  
// Parameters in the constructor:
// const G4String& pName             - solid name
// std::vector<G4TwoVector> polygon  - the vertices of the outlined polygon
//                                     defined in clockwise or anti-clockwise order     
// std::vector<ZSection>             - the z-sections defined by
//                                     z position, offset and scale
//                                     in increasing z-position order
//
// Parameters in the special constructor (for solid with 2 z-sections:
// G4double halfZ                    - the solid half length in Z
// G4TwoVector off1                  - offset of the side in -halfZ
// G4double scale1                   - scale of the side in -halfZ
// G4TwoVector off2                  - offset of the side in +halfZ
// G4double scale2                   - scale of the side in -halfZ

// Author:
//   Ivana Hrivnacova, IPN Orsay
//
// --------------------------------------------------------------------

#ifndef G4ExtrudedSolid_HH
#define G4ExtrudedSolid_HH

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UEXTRUDEDSOLID 1
#endif

#if defined(G4GEOM_USE_UEXTRUDEDSOLID)
  #define G4UExtrudedSolid G4ExtrudedSolid
  #include "G4UExtrudedSolid.hh"
#else

#include <vector>

#include "G4TwoVector.hh"

#include "G4TessellatedSolid.hh"

class G4VFacet;

class G4ExtrudedSolid : public G4TessellatedSolid
{

  public:  // without description

    struct ZSection
    {
      ZSection(G4double z, const G4TwoVector& offset, G4double scale)
        : fZ(z), fOffset(offset), fScale(scale) {}

      G4double    fZ;
      G4TwoVector fOffset;
      G4double    fScale;
    };

  public:  // with description

     G4ExtrudedSolid( const G4String&                 pName,
                      const std::vector<G4TwoVector>& polygon,
                      const std::vector<ZSection>&    zsections);
       // General constructor

     G4ExtrudedSolid( const G4String&                 pName,
                      const std::vector<G4TwoVector>& polygon,
                            G4double                  halfZ,
                      const G4TwoVector& off1, G4double scale1,
                      const G4TwoVector& off2, G4double scale2 );
       // Special constructor for solid with 2 z-sections

     virtual ~G4ExtrudedSolid();
       // Destructor

    // Accessors

    inline G4int       GetNofVertices() const;
    inline G4TwoVector GetVertex(G4int index) const;
    inline std::vector<G4TwoVector> GetPolygon() const;

    inline G4int       GetNofZSections() const;
    inline ZSection    GetZSection(G4int index) const;
    inline std::vector<ZSection> GetZSections() const;

    // Solid methods                                

    EInside  Inside (const G4ThreeVector &p) const;
    G4double DistanceToOut(const G4ThreeVector &p,
                           const G4ThreeVector &v,
                           const G4bool calcNorm=false,
                                 G4bool *validNorm=0, G4ThreeVector *n=0) const;
    G4double DistanceToOut (const G4ThreeVector &p) const;
    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const;
    G4GeometryType GetEntityType () const;
    G4VSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream &os) const;

  public:  // without description

    G4ExtrudedSolid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4ExtrudedSolid(const G4ExtrudedSolid& rhs);
    G4ExtrudedSolid& operator=(const G4ExtrudedSolid& rhs);
      // Copy constructor and assignment operator.

  private:

    void ComputeProjectionParameters();
    
    G4ThreeVector GetVertex(G4int iz, G4int ind) const;
    G4TwoVector ProjectPoint(const G4ThreeVector& point) const;

    G4bool IsSameLine(const G4TwoVector& p,
                      const G4TwoVector& l1,
                      const G4TwoVector& l2) const;
    G4bool IsSameLineSegment(const G4TwoVector& p,
                             const G4TwoVector& l1,
                             const G4TwoVector& l2) const;
    G4bool IsSameSide(const G4TwoVector& p1,
                      const G4TwoVector& p2,
                      const G4TwoVector& l1,
                      const G4TwoVector& l2) const;
    G4bool IsPointInside(const G4TwoVector& a,
                         const G4TwoVector& b,
                         const G4TwoVector& c,
                         const G4TwoVector& p) const;
    G4double GetAngle(const G4TwoVector& p0,
                      const G4TwoVector& pa,
                      const G4TwoVector& pb) const;                      
      
    G4VFacet* MakeDownFacet(G4int ind1, G4int ind2, G4int ind3) const;      
    G4VFacet* MakeUpFacet(G4int ind1, G4int ind2, G4int ind3) const;      

    G4bool AddGeneralPolygonFacets();
    G4bool MakeFacets();

  private:

    G4int       fNv;
    G4int       fNz;
    std::vector<G4TwoVector> fPolygon;
    std::vector<ZSection>    fZSections;
    std::vector< std::vector<G4int> > fTriangles;
    G4bool          fIsConvex;
    G4GeometryType  fGeometryType;

    std::vector<G4double>      fKScales;
    std::vector<G4double>      fScale0s;
    std::vector<G4TwoVector>   fKOffsets;
    std::vector<G4TwoVector>   fOffset0s;
};    

#include "G4ExtrudedSolid.icc"

#endif

#endif
