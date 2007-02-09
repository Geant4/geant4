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
// $Id: G4ExtrudedSolid.hh,v 1.1 2007-02-09 12:05:50 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// polygon with fixed outline in the two Z sections.
// The z-sides of the solid are the scaled versions of the same polygon.
// The solid is implemented as a specification of G4TessellatedSolid.
//  
// Parameters in the constructor:
// const G4String& pName             - solid name
// std::vector<G4TwoVector> polygon  - the vertices of the outlined polygon
//                                     defined in clock-wise order     
// G4double hz                       - the solid half length in Z
// G4TwoVector off1                  - offset of the side in -hz
// G4double scale1                   - scale of the side in -hz
// G4TwoVector off2                  - offset of the side in +hz
// G4double scale2                   - scale of the side in -hz

// Author:
//   Ivana Hrivnacova, IPN Orsay
// --------------------------------------------------------------------

#ifndef G4ExtrudedSolid_HH
#define G4ExtrudedSolid_HH

#include <vector>

#include "G4TwoVector.hh"

#include "G4TessellatedSolid.hh"

class G4ExtrudedSolid : public G4TessellatedSolid
{

  public:  // with description

     G4ExtrudedSolid( const G4String&                pName,
                            std::vector<G4TwoVector> polygon,       
                            G4double                 hz,
                            G4TwoVector off1, G4double scale1,
                            G4TwoVector off2, G4double scale2 );
       // Constructor

     virtual ~G4ExtrudedSolid();
       // Destructor

    // Accessors

    inline G4double    GetZHalfLength()  const;
    inline G4TwoVector GetOffset1() const;
    inline G4double    GetScale1() const;
    inline G4TwoVector GetOffset2() const;
    inline G4double    GetScale2() const;
    inline G4int       GetNofVertices() const;
    inline G4TwoVector GetVertex(G4int index) const;
    inline std::vector<G4TwoVector> GetPolygon() const;

    // Solid methods                                

    virtual G4GeometryType GetEntityType () const;
    virtual EInside  Inside (const G4ThreeVector &p) const;
    virtual G4double DistanceToOut(const G4ThreeVector &p,
                                   const G4ThreeVector &v,
                                   const G4bool calcNorm,
                                         G4bool *validNorm,
                                         G4ThreeVector *n) const;
    virtual G4double DistanceToOut (const G4ThreeVector &p) const;
    std::ostream& StreamInfo(std::ostream &os) const;

  public:  // without description

    G4ExtrudedSolid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

  private:

    G4ThreeVector GetDownVertex(G4int ind) const;
    G4ThreeVector GetUpVertex(G4int ind) const;

    G4bool IsSameLine(G4TwoVector p,
                      G4TwoVector l1, G4TwoVector l2) const;
    G4bool IsSameSide(G4TwoVector p1, G4TwoVector p2, 
                      G4TwoVector l1, G4TwoVector l2) const;
    G4bool IsPointInside(G4TwoVector a, G4TwoVector b, G4TwoVector c, 
                      G4TwoVector p) const;
      
    G4VFacet* MakeDownFacet(G4int ind1, G4int ind2, G4int ind3) const;      
    G4VFacet* MakeUpFacet(G4int ind1, G4int ind2, G4int ind3) const;      

    G4bool AddGeneralPolygonFacets();
    G4bool MakeFacets();
    G4bool IsConvex() const;

  private:

    G4int       fNv;
    G4double    fHz;
    G4TwoVector fOffset1, fOffset2;
    G4double    fScale1, fScale2;
    std::vector<G4TwoVector> fPolygon;
    std::vector< std::vector<G4int> > fTriangles;
    G4bool          fIsConvex;
    G4GeometryType  fGeometryType;
};    

#include "G4ExtrudedSolid.icc"

#endif
