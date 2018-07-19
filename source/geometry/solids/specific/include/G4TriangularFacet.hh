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
// * subject to DEFCON 705 IPR conditions.                            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4TriangularFacet.hh 95801 2016-02-25 10:59:41Z gcosmo $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class G4TriangularFacet
//
// Class description:
//
//   The G4TriangularFacet class is used for the contruction of
//   G4TessellatedSolid.
//   It is defined by three fVertices, which shall be supplied in anti-clockwise
//   order looking from the outsider of the solid where it belongs.
//   Its constructor:
//   
//      G4TriangularFacet (const G4ThreeVector Pt0, const G4ThreeVector vt1,
//                         const G4ThreeVector vt2, G4FacetVertexType);
//
//   takes 4 parameters to define the three fVertices:
//      1) G4FacetvertexType = "ABSOLUTE": in this case Pt0, vt1 and vt2 are 
//         the 3 fVertices in anti-clockwise order looking from the outsider.
//      2) G4FacetvertexType = "RELATIVE": in this case the first vertex is Pt0,
//         the second vertex is Pt0+vt1 and the third vertex is Pt0+vt2, all  
//         in anti-clockwise order when looking from the outsider.

// CHANGE HISTORY
// --------------
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK - Created.
// 12 October 2012, M Gayer, CERN, - Reviewed optimized implementation.
//
///////////////////////////////////////////////////////////////////////////////
#ifndef G4TriangularFacet_hh
#define G4TriangularFacet_hh 1

#include "G4VFacet.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"

class G4TriangularFacet : public G4VFacet
{
  public:  // with desctiption

    G4TriangularFacet ();
   ~G4TriangularFacet ();

    G4TriangularFacet (const G4ThreeVector &vt0, const G4ThreeVector &vt1,
                       const G4ThreeVector &vt2, G4FacetVertexType);
    G4TriangularFacet (const G4TriangularFacet &right);

    G4TriangularFacet &operator=(const G4TriangularFacet &right);    

    G4VFacet *GetClone ();
    G4TriangularFacet *GetFlippedFacet ();

    G4ThreeVector Distance (const G4ThreeVector &p);
    G4double Distance (const G4ThreeVector &p, G4double minDist);
    G4double Distance (const G4ThreeVector &p, G4double minDist,
                       const G4bool outgoing);
    G4double Extent   (const G4ThreeVector axis);
    G4bool Intersect  (const G4ThreeVector &p, const G4ThreeVector &v,
                       const G4bool outgoing, G4double &distance,
                             G4double &distFromSurface, G4ThreeVector &normal);
    G4double GetArea () const;
    G4ThreeVector GetPointOnFace () const;

    G4ThreeVector GetSurfaceNormal () const;
    void SetSurfaceNormal (G4ThreeVector normal);

    G4GeometryType GetEntityType () const;

    inline G4bool IsDefined () const;
    inline G4int GetNumberOfVertices () const;
    inline G4ThreeVector GetVertex (G4int i) const;
    inline void SetVertex (G4int i, const G4ThreeVector &val);

    inline G4ThreeVector GetCircumcentre () const;
    inline G4double GetRadius () const;

    inline G4int AllocatedMemory();

    inline G4int GetVertexIndex (G4int i) const;
    inline void SetVertexIndex (G4int i, G4int j);
    inline void SetVertices(std::vector<G4ThreeVector> *v);

  private:

    void CopyFrom(const G4TriangularFacet &rhs);

    G4ThreeVector fSurfaceNormal;
    G4double fArea;
    G4ThreeVector fCircumcentre;
    G4double fRadius;
    G4int fIndices[3];

    std::vector<G4ThreeVector> *fVertices;

    G4double fA, fB, fC;
    G4double fDet;
    G4double fSqrDist;
    G4ThreeVector fE1, fE2;
    G4bool fIsDefined;
};

///////////////////////////////////////////////////////////////////////////////
// Inlined Methods
///////////////////////////////////////////////////////////////////////////////

inline G4bool G4TriangularFacet::IsDefined () const
{
  return fIsDefined;
}

inline G4int G4TriangularFacet::GetNumberOfVertices () const
{
  return 3;
}

inline G4ThreeVector G4TriangularFacet::GetVertex (G4int i) const
{      
  G4int indice = fIndices[i];
  return indice < 0 ? (*fVertices)[i] : (*fVertices)[indice];
}

inline void G4TriangularFacet::SetVertex (G4int i, const G4ThreeVector &val)
{
  (*fVertices)[i] = val;
}

inline G4ThreeVector G4TriangularFacet::GetCircumcentre () const
{
  return fCircumcentre;
}

inline G4double G4TriangularFacet::GetRadius () const
{
  return fRadius;
}

inline G4int G4TriangularFacet::AllocatedMemory()
{
  G4int size = sizeof(*this);
  size += GetNumberOfVertices() * sizeof(G4ThreeVector);
  return size;
}

inline G4int G4TriangularFacet::GetVertexIndex (G4int i) const
{
  if (i < 3) return fIndices[i];
  else       return 999999999;
}

inline void G4TriangularFacet::SetVertexIndex (G4int i, G4int j)
{
  fIndices[i] = j;
}

inline void G4TriangularFacet::SetVertices(std::vector<G4ThreeVector> *v)
{
  if (fIndices[0] < 0 && fVertices)
  {
    delete fVertices;
    fVertices = 0;
  }
  fVertices = v;
}

#endif
