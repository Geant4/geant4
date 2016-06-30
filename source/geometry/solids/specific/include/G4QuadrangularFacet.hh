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
// $Id: G4QuadrangularFacet.hh 95801 2016-02-25 10:59:41Z gcosmo $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class G4QuadrangularFacet
//
// Class description:
//
//   The G4QuadrangularFacet class is used for the contruction of
//   G4TessellatedSolid.
//   It is defined by four fVertices, which shall be in the same plane and be
//   supplied in anti-clockwise order looking from the outsider of the solid
//   where it belongs. Its constructor
//   
//     G4QuadrangularFacet (const G4ThreeVector Pt0, const G4ThreeVector vt1,
//                          const G4ThreeVector vt2, const G4ThreeVector vt3,
//                          G4FacetVertexType);
//
//   takes 5 parameters to define the four fVertices:
//     1) G4FacetvertexType = "ABSOLUTE": in this case Pt0, vt1, vt2 and vt3
//        are the four fVertices required in anti-clockwise order when looking
//        from the outsider.
//     2) G4FacetvertexType = "RELATIVE": in this case the first vertex is Pt0,
//        the second vertex is Pt0+vt, the third vertex is Pt0+vt2 and 
//        the fourth vertex is Pt0+vt3, in anti-clockwise order when looking 
//        from the outsider.

// CHANGE HISTORY
// --------------
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK - Created.
// 12 October 2012, M Gayer, CERN, - Reviewed optimized implementation.
//
///////////////////////////////////////////////////////////////////////////////
#ifndef G4QuadrangularFacet_HH
#define G4QuadrangularFacet_HH 1

#include "G4VFacet.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4TriangularFacet.hh"

class G4QuadrangularFacet : public G4VFacet
{
  public:  // with description

    G4QuadrangularFacet (const G4ThreeVector &Pt0, const G4ThreeVector &vt1,
                         const G4ThreeVector &vt2, const G4ThreeVector &vt3,
                               G4FacetVertexType);
    G4QuadrangularFacet (const G4QuadrangularFacet &right);
   ~G4QuadrangularFacet ();

    G4QuadrangularFacet &operator=(const G4QuadrangularFacet &right);    

    G4VFacet *GetClone ();

    G4ThreeVector Distance (const G4ThreeVector &p);
    G4double Distance (const G4ThreeVector &p, G4double minDist);
    G4double Distance (const G4ThreeVector &p, G4double minDist,
                       const G4bool outgoing);
    G4double Extent   (const G4ThreeVector axis);
    G4bool Intersect  (const G4ThreeVector &p, const G4ThreeVector &v,
                       const G4bool outgoing, G4double &distance,
                             G4double &distFromSurface, G4ThreeVector &normal);
    G4ThreeVector GetSurfaceNormal () const;

    G4double GetArea () const;
    G4ThreeVector GetPointOnFace () const;

    G4GeometryType GetEntityType () const;

    inline G4bool IsDefined () const;
    inline G4int GetNumberOfVertices () const;
    inline G4ThreeVector GetVertex (G4int i) const;
    inline void SetVertex (G4int i, const G4ThreeVector &val);
    inline void SetVertices(std::vector<G4ThreeVector> *v);

    inline G4double GetRadius () const;
    inline G4ThreeVector GetCircumcentre () const;

  private:

    inline G4int GetVertexIndex (G4int i) const;
    inline void SetVertexIndex (G4int i, G4int val);

    inline G4int AllocatedMemory();

  private:

    G4double fRadius;
    G4ThreeVector fCircumcentre;

    G4TriangularFacet fFacet1, fFacet2;
};

///////////////////////////////////////////////////////////////////////////////
// Inlined Methods
///////////////////////////////////////////////////////////////////////////////

inline G4int G4QuadrangularFacet::GetNumberOfVertices () const
{
  return 4;
}

inline G4ThreeVector G4QuadrangularFacet::GetVertex (G4int i) const
{
  return i == 3 ? fFacet2.GetVertex(2) : fFacet1.GetVertex(i);
}


inline G4double G4QuadrangularFacet::GetRadius () const
{
  return fRadius;
}

inline G4ThreeVector G4QuadrangularFacet::GetCircumcentre () const
{
  return fCircumcentre;
}

inline void G4QuadrangularFacet::SetVertex (G4int i, const G4ThreeVector &val)
{
  switch (i)
  {
    case 0:
      fFacet1.SetVertex(0, val);
      fFacet2.SetVertex(0, val);
      break;
    case 1:
      fFacet1.SetVertex(1, val);
      break;
    case 2:
      fFacet1.SetVertex(2, val);
      fFacet2.SetVertex(1, val);
      break;
    case 3:
      fFacet2.SetVertex(2, val);
      break;
    }
}

inline void G4QuadrangularFacet::SetVertices(std::vector<G4ThreeVector> *v)
{
  fFacet1.SetVertices(v);
  fFacet2.SetVertices(v);
}

inline G4bool G4QuadrangularFacet::IsDefined () const
{
  return fFacet1.IsDefined();
}

inline G4int G4QuadrangularFacet::GetVertexIndex (G4int i) const
{
  return i == 3 ? fFacet2.GetVertexIndex(2) : fFacet1.GetVertexIndex(i);
}


inline void G4QuadrangularFacet::SetVertexIndex (G4int i, G4int val)
{
  switch (i)
  {
    case 0:
      fFacet1.SetVertexIndex(0, val);
      fFacet2.SetVertexIndex(0, val);
      break;
    case 1:
      fFacet1.SetVertexIndex(1, val);
      break;
    case 2:
      fFacet1.SetVertexIndex(2, val);
      fFacet2.SetVertexIndex(1, val);
      break;
    case 3:
      fFacet2.SetVertexIndex(2, val);
      break;
  }
}

inline G4int G4QuadrangularFacet::AllocatedMemory()
{
  return sizeof(*this) + fFacet1.AllocatedMemory() + fFacet2.AllocatedMemory();
}

#endif
