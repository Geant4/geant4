// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		G4TriangularFacet.hh
//
// Date:		15/06/2005
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		UK Ministry of Defence : RAO CRP TD Electronic Systems
// Contract:		C/MAT/N03517
//
// This software is the intelectual property of QinetiQ Ltd, subject
// DEFCON 705 IPR conditions.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DISCLAIMER
// ----------
//
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
//    The G4TriangularFacet class is used for the contruction of G4TessellatedSolid.
//    It is defined by three vertices, which shall be supplied 
//    in anti-clockwise order looking from the outsider of the solid where 
//    it belongs.Its constructor
//   
//      G4TriangularFacet (const G4ThreeVector Pt0, const G4ThreeVector vt1,
//          const G4ThreeVector vt2, G4FacetVertexType);
//
//    takes 4 parameters to define the three vertices:
//          1) G4FacetvertexType = "ABSOLUTE": in this case Pt0, vt1 and vt2 are 
//             the three vertices in anti-clockwise order looking from the outsider.
//          2) G4FacetvertexType = "RELATIVE": in this case the first vertex is Pt0,
//             the second vertex is Pt0+vt1 and the third vertex is Pt0+vt2, all  
//             in anti-clockwise order when looking from the outsider.
//
///////////////////////////////////////////////////////////////////////////////
//
//
#ifndef G4TriangularFacet_hh
#define G4TriangularFacet_hh 1

#include "G4VFacet.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
using namespace std;

class G4TriangularFacet : public G4VFacet
{
  public: 
    G4TriangularFacet (const G4ThreeVector Pt0, const G4ThreeVector vt1,
      const G4ThreeVector vt2, G4FacetVertexType);
    ~G4TriangularFacet ();
    
    G4TriangularFacet (const G4TriangularFacet &right);
    const G4TriangularFacet &operator=(G4TriangularFacet &right);    

    G4VFacet *GetClone();
    G4TriangularFacet *GetFlippedFacet ();
    
    G4ThreeVector Distance (const G4ThreeVector &p);
//    G4double Distance (const G4ThreeVector &p, const G4bool outgoing);
    G4double Distance (const G4ThreeVector &p, const G4double minDist);
    G4double Distance (const G4ThreeVector &p, const G4double minDist,
      const G4bool outgoing);
    G4double Extent (const G4ThreeVector axis);
    G4bool Intersect (const G4ThreeVector &p, const G4ThreeVector &v,
      const G4bool outgoing, G4double &distance, G4double &distFromSurface,
      G4ThreeVector &normal);
    G4bool IsInside(const G4ThreeVector &p) const;
      
  private:
    G4double a;
    G4double b;
    G4double c;
    G4double det;
    
    G4double sMin, sMax;
    G4double tMin;
};
#endif
///////////////////////////////////////////////////////////////////////////////
//

