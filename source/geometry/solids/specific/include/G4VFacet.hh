// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		G4VFacet.hh
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
//
//
///////////////////////////////////////////////////////////////////////////////
//
//
#ifndef G4VFacet_hh
#define G4VFacet_hh 1

#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <iostream>
#include <vector>
using namespace std;
enum G4FacetVertexType {ABSOLUTE, RELATIVE};

class G4VFacet
{
  public:
    G4VFacet ();
    virtual ~G4VFacet ();

    G4VFacet (const G4VFacet &right);
    const G4VFacet &operator=(G4VFacet &right);
    
    G4bool operator== (const G4VFacet &right) const;
                                                                                    
  public:
    size_t GetNumberOfVertices () const;
    G4ThreeVector GetVertex (size_t i) const;
    
    G4GeometryType GetEntityType () const;
    G4ThreeVector GetSurfaceNormal () const;
    G4bool IsInside(const G4ThreeVector &p) const;
    G4bool IsDefined () const;
    void SetVertexIndex (const size_t i, const size_t j);
    size_t GetVertexIndex (const size_t i) const;
    G4ThreeVector GetCentroid () const;
    G4double GetRadius () const;
    G4double GetRadiusSquared() const;
    
    void ApplyTranslation (const G4ThreeVector v);
    
    std::ostream &StreamInfo(std::ostream &os) const;

    virtual G4VFacet *GetClone () {return NULL;};
//    virtual G4double Distance (const G4ThreeVector &p, const G4double minDist) {return kInfinity;};
    virtual G4double Distance (const G4ThreeVector&, const G4double) {return kInfinity;};
//    virtual G4double Distance (const G4ThreeVector &p, const G4bool outgoing) {return kInfinity;};
//    virtual G4double Distance (const G4ThreeVector &p, const G4double minDist,
//      const G4bool outgoing) {return kInfinity;};
    virtual G4double Distance (const G4ThreeVector&, const G4double,
      const G4bool) {return kInfinity;};
//    virtual G4double Extent (const G4ThreeVector axis) {return 0.0;};
    virtual G4double Extent (const G4ThreeVector) {return 0.0;};
//    virtual G4bool Intersect (const G4ThreeVector &p, const G4ThreeVector &v,
//      const G4bool outgoing, G4double &distance, G4double &distFromSurface,
//      G4ThreeVector &normal) {return false;};
    virtual G4bool Intersect (const G4ThreeVector&, const G4ThreeVector &,
      const G4bool , G4double &, G4double &,
      G4ThreeVector &) {return false;};
  
//  public:
  protected:
    G4GeometryType       geometryType;
    G4bool               isDefined;
    size_t               nVertices;
    G4ThreeVector        P0;
    G4ThreeVectorList    P;
    G4ThreeVectorList    E;
    std::vector<size_t>  I;
    G4ThreeVector        surfaceNormal;
    G4ThreeVector        centroid;
    G4double             radius;
    G4double             radiusSqr;
    
    G4double             dirTolerance;
};

typedef std::vector<G4VFacet*>::iterator       FacetI;
typedef std::vector<G4VFacet*>::const_iterator FacetCI;

#include "G4VFacet.icc"

#endif
///////////////////////////////////////////////////////////////////////////////
//

