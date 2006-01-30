//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration and of QinetiQ Ltd,  subject DEFCON 705 IPR *
// * conditions.                                                      *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4VFacet.hh,v 1.2 2006-01-30 14:39:53 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4VFacet.hh
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
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK - Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class description:
//
//   Base class defining the facets which are components of a
//   G4TessellatedSolid shape.
//

///////////////////////////////////////////////////////////////////////////////
#ifndef G4VFacet_hh
#define G4VFacet_hh 1

#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <iostream>
#include <vector>

enum G4FacetVertexType {ABSOLUTE, RELATIVE};

class G4VFacet
{
  public:  // with description

    G4VFacet ();
    virtual ~G4VFacet ();

    G4VFacet (const G4VFacet &right);
    const G4VFacet &operator=(G4VFacet &right);
    
    G4bool operator== (const G4VFacet &right) const;

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

    virtual G4VFacet *GetClone ();
    virtual G4double Distance (const G4ThreeVector&, const G4double);
    virtual G4double Distance (const G4ThreeVector&, const G4double,
                               const G4bool);
    virtual G4double Extent   (const G4ThreeVector);
    virtual G4bool Intersect  (const G4ThreeVector&, const G4ThreeVector &,
                               const G4bool , G4double &, G4double &,
                                     G4ThreeVector &);
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
