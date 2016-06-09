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
// $Id: G4VFacet.hh,v 1.8 2010-09-23 10:27:25 gcosmo Exp $
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

    G4bool operator== (const G4VFacet &right) const;

    inline size_t GetNumberOfVertices () const;
    inline G4ThreeVector GetVertex (size_t i) const;
    
    inline G4GeometryType GetEntityType () const;
    inline G4ThreeVector GetSurfaceNormal () const;
    inline G4bool IsInside(const G4ThreeVector &p) const;
    inline G4bool IsDefined () const;
    inline void SetVertexIndex (const size_t i, const size_t j);
    inline size_t GetVertexIndex (const size_t i) const;
    inline G4ThreeVector GetCircumcentre () const;
    inline G4double GetRadius () const;
    inline G4double GetRadiusSquared() const;
    
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
    virtual G4double GetArea() = 0;
    virtual G4ThreeVector GetPointOnFace() const = 0;

  public:  // without description

    G4VFacet (const G4VFacet &right);
    const G4VFacet &operator=(G4VFacet &right);

  protected:

    G4GeometryType       geometryType;
    G4bool               isDefined;
    size_t               nVertices;
    G4ThreeVector        P0;
    G4ThreeVectorList    P;
    G4ThreeVectorList    E;
    std::vector<size_t>  I;
    G4ThreeVector        surfaceNormal;
    G4ThreeVector        circumcentre;
    G4double             radius;
    G4double             radiusSqr;

    G4double             dirTolerance;
    G4double             kCarTolerance;
    G4double             area;
};

typedef std::vector<G4VFacet*>::iterator       FacetI;
typedef std::vector<G4VFacet*>::const_iterator FacetCI;

#include "G4VFacet.icc"

#endif
