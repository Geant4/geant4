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
// $Id: G4VFacet.hh 95801 2016-02-25 10:59:41Z gcosmo $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class G4VFacet
//
// Class description:
//
//   Base class defining the facets which are components of a
//   G4TessellatedSolid shape.

// CHANGE HISTORY
// --------------
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK - Created.
// 12 October 2012, M Gayer, CERN, - Reviewed optimized implementation.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef G4VFacet_hh
#define G4VFacet_hh

#include <iostream>
#include <vector>

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"

enum G4FacetVertexType {ABSOLUTE, RELATIVE};

class G4TessellatedSolid;

class G4VFacet
{
  public:

    virtual ~G4VFacet ();

    G4bool operator== (const G4VFacet &right) const;

    virtual G4int GetNumberOfVertices () const = 0;
    virtual G4ThreeVector GetVertex (G4int i) const = 0;
    virtual void SetVertex (G4int i, const G4ThreeVector &val) = 0;
    virtual G4GeometryType GetEntityType () const = 0;
    virtual G4ThreeVector GetSurfaceNormal () const = 0;
    virtual G4bool IsDefined () const = 0;
    virtual G4ThreeVector GetCircumcentre () const = 0;
    virtual G4double GetRadius () const = 0;
    virtual G4VFacet *GetClone () = 0;
    virtual G4double Distance (const G4ThreeVector&, G4double) = 0;
    virtual G4double Distance (const G4ThreeVector&, G4double,
                               const G4bool) = 0;
    virtual G4double Extent (const G4ThreeVector) = 0;
    virtual G4bool Intersect (const G4ThreeVector&, const G4ThreeVector &,
                              const G4bool , G4double &, G4double &,
                                    G4ThreeVector &) = 0;
    virtual G4double GetArea() const = 0;
    virtual G4ThreeVector GetPointOnFace() const = 0;

    void ApplyTranslation (const G4ThreeVector v);

    std::ostream &StreamInfo(std::ostream &os) const;

    G4bool IsInside(const G4ThreeVector &p) const;

    virtual G4int AllocatedMemory() = 0;
    virtual void SetVertexIndex (G4int i, G4int j) = 0;
    virtual G4int GetVertexIndex (G4int i) const = 0;

    virtual void SetVertices(std::vector<G4ThreeVector> *vertices) = 0;

  protected:

    static const G4double dirTolerance;
    static const G4double kCarTolerance;
};

#endif
