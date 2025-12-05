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
// G4VFacet
//
// Class description:
//
// Base class defining the facets which are components of a
// G4TessellatedSolid shape.

// Author: P.R.Truscott (QinetiQ Ltd, UK), 31.10.2004 - Created.
//         M.Gayer (CERN), 12.10.2012 - Reviewed optimised implementation.
// --------------------------------------------------------------------
#ifndef G4VFACET_HH
#define G4VFACET_HH

#include <iostream>
#include <vector>

#include "globals.hh"
#include "windefs.hh"
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"

enum G4FacetVertexType { ABSOLUTE, RELATIVE };

/**
 * @brief G4VFacet is a base class defining the facets which are components
 * of a G4TessellatedSolid shape.
 */

class G4VFacet
{
  public:

    /**
     * Constructor and default Destructor.
     */
    G4VFacet();
    virtual ~G4VFacet() = default;

    /**
     * Equality operator.
     */
    G4bool operator== (const G4VFacet& right) const;

    /**
     * Returns the number of vertices of the facet.
     */
    virtual G4int GetNumberOfVertices () const = 0;

    /**
     * Returns the vertex based on the index 'i'.
     */
    virtual G4ThreeVector GetVertex (G4int i) const = 0;

    /**
     * Methods to set the vertices.
     */
    virtual void SetVertex (G4int i, const G4ThreeVector& val) = 0;
    virtual void SetVertices(std::vector<G4ThreeVector>* vertices) = 0;

    /**
     * Returns the type ID of the facet.
     */
    virtual G4GeometryType GetEntityType () const = 0;

    /**
     * Returns the normal vector to the facet.
     */
    virtual G4ThreeVector GetSurfaceNormal () const = 0;

    /**
     * Returns true if the facet is defined.
     */
    virtual G4bool IsDefined () const = 0;

    /**
     * Returns the circumcentre point of the facet.
     */
    virtual G4ThreeVector GetCircumcentre () const = 0;

    /**
     * Returns the radius to the anchor point and centered on the circumcentre.
     */
    virtual G4double GetRadius () const = 0;

    /**
     * Returns a pointer to a newly allocated duplicate copy of the facet.
     */
    virtual G4VFacet* GetClone () = 0;

    /**
     * Determines the closest distance between point p and the facet.
     */
    virtual G4double Distance (const G4ThreeVector&, G4double minDist) = 0;

    /**
     * Determines the distance to point 'p'. kInfinity is returned if either:
     * (1) outgoing is TRUE and the dot product of the normal vector to the
     * facet and the displacement vector from p to the triangle is negative.
     * (2) outgoing is FALSE and the dot product of the normal vector to the
     * facet and the displacement vector from p to the triangle is positive.
     */
    virtual G4double Distance (const G4ThreeVector&, G4double minDist,
                               const G4bool) = 0;

    /**
     * Calculates the furthest the triangle extends in fA particular
     * direction defined by the vector axis.
     */
    virtual G4double Extent (const G4ThreeVector axis) = 0;

    /**
     * Finds the next intersection when going from 'p' in the direction of 'v'.
     * If 'outgoing' is true, only consider the face if we are going out
     * through the face; otherwise, if false, only consider the face if we are
     * going in through the face.
     *  @returns true if there is an intersection, false otherwise.
     */
    virtual G4bool Intersect (const G4ThreeVector& p, const G4ThreeVector& v,
                              const G4bool outgoing, G4double& distance,
                                    G4double& distFromSurface,
                                    G4ThreeVector& normal) = 0;

    /**
     * Auxiliary method for returning the surface area.
     */
    virtual G4double GetArea() const = 0;

    /**
     * Auxiliary method to get a uniform random point on the facet.
     */
    virtual G4ThreeVector GetPointOnFace() const = 0;

    /**
     * Adds a translation 'v' to the vertices of the facet.
     */
    void ApplyTranslation (const G4ThreeVector& v);

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo(std::ostream& os) const;

    /**
     * Returns true if point 'p' is inside the facet.
     */
    G4bool IsInside(const G4ThreeVector& p) const;

    /**
     * Logger methods for allocated memory of facets.
     */
    virtual G4int AllocatedMemory() = 0;
    virtual void SetVertexIndex (G4int i, G4int j) = 0;
    virtual G4int GetVertexIndex (G4int i) const = 0;


  protected:

    static const G4double dirTolerance;
    G4double kCarTolerance;
};

#endif
