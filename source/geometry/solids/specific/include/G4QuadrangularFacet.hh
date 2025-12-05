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
// G4QuadrangularFacet
//
// Class description:
//
// The G4QuadrangularFacet class is used for the contruction of
// G4TessellatedSolid.
// It is defined by four fVertices, which shall be in the same plane and be
// supplied in anti-clockwise order looking from the outsider of the solid
// where it belongs. Its constructor:
//   
//   G4QuadrangularFacet (const G4ThreeVector& Pt0, const G4ThreeVector& vt1,
//                        const G4ThreeVector& vt2, const G4ThreeVector& vt3,
//                        G4FacetVertexType);
//
// takes 5 parameters to define the four fVertices:
//   1) G4FacetvertexType = "ABSOLUTE": in this case Pt0, vt1, vt2 and vt3
//      are the four fVertices required in anti-clockwise order when looking
//      from the outsider.
//   2) G4FacetvertexType = "RELATIVE": in this case the first vertex is Pt0,
//      the second vertex is Pt0+vt1, the third vertex is Pt0+vt2 and 
//      the fourth vertex is Pt0+vt3, in anti-clockwise order when looking 
//      from the outsider.

// Author: P.R.Truscott (QinetiQ Ltd, UK), 31.10.2004 - Created
//         M.Gayer (CERN), 12.10.2012 - Reviewed optimised implementation
// --------------------------------------------------------------------
#ifndef G4QUADRANGULARFACET_HH
#define G4QUADRANGULARFACET_HH

#include "G4VFacet.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4TriangularFacet.hh"

/**
 * @brief G4QuadrangularFacet defines a facet with 4 vertices, used for the
 * contruction of G4TessellatedSolid. Vertices shall be in the same plane and
 * be supplied in anti-clockwise order looking from the outsider of the solid
 * where it belongs.
 */

class G4QuadrangularFacet : public G4VFacet
{
  public:

    /**
     * Constructs a facet with 4 vertices, given its parameters.
     *  @param[in] Pt0 The anchor point, first vertex.
     *  @param[in] vt1 Second vertex.
     *  @param[in] vt2 Third vertex.
     *  @param[in] vt3 Fourth vertex.
     *  @param[in] vType The positioning type for the vertices, either:
     *             "ABSOLUTE" - vertices set in anti-clockwise order
     *                          when looking from the outsider.
     *             "RELATIVE" - first vertex is Pt0, second is Pt0+vt1,
     *                          third vertex is Pt0+vt2 and fourth is Pt0+vt3,
     *                          still in anti-clockwise order.
     */
    G4QuadrangularFacet (const G4ThreeVector& Pt0, const G4ThreeVector& vt1,
                         const G4ThreeVector& vt2, const G4ThreeVector& vt3,
                               G4FacetVertexType vType);

    /**
     * Default Destructor.
     */
    ~G4QuadrangularFacet () override = default;

    /**
     * Copy constructor and assignment operator.
     */
    G4QuadrangularFacet (const G4QuadrangularFacet& right);
    G4QuadrangularFacet& operator=(const G4QuadrangularFacet& right);    

    /**
     * Returns a pointer to a newly allocated duplicate copy of the facet.
     */
    G4VFacet* GetClone () override;

    /**
     * Determines the vector between p and the closest point on the facet to p.
     */
    G4ThreeVector Distance (const G4ThreeVector& p);

    /**
     * Determines the closest distance between point p and the facet.
     */
    G4double Distance (const G4ThreeVector& p, G4double minDist) override;

    /**
     * Determines the distance to point 'p'. kInfinity is returned if either:
     * (1) outgoing is TRUE and the dot product of the normal vector to the
     * facet and the displacement vector from p to the triangle is negative.
     * (2) outgoing is FALSE and the dot product of the normal vector to the
     * facet and the displacement vector from p to the triangle is positive.
     */
    G4double Distance (const G4ThreeVector& p, G4double minDist,
                       const G4bool outgoing) override;

    /**
     * Calculates the furthest the quadrangle extends in fA particular
     * direction defined by the vector axis.
     */
    G4double Extent (const G4ThreeVector axis) override;

    /**
     * Finds the next intersection when going from 'p' in the direction of 'v'.
     * If 'outgoing' is true, only consider the face if we are going out
     * through the face; otherwise, if false, only consider the face if we are
     * going in through the face.
     *  @returns true if there is an intersection, false otherwise.
     */
    G4bool Intersect  (const G4ThreeVector& p, const G4ThreeVector& v,
                       const G4bool outgoing, G4double& distance,
                             G4double& distFromSurface,
                             G4ThreeVector& normal) override;
    /**
     * Returns the normal vector to the face.
     */
    G4ThreeVector GetSurfaceNormal () const override;

    /**
     * Auxiliary method for returning the surface area.
     */
    G4double GetArea () const override;

    /**
     * Auxiliary method to get a uniform random point on the facet.
     */
    G4ThreeVector GetPointOnFace () const override;

    /**
     * Returns the type ID, "G4QuadrangularFacet" of the facet.
     */
    G4GeometryType GetEntityType () const override;

    /**
     * Returns true if the facet is defined.
     */
    inline G4bool IsDefined () const override;

    /**
     * Returns the number of vertices, i.e. 4.
     */
    inline G4int GetNumberOfVertices () const override;

    /**
     * Returns the vertex based on the index 'i'.
     */
    inline G4ThreeVector GetVertex (G4int i) const override;

    /**
     * Methods to set the vertices.
     */
    inline void SetVertex (G4int i, const G4ThreeVector& val) override;
    inline void SetVertices(std::vector<G4ThreeVector>* v) override;

    /**
     * Returns the radius to the anchor point and centered on the circumcentre.
     */
    inline G4double GetRadius () const override;

    /**
     * Returns the circumcentre point of the facet.
     */
    inline G4ThreeVector GetCircumcentre () const override;

  private:

    /**
     * Private accessor/setter for the vertex index.
     */
    inline G4int GetVertexIndex (G4int i) const override;
    inline void SetVertexIndex (G4int i, G4int val) override;

    /**
     * Returns the allocated memory (sizeof) by the facet.
     */
    inline G4int AllocatedMemory() override;

  private:

    G4double fRadius = 0.0;

    G4ThreeVector fCircumcentre;

    G4TriangularFacet fFacet1, fFacet2;
};

// --------------------------------------------------------------------
// Inline Methods
// --------------------------------------------------------------------

#include "G4QuadrangularFacet.icc"

#endif
