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
// G4TriangularFacet
//
// Class description:
//
// The G4TriangularFacet class is used for the contruction of G4TessellatedSolid.
// It is defined by three fVertices, which shall be supplied in anti-clockwise
// order looking from the outsider of the solid where it belongs.
// Its constructor:
//   
//    G4TriangularFacet (const G4ThreeVector Pt0, const G4ThreeVector vt1,
//                       const G4ThreeVector vt2, G4FacetVertexType);
//
// takes 4 parameters to define the three fVertices:
//    1) G4FacetvertexType = "ABSOLUTE": in this case Pt0, vt1 and vt2 are 
//       the 3 fVertices in anti-clockwise order looking from the outsider.
//    2) G4FacetvertexType = "RELATIVE": in this case the first vertex is Pt0,
//       the second vertex is Pt0+vt1 and the third vertex is Pt0+vt2, all  
//       in anti-clockwise order when looking from the outsider.

// Author: P.R.Truscott (QinetiQ Ltd, UK), 31.10.2004 - Created
//         M.Gayer (CERN), 12.10.2012 - Reviewed optimised implementation
// --------------------------------------------------------------------
#ifndef G4TRIANGULARFACET_HH
#define G4TRIANGULARFACET_HH

#include "G4VFacet.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"

#include <vector>
#include <array>

/**
 * @brief G4TriangularFacet defines a facet with 3 vertices, used for the
 * contruction of G4TessellatedSolid. Vertices shall be supplied in
 * anti-clockwise order looking from the outsider of the solid where it belongs.
 */

class G4TriangularFacet : public G4VFacet
{
  public:

    /**
     * Default Constructor.
     */
    G4TriangularFacet ();

    /**
     * Constructs a facet with 3 vertices, given its parameters.
     *  @param[in] Pt0 The anchor point, first vertex.
     *  @param[in] vt1 Second vertex.
     *  @param[in] vt2 Third vertex.
     *  @param[in] vType The positioning type for the vertices, either:
     *             "ABSOLUTE" - vertices set in anti-clockwise order
     *                          when looking from the outsider.
     *             "RELATIVE" - first vertex is Pt0, second is Pt0+vt1,
     *                          and the third vertex is Pt0+vt2,
     *                          still in anti-clockwise order.
     */
    G4TriangularFacet (const G4ThreeVector& Pt0, const G4ThreeVector& vt1,
                       const G4ThreeVector& vt2, G4FacetVertexType vType);

    /**
     * Destructor.
     */
    ~G4TriangularFacet () override;

    /**
     * Copy and move constructors.
     */
    G4TriangularFacet (const G4TriangularFacet& right);
    G4TriangularFacet (      G4TriangularFacet&& right) noexcept ;

    /**
     * Assignment and move assignment operators.
     */
    G4TriangularFacet& operator=(const G4TriangularFacet& right);    
    G4TriangularFacet& operator=(      G4TriangularFacet&& right) noexcept ;    

    /**
     * Returns a pointer to a newly allocated duplicate copy of the facet.
     */
    G4VFacet* GetClone () override;

    /**
     * Generates and returns an identical facet, but with the normal vector
     * pointing at 180 degrees.
     */
    G4TriangularFacet* GetFlippedFacet ();

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
     * Calculates the furthest the triangle extends in fA particular
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
     * Auxiliary method for returning the surface area.
     */
    G4double GetArea () const override;

    /**
     * Auxiliary method to get a uniform random point on the facet.
     */
    G4ThreeVector GetPointOnFace () const override;

    /**
     * Returns/sets the normal vector to the facet.
     */
    G4ThreeVector GetSurfaceNormal () const override;
    void SetSurfaceNormal (const G4ThreeVector& normal);

    /**
     * Returns the type ID, "G4TriangularFacet" of the facet.
     */
    G4GeometryType GetEntityType () const override;

    /**
     * Returns true if the facet is defined.
     */
    inline G4bool IsDefined () const override;

    /**
     * Returns the number of vertices, i.e. 3.
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

    /**
     * Returns the allocated memory (sizeof) by the facet.
     */
    inline G4int AllocatedMemory() override;

    /**
     * Accessor/setter for the vertex index.
     */
    inline G4int GetVertexIndex (G4int i) const override;
    inline void SetVertexIndex (G4int i, G4int j) override;

  private:

    /**
     * Utilities for copying/moving data.
     */
    void CopyFrom(const G4TriangularFacet& rhs);
    void MoveFrom(G4TriangularFacet& rhs);

    G4ThreeVector fSurfaceNormal;
    G4double fArea = 0.0;
    G4ThreeVector fCircumcentre;
    G4double fRadius = 0.0;
    std::array<G4int, 3> fIndices;
    std::vector<G4ThreeVector>* fVertices = nullptr;

    G4double fA, fB, fC;
    G4double fDet;
    G4double fSqrDist = 0.0;
    G4ThreeVector fE1, fE2;
    G4bool fIsDefined = false;
};

// --------------------------------------------------------------------
// Inline Methods
// --------------------------------------------------------------------

#include "G4TriangularFacet.icc"

#endif
