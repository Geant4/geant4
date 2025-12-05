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
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4BooleanSolid
//
// Class description:
//
// Base class for solids created by Boolean operations between other solids.

// Author: Vladimir Grichine (CERN), 10.09.1998 - Created.
// --------------------------------------------------------------------
#ifndef G4BOOLEANSOLID_HH
#define G4BOOLEANSOLID_HH

#include "G4DisplacedSolid.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

#include "G4VBooleanProcessor.hh"

class HepPolyhedronProcessor;

/**
 * @brief G4BooleanSolid is the base class for solids created by Boolean
 * operations between other solids.
 */

class G4BooleanSolid : public G4VSolid
{
  public:
 
    /**
     * Constructor of a Boolean composition between two solids with no
     * displacement.
     *  @param[in] pName The name of the Boolean composition.
     *  @param[in] pSolidA Pointer to the first reference solid.
     *  @param[in] pSolidB Pointer to the second solid to form the composition.
     */
    G4BooleanSolid( const G4String& pName,
                          G4VSolid* pSolidA ,
                          G4VSolid* pSolidB   );

    /**
     * Constructor of a Boolean composition between two solids with rotation
     * and translation, used to transform the coordinate system of the second
     * solid to the coordinate system of the first solid.
     *  @param[in] pName The name of the Boolean composition.
     *  @param[in] pSolidA Pointer to the first reference solid.
     *  @param[in] pSolidB Pointer to the second solid to form the composition.
     *  @param[in] rotMatrix Pointer to the rotation vector.
     *  @param[in] transVector The translation vector.
     */
    G4BooleanSolid( const G4String& pName,
                          G4VSolid* pSolidA ,
                          G4VSolid* pSolidB,
                          G4RotationMatrix* rotMatrix,
                    const G4ThreeVector& transVector    );

    /**
     * Constructor of a Boolean composition between two solids with a
     * transformation that moves the second solid from its desired position
     * to its standard position.
     *  @param[in] pName The name of the Boolean composition.
     *  @param[in] pSolidA Pointer to the first reference solid.
     *  @param[in] pSolidB Pointer to the second solid to form the composition.
     *  @param[in] transform The composed 3D transformation.
     */
    G4BooleanSolid( const G4String& pName,
                          G4VSolid* pSolidA ,
                          G4VSolid* pSolidB , 
                    const G4Transform3D& transform   );

    /**
     * Destructor. If using a displaced solid, deletes all cached
     * transformations.
     */
    ~G4BooleanSolid() override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4BooleanSolid(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4BooleanSolid(const G4BooleanSolid& rhs);
    G4BooleanSolid& operator=(const G4BooleanSolid& rhs);

    /**
     * Methods returning the component solids of the Boolean composition.
     * If the solid is made up from a Boolean operation of two solids,
     * return the corresponding solid (for no=0 and 1). A fatal exception
     * is thrown if the index provided is different from 0 or 1.
     *  @param[in] no Index 0/1 of the components.
     *  @returns The pointer to (const or not const) of the component solid.
     *           If the solid is not a "Boolean", returns nullptr.
     */
    const G4VSolid* GetConstituentSolid(G4int no) const override;
          G4VSolid* GetConstituentSolid(G4int no) override;

    /**
     * Methods returning the computed capacity and surface area of the
     * composition. The quantities returned are an estimate obtained by
     * randomly sampling the Boolean composition and caching them for reuse.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    /**
     * Returns the type ID, "G4BooleanSolid" of the solid.
     */
    G4GeometryType GetEntityType() const override;

    /**
     * Returns a pointer to the generated polyhedron representation of the
     * composition, for use in visualisation.
     */
    G4Polyhedron* GetPolyhedron() const override;

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo(std::ostream& os) const override;

    /**
     * Accessor and setter for controlling/tuning the number of random points
     * to be used for computing the cubic volume.
     */
    inline G4int GetCubVolStatistics() const;
    void SetCubVolStatistics(G4int st);

    /**
     * Accessor and setter for controlling/tuning the epsilon used for
     * computing the cubic volume.
     */
    inline G4double GetCubVolEpsilon() const;
    void SetCubVolEpsilon(G4double ep);
   
    /**
     * Accessor and setter for controlling/tuning the number of random points
     * to be used for computing the surface area.
     */
    inline G4int GetAreaStatistics() const;
    inline void SetAreaStatistics(G4int st);

    /**
     * Accessor and setter for controlling/tuning the level of accuracy used
     * for computing the surface area.
     */
    inline G4double GetAreaAccuracy() const;
    inline void SetAreaAccuracy(G4double ep);

    /**
     * Returns a point (G4ThreeVector) randomly and uniformly generated
     * on the surface of the solid.
     */
    G4ThreeVector GetPointOnSurface() const override;

    /**
     * Returns the total number of constituent solids forming the Boolean
     * composition.
     */
    G4int GetNumOfConstituents() const override;

    /**
     * Return true if the resulting solid has only planar faces.
     */
    G4bool IsFaceted() const override;

    /**
     * Gets/sets the Boolean processor for polyhedron to replace the default
     * processor.
     */
    static G4VBooleanProcessor* GetExternalBooleanProcessor();
    static void SetExternalBooleanProcessor(G4VBooleanProcessor* extProcessor);

  protected:
  
    /**
     * Gets the list of constituent primitives of the solid and their placements.
     */
    void GetListOfPrimitives(std::vector<std::pair<G4VSolid *,G4Transform3D>>&,
                             const G4Transform3D&) const;
      // 

    /**
     * Stacks the polyhedra for processing.
     *  @returns A pointer to the top polyhedron.
     */
    G4Polyhedron* StackPolyhedron(HepPolyhedronProcessor&,
                                  const G4VSolid*) const;

  protected:
  
    /** Pointers to the costituent solids. */
    G4VSolid* fPtrSolidA = nullptr;
    G4VSolid* fPtrSolidB = nullptr;

    /** Cached value of the capacity. */
    G4double fCubicVolume = -1.0;

    /** Cached value of the surface area. */
    G4double fSurfaceArea = -1.0;

    /** Static pointer to the external Boolean processor. */
    static G4VBooleanProcessor* fExternalBoolProcessor;

  private:

    G4int fCubVolStatistics = 1000000;
    G4int fAreaStatistics = 1000000;
    G4double fCubVolEpsilon = 0.001;
    G4double fAreaAccuracy = -1;

    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;

    mutable std::vector<std::pair<G4VSolid *,G4Transform3D>> fPrimitives;
    mutable G4double fPrimitivesSurfaceArea = 0.0;

    /** If & only if this object created it, it must delete it. */
    G4bool createdDisplacedSolid = false;
};

#include "G4BooleanSolid.icc"

#endif
