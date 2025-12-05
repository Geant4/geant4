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
// G4VTwistSurface
//
// Class description:
//
// Abstract base class for boundary surface of G4VSolid.

// Author: Kotoyo Hoshina (Chiba University), 01.08.2002 - Created.
//         Oliver Link (CERN), 13.11.2003 - Integration in Geant4
//               from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------
#ifndef G4VTWISTSURFACE_HH
#define G4VTWISTSURFACE_HH

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VSolid.hh"
#include "geomdefs.hh"

#include "G4RotationMatrix.hh"

#define G4VSURFACENXX 10

/**
 * @brief G4VTwistSurface is a base class for boundary surface of a G4VSolid.
 */

class G4VTwistSurface
{
  public:

    enum EValidate { kDontValidate = 0, kValidateWithTol = 1, 
                     kValidateWithoutTol = 2, kUninitialized = 3 };

    /**
     * Base constructor, given a name.
     *  @param[in] name The surface name.
     */
    G4VTwistSurface (const G4String& name);

    /**
     * Constructs a surface, given its parameters.
     *  @param[in] name The surface name.
     *  @param[in] rot Rotation: 0.5*(phi-width segment).
     *  @param[in] tlate Translation.
     *  @param[in] handedness Orientation: R-hand = 1, L-hand = -1.
     *  @param[in] axis0 First axis.
     *  @param[in] axis1 Second axis.
     *  @param[in] axis0min Minimum in first axis.
     *  @param[in] axis1min Minimum in second axis.
     *  @param[in] axis0max Maximum in first axis.
     *  @param[in] axis1max Maximum in second axis.
     */
    G4VTwistSurface (const G4String& name,
                     const G4RotationMatrix& rot,
                     const G4ThreeVector&    tlate,
                           G4int             handedness,
                     const EAxis             axis0,
                     const EAxis             axis1,
                           G4double          axis0min = -kInfinity,
                           G4double          axis1min = -kInfinity,
                           G4double          axis0max = kInfinity,
                           G4double          axis1max = kInfinity);
 
    /**
     * Default destructor.
     */
    virtual ~G4VTwistSurface() = default;

    /**
     * Returns the phi-location of point 'me' (phi relation between 'me' and
     * 'vec' projected on z=0 plane). If 'me' is on -ve-phi-side of 'vec', it
     * returns 1; on the other hand, if 'me' is on +ve-phi-side of 'vec',
     * it returns -1. If 'me' is on boundary of 'vec', it returns 0.
     *  @returns The z-coordinate of normal vector of me.cross(vec).
     */
    virtual G4int AmIOnLeftSide(const G4ThreeVector& me, 
                                const G4ThreeVector& vec, 
                                      G4bool withTol = true);

    /**
     * Returns the distance to the nearest boundary from an arbitrary point 'p' 
     * in local coodinate.
     *  @param[in] areacode Point location code according to axes.
     *  @param[out] xx The returned reference point on boundary.
     *  @param[in] p The arbitrary point in input.
     */
    virtual G4double DistanceToBoundary(G4int areacode,
                                        G4ThreeVector& xx,
                                        const G4ThreeVector& p) ;

    /**
     * Returns the distance from outside to a boundary from an arbitrary
     * point 'gp' and direction 'gv'.
     *  @param[in] gp The arbitrary point.
     *  @param[in] gv The direction vector.
     *  @param[out] gxxbest The returned reference point on boundary.
     */
    virtual G4double DistanceToIn(const G4ThreeVector& gp,
                                  const G4ThreeVector& gv,
                                        G4ThreeVector& gxxbest);

    /**
     * Returns the distance from inside to a boundary from an arbitrary
     * point 'gp' and direction 'gv'.
     *  @param[in] gp The arbitrary point.
     *  @param[in] gv The direction vector.
     *  @param[out] gxxbest The returned reference point on boundary.
     */
    virtual G4double  DistanceToOut(const G4ThreeVector& gp,
                                    const G4ThreeVector& gv,
                                          G4ThreeVector& gxxbest);

    /**
     * Returns the safety distance to a boundary from an arbitrary point 'gp'.
     *  @param[in] gp The arbitrary point.
     *  @param[out] gxxbest The returned reference point on boundary.
     */
    virtual G4double DistanceTo(const G4ThreeVector& gp,
                                      G4ThreeVector& gxx);
      
    /**
     * Returns the distance to surface, given point 'gp' and direction 'gv'.
     *  @param[in] gp The point from where computing the distance.
     *  @param[in] gv The direction along which computing the distance.
     *  @param[out] gxx Vector of global points based on number of solutions.
     *  @param[out] distance The distance vector based on number of solutions.
     *  @param[out] areacode The location vector based on number of solutions.
     *  @param[out] isvalid Validity vector based on number of solutions.
     *  @param[in] validate Adopted validation criteria.
     *  @returns The number of solutions.
     */
    virtual G4int DistanceToSurface(const G4ThreeVector& gp,
                                    const G4ThreeVector& gv,
                                          G4ThreeVector  gxx[],
                                          G4double       distance[],
                                          G4int          areacode[],
                                          G4bool         isvalid[],
                                    EValidate validate = kValidateWithTol) = 0;

    /**
     * Returns the safety distance to surface, given point 'gp'.
     *  @param[in] gp The point from where computing the safety distance.
     *  @param[out] gxx Vector of global points based on number of solutions.
     *  @param[out] distance The distance vector based on number of solutions.
     *  @param[out] areacode The location vector based on number of solutions.
     *  @returns The number of solutions.
     */
    virtual G4int DistanceToSurface(const G4ThreeVector& gp,
                                          G4ThreeVector  gxx[],
                                          G4double       distance[],
                                          G4int          areacode[]) = 0;
                                             
    /**
     * Returns a normal vector at a surface (or very close to the surface)
     * point at 'p'.
     *  @param[in] p The point where computing the normal.
     *  @param[in] isGlobal If true, it returns the normal in global coordinates.
     *  @returns The normal vector.
     */
    virtual G4ThreeVector GetNormal(const G4ThreeVector& p,
                                          G4bool isGlobal) = 0;
   
    /**
     * Returns parameters of boundaries.
     *  @param[in] areacode The location code.
     *  @param[out] d Boundary direction.
     *  @param[out] x0 The point on boundary.
     *  @param[out] boundarytype The boundary code.
     */
    virtual void GetBoundaryParameters(const G4int& areacode,
                                             G4ThreeVector& d,
                                             G4ThreeVector& x0,
                                             G4int& boundarytype) const;

    /**
     * Returns Z projection of point 'p' on boundary.
     *  @param[in] areacode The location code.
     *  @param[out] p The arbitrary point.
     *  @returns The Z projection on boundary.
     */
    virtual G4ThreeVector GetBoundaryAtPZ(G4int areacode,
                                          const G4ThreeVector& p) const;

    /**
     * Utility methods to compute the distance from a plane/line.
     */
    inline G4double DistanceToPlaneWithV(const G4ThreeVector& p,
                                         const G4ThreeVector& v,
                                         const G4ThreeVector& x0,
                                         const G4ThreeVector& n0,
                                               G4ThreeVector& xx);
    inline G4double DistanceToPlane(const G4ThreeVector& p,
                                    const G4ThreeVector& x0,
                                    const G4ThreeVector& n0,
                                          G4ThreeVector& xx);
    inline G4double DistanceToPlane(const G4ThreeVector& p,
                                    const G4ThreeVector& x0,
                                    const G4ThreeVector& t1,
                                    const G4ThreeVector& t2,
                                          G4ThreeVector& xx,
                                          G4ThreeVector& n);
    inline G4double DistanceToLine (const G4ThreeVector& p,
                                    const G4ThreeVector& x0,
                                    const G4ThreeVector& d,
                                          G4ThreeVector& xx);

    /**
     * Test functions/accessors, based on 'areacode' location.
     */
    inline G4bool IsAxis0 (G4int areacode) const;
    inline G4bool IsAxis1 (G4int areacode) const;
    inline G4bool IsOutside (G4int areacode) const;
    inline G4bool IsInside (G4int areacode, G4bool testbitmode = false) const;
    inline G4bool IsBoundary (G4int areacode, G4bool testbitmode = false) const;
    inline G4bool IsCorner (G4int areacode, G4bool testbitmode = false) const;
    inline G4bool IsValidNorm() const { return fIsValidNorm; }
    G4bool IsSameBoundary (G4VTwistSurface* surface1, G4int areacode1,
                           G4VTwistSurface* surface2, G4int areacode2 ) const;
    inline G4int GetAxisType(G4int areacode, G4int whichaxis) const;

    /**
     * Utility functions computing global/local point/direction.
     */
    inline G4ThreeVector ComputeGlobalPoint (const G4ThreeVector& lp) const;
    inline G4ThreeVector ComputeLocalPoint (const G4ThreeVector& gp) const;
    inline G4ThreeVector ComputeGlobalDirection (const G4ThreeVector& lp) const;
    inline G4ThreeVector ComputeLocalDirection (const G4ThreeVector& gp) const;
  
    /**
     * Modifiers for axes.
     */
    inline void SetAxis(G4int i, const EAxis axis) { fAxis[i] = axis; }
    inline void SetNeighbours(G4VTwistSurface* ax0min, G4VTwistSurface* ax1min, 
                              G4VTwistSurface* ax0max, G4VTwistSurface* ax1max);

    /**
     * Returns point on surface given the axes.
     */
    virtual G4ThreeVector SurfacePoint(G4double, G4double,
                                       G4bool isGlobal = false ) = 0 ;

    /**
     * Internal accessors to be implemented in concrete derived classes.
     */
    virtual G4double GetBoundaryMin(G4double) = 0 ;
    virtual G4double GetBoundaryMax(G4double) = 0 ;
    virtual G4double GetSurfaceArea() = 0 ;
    virtual void GetFacets(G4int m, G4int n, G4double xyz[][3],
                           G4int faces[][4], G4int iside) = 0 ;

    /**
     * Node/face mapping functions.
     */
    G4int GetNode( G4int i, G4int j, G4int m, G4int n, G4int iside ) ;
    G4int GetFace( G4int i, G4int j, G4int m, G4int n, G4int iside ) ;

    /**
     * Checks face visibility, based on verteces and orientation.
     *  @returns -1 if invisible; 1 if visible.
     */
    G4int GetEdgeVisibility( G4int i, G4int j, G4int m, G4int n,
                             G4int number, G4int orientation) ;

    /**
     * Returns the surface name.
     */
    inline const G4String& GetName() const { return fName; }

    /**
     * Prints on standard output surface data.
     */
    void DebugPrint() const;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4VTwistSurface(__void__&);

  protected:
 
    /**
     * Internal accessors.
     */
    inline G4VTwistSurface** GetNeighbours() { return fNeighbours; } 
    inline G4int GetNeighbours(G4int areacode, G4VTwistSurface* surfaces[]);
    inline G4ThreeVector GetCorner(G4int areacode) const;
    void GetBoundaryAxis(G4int areacode, EAxis axis[]) const;
    void GetBoundaryLimit(G4int areacode, G4double limit[]) const;
    virtual G4int GetAreaCode(const G4ThreeVector& xx, G4bool withtol=true) = 0;

    /**
     * Modifier for boundaries.
     *  @param[in] axiscode The location code, either: sAxis0 & sAxisMin,
     *             sAxis0 & sAxisMax, sAxis1 & sAxisMin, sAxis1 & sAxisMax.
     *  @param[in] direction Direction.
     *  @param[in] x0 Point on boundary.
     *  @param[in] boundarytype The shape of locus from the start point to
     *             end point of boundary (ex. sAxisRho = linear line which
     *             start point is fixed at origin; sAxisPhi = part of circle
     *             which center placed at the origin.
     */
    virtual void SetBoundary(const G4int& axiscode, 
                             const G4ThreeVector& direction,
                             const G4ThreeVector& x0, 
                             const G4int& boundarytype);

    /**
     * Modifier for corners.
     */
    void SetCorner(G4int areacode, G4double x, G4double y, G4double z);

  private:

    /**
     * Internal setters to be implemented in concrete dereived classes.
     */
    virtual void SetBoundaries() = 0;
    virtual void SetCorners()    = 0;
   
  // data members ---------------------------------------------------------

  public:

    /** Boundary types. */
    static const G4int sOutside ;
    static const G4int sInside  ;
    static const G4int sBoundary;
    static const G4int sCorner;
    static const G4int sC0Min1Min;
    static const G4int sC0Max1Min;
    static const G4int sC0Max1Max;
    static const G4int sC0Min1Max;
    static const G4int sAxisMin;
    static const G4int sAxisMax;
    static const G4int sAxisX;
    static const G4int sAxisY;
    static const G4int sAxisZ;
    static const G4int sAxisRho;
    static const G4int sAxisPhi;
    static const G4int sAxis0;
    static const G4int sAxis1;
    static const G4int sSizeMask;
    static const G4int sAxisMask;
    static const G4int sAreaMask;

  protected:
 
    class CurrentStatus 
    {
      /**
       * @brief Internal class defining the surface status.
       */
      public:

        CurrentStatus();
        virtual ~CurrentStatus();
      
        inline G4ThreeVector GetXX(G4int i) const { return fXX[i]; }
        inline G4double GetDistance(G4int i) const { return fDistance[i]; }
        inline G4int GetAreacode(G4int i) const { return fAreacode[i]; }
        inline G4int GetNXX() const { return fNXX; }
        inline G4bool IsDone() const { return fDone; }
        inline G4bool IsValid(G4int i) const { return fIsValid[i]; }

        void SetCurrentStatus(G4int  i, 
                              G4ThreeVector& xx, 
                              G4double& dist, 
                              G4int& areacode, 
                              G4bool& isvalid,
                              G4int nxx,
                              EValidate validate,
                        const G4ThreeVector* p, 
                        const G4ThreeVector* v = nullptr);

        void ResetfDone(EValidate validate,
                  const G4ThreeVector* p, 
                  const G4ThreeVector* v = nullptr);


        void DebugPrint() const;

      private:

        G4double fDistance[G4VSURFACENXX];
        G4ThreeVector fXX[G4VSURFACENXX];
        G4int fAreacode[G4VSURFACENXX];
        G4bool fIsValid[G4VSURFACENXX];
        G4int fNXX;
        G4ThreeVector fLastp;
        G4ThreeVector fLastv;
        EValidate fLastValidate;
        G4bool fDone;
    };
      
    class Boundary 
    {
      /**
       * @brief Internal class defining a surface boundary type.
       */
      public:

        Boundary() = default;
        virtual ~Boundary() = default;
      
        void SetFields(const G4int& areacode, 
                       const G4ThreeVector& d, 
                       const G4ThreeVector& x0, 
                       const G4int& boundarytype);
      
        G4bool IsEmpty() const;
      
        G4bool GetBoundaryParameters(const G4int& areacode, 
                                           G4ThreeVector& d,
                                           G4ThreeVector& x0, 
                                           G4int& boundarytype) const;

      private:

        G4int fBoundaryAcode{-1};
        G4ThreeVector fBoundaryDirection;
        G4ThreeVector fBoundaryX0;
        G4int fBoundaryType{0};
    };

    EAxis fAxis[2];
    G4double fAxisMin[2];
    G4double fAxisMax[2];
    CurrentStatus fCurStatWithV;
    CurrentStatus fCurStat;
    G4RotationMatrix fRot;
    G4ThreeVector fTrans;
    G4int fHandedness;

    class G4SurfCurNormal
    {
      public:

        G4ThreeVector p;
        G4ThreeVector normal;
    };

    G4SurfCurNormal fCurrentNormal;
    G4bool fIsValidNorm;
    G4double kCarTolerance;

  private:
                     
    G4VTwistSurface* fNeighbours[4]; // {0,1,2,3} = sAxis0min, sAxis1min, 
                                     //             sAxis0max, sAxis1max
    G4ThreeVector fCorners[4]; // corners of the surface in local coordinate
    Boundary fBoundaries[4];   // boundaries of the surface.
    G4String fName;
   
    class G4SurfSideQuery
    {
      public:

        G4ThreeVector me;
        G4ThreeVector vec;
        G4bool  withTol;
        G4int amIOnLeftSide;
    };
    G4SurfSideQuery fAmIOnLeftSide;
};

//========================================================
// inline functions
//========================================================

#include "G4VTwistSurface.icc"

#endif
