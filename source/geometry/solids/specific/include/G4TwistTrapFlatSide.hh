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
// G4TwistTrapFlatSide
//
// Class description:
//
// Class describing a flat boundary surface for a trapezoid.

// Author: Oliver Link (CERN), 27.10.2004 - Created
// --------------------------------------------------------------------
#ifndef G4TWISTTRAPFLATSIDE_HH
#define G4TWISTTRAPFLATSIDE_HH

#include "G4VTwistSurface.hh"

/**
 * @brief G4TwistTrapFlatSide describes a flat boundary surface for a trapezoid.
 */

class G4TwistTrapFlatSide : public G4VTwistSurface
{
  public:

    /**
     * Constructs a trapezoid flat boundary surface, given its parameters.
     *  @param[in] name The surface name.
     *  @param[in] PhiTwist The twist angle.
     *  @param[in] pDx1 Half x length at -pDz,-pDy.
     *  @param[in] pDx2 Half x length at -pDz,+pDy.
     *  @param[in] pDy Half y length.
     *  @param[in] pDz Half z length.
     *  @param[in] pAlpha Tilt angle at +pDz.
     *  @param[in] pPhi Direction between end planes - azimuthal angle.
     *  @param[in] pTheta Direction between end planes - polar angle.
     *  @param[in] handedness Orientation: +z = +ve, -z = -ve.
     */
    G4TwistTrapFlatSide( const G4String& name,
                               G4double  PhiTwist,
                               G4double  pDx1,
                               G4double  pDx2,
                               G4double  pDy,
                               G4double  pDz,
                               G4double  pAlpha,
                               G4double  pPhi,
                               G4double  pTheta,
                               G4int     handedness  );

    /**
     * Default destructor.
     */
    ~G4TwistTrapFlatSide() override = default;

    /**
     * Returns a normal vector at a surface (or very close to the surface)
     * point at 'p'.
     *  @param[in] p Not used. Using current normal.
     *  @param[in] isGlobal If true, it returns the normal in global coordinates.
     *  @returns The current normal vector.
     */
    G4ThreeVector GetNormal(const G4ThreeVector& /* p */ ,
                                  G4bool isGlobal = false) override;

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
    G4int DistanceToSurface(const G4ThreeVector& gp,
                            const G4ThreeVector& gv,
                                  G4ThreeVector  gxx[],
                                  G4double       distance[],
                                  G4int          areacode[],
                                  G4bool         isvalid[],
                            EValidate validate = kValidateWithTol) override;

    /**
     * Returns the safety distance to surface, given point 'gp'.
     *  @param[in] gp The point from where computing the safety distance.
     *  @param[out] gxx Vector of global points based on number of solutions.
     *  @param[out] distance The distance vector based on number of solutions.
     *  @param[out] areacode The location vector based on number of solutions.
     *  @returns The number of solutions.
     */
    G4int DistanceToSurface(const G4ThreeVector& gp,
                                  G4ThreeVector  gxx[],
                                  G4double       distance[],
                                  G4int          areacode[]) override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4TwistTrapFlatSide(__void__&);

  private:

    /**
     * Returns point on surface given 'phi' and 'u'.
     */
    inline G4ThreeVector SurfacePoint(G4double x, G4double y,
                                      G4bool isGlobal = false) override;  

    /**
     * Internal accessors.
     */
    inline G4double GetBoundaryMin(G4double u) override;
    inline G4double GetBoundaryMax(G4double u) override;
    inline G4double GetSurfaceArea() override;
    void GetFacets( G4int m, G4int n, G4double xyz[][3],
                    G4int faces[][4], G4int iside ) override;
    inline G4double xAxisMax(G4double u, G4double fTanAlpha) const;

    /**
     * Setters.
     */
    void SetCorners() override;
    void SetBoundaries() override;

    /**
     * Returns the area code for point 'xx' using or not surface tolerance.
     */
    G4int GetAreaCode(const G4ThreeVector& xx, 
                            G4bool withTol = true) override;

  private:
  
    G4double fDx1;
    G4double fDx2;
    G4double fDy;
    G4double fDz;
    G4double fPhiTwist;
    G4double fAlpha;
    G4double fTAlph;
    G4double fPhi;
    G4double fTheta;
    G4double fdeltaX;
    G4double fdeltaY;
};

//========================================================
// inline functions
//========================================================

inline 
G4double G4TwistTrapFlatSide::xAxisMax(G4double u, G4double fTanAlpha) const
{
  return (  ( fDx2 + fDx1 )/2. + u*(fDx2 - fDx1)/(2.*fDy) + u *fTanAlpha  ) ;
}

inline G4ThreeVector
G4TwistTrapFlatSide::SurfacePoint(G4double x, G4double y, G4bool isGlobal)
{
  G4ThreeVector SurfPoint ( x,y,0);

  if (isGlobal) { return (fRot*SurfPoint + fTrans); }
  return SurfPoint;
}

inline
G4double G4TwistTrapFlatSide::GetBoundaryMin(G4double y)
{
  return -xAxisMax(y, -fTAlph) ;
}

inline
G4double G4TwistTrapFlatSide::GetBoundaryMax(G4double y)
{
  return xAxisMax(y, fTAlph) ; 
}

inline
G4double G4TwistTrapFlatSide::GetSurfaceArea()
{
  return 2*(fDx1 + fDx2)*fDy ;
}

#endif
