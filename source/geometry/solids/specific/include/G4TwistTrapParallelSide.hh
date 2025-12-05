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
// G4TwistTrapParallelSide
//
// Class description:
//
// Class describing a twisted boundary surface for a trapezoid.

// Author: Oliver Link (CERN), 27.10.2004 - Created
// --------------------------------------------------------------------
#ifndef G4TWISTTRAPPARALLELSIDE_HH
#define G4TWISTTRAPPARALLELSIDE_HH

#include "G4VTwistSurface.hh"

#include <vector>

/**
 * @brief G4TwistTrapParallelSide describes a twisted boundary surface for
 * a trapezoid.
 */

class G4TwistTrapParallelSide : public G4VTwistSurface
{
  public:
   
    /**
     * Constructs a trapezoid twisted boundary surface, given its parameters.
     *  @param[in] name The surface name.
     *  @param[in] PhiTwist The twist angle.
     *  @param[in] pDz Half z length.
     *  @param[in] pTheta Direction between end planes - polar angle.
     *  @param[in] pPhi Direction between end planes - azimuthal angle.
     *  @param[in] pDy1 Half y length at -pDz.
     *  @param[in] pDx1 Half x length at -pDz,-pDy.
     *  @param[in] pDx2 Half x length at -pDz,+pDy.
     *  @param[in] pDy2 Half y length at +pDz.
     *  @param[in] pDx3 Half x length at +pDz,-pDy.
     *  @param[in] pDx4 Half x length at +pDz,+pDy.
     *  @param[in] pAlph Tilt angle at +pDz.
     *  @param[in] AngleSide Parity.
     */
    G4TwistTrapParallelSide(const G4String& name,
                              G4double  PhiTwist, // twist angle
                              G4double  pDz,      // half z lenght
                              G4double  pTheta, // direction between end planes
                              G4double  pPhi,   // by polar and azimutal angles
                              G4double  pDy1,     // half y length at -pDz
                              G4double  pDx1,     // half x length at -pDz,-pDy
                              G4double  pDx2,     // half x length at -pDz,+pDy
                              G4double  pDy2,     // half y length at +pDz
                              G4double  pDx3,     // half x length at +pDz,-pDy
                              G4double  pDx4,     // half x length at +pDz,+pDy
                              G4double  pAlph,    // tilt angle at +pDz
                              G4double  AngleSide // parity
                            );
  
    /**
     * Default destructor.
     */
    ~G4TwistTrapParallelSide() override = default;
   
    /**
     * Returns a normal vector at a surface (or very close to the surface)
     * point at 'p'.
     *  @param[in] p The point where computing the normal.
     *  @param[in] isGlobal If true, it returns the normal in global coordinates.
     *  @returns The normal vector.
     */
    G4ThreeVector GetNormal(const G4ThreeVector& p,
                                  G4bool isGlobal = false) override ;   
   
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
                                  G4double  distance[],
                                  G4int     areacode[],
                                  G4bool    isvalid[],
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
    G4TwistTrapParallelSide(__void__&);

  private:

    /**
     * Returns the area code for point 'xx' using or not surface tolerance.
     */
    G4int GetAreaCode(const G4ThreeVector& xx, 
                            G4bool withTol = true) override;

    /**
     * Setters.
     */
    void SetCorners() override;
    void SetBoundaries() override;

    /**
     * Finds the closest point on surface for a given point 'p', returning
     * 'phi' and 'u'.
     */
    void GetPhiUAtX(const G4ThreeVector& p, G4double& phi, G4double& u);

    /**
     * Returns projection on surface of a given point 'p'.
     */
    G4ThreeVector ProjectPoint(const G4ThreeVector& p,
                                     G4bool isglobal = false);

    /**
     * Returns point on surface given 'phi' and 'u'.
     */
    inline G4ThreeVector SurfacePoint(G4double phi, G4double u,
                                      G4bool isGlobal = false) override;

    /**
     * Internal accessors.
     */
    inline G4double GetBoundaryMin(G4double phi) override;
    inline G4double GetBoundaryMax(G4double phi) override;
    inline G4double GetSurfaceArea() override;
    void GetFacets( G4int m, G4int n, G4double xyz[][3],
                    G4int faces[][4], G4int iside ) override;

    inline G4ThreeVector NormAng(G4double phi, G4double u);
    inline G4double GetValueB(G4double phi) ;
    inline G4double Xcoef(G4double u);
      // To calculate the w(u) function

  private:

    G4double fTheta;   
    G4double fPhi ;

    G4double fDy1;   
    G4double fDx1;     
    G4double fDx2;     

    G4double fDy2;   
    G4double fDx3;     
    G4double fDx4;     

    G4double fDz;         // Half-length along the z axis

    G4double fAlph;
    G4double fTAlph;      // std::tan(fAlph)
    
    G4double fPhiTwist;   // twist angle ( dphi in surface equation)

    G4double fAngleSide;

    G4double fdeltaX;
    G4double fdeltaY;

    G4double fDx4plus2;  // fDx4 + fDx2  == a2/2 + a1/2
    G4double fDx4minus2; // fDx4 - fDx2          -
    G4double fDx3plus1;  // fDx3 + fDx1  == d2/2 + d1/2
    G4double fDx3minus1; // fDx3 - fDx1          -
    G4double fDy2plus1;  // fDy2 + fDy1  == b2/2 + b1/2
    G4double fDy2minus1; // fDy2 - fDy1          -
    G4double fa1md1;     // 2 fDx2 - 2 fDx1  == a1 - d1
    G4double fa2md2;     // 2 fDx4 - 2 fDx3 
};   

//========================================================
// inline functions
//========================================================

#include "G4TwistTrapParallelSide.icc"

#endif
