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
// G4TwistTubsFlatSide
//
// Class description:
//
// Class describing a flat boundary surface for a cylinder.

// Author: Kotoyo Hoshina (Chiba University), 01.08.2002 - Created.
//         Oliver Link (CERN), 13.11.2003 - Integration in Geant4
//               from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------
#ifndef G4TWISTTUBSFLATSIDE_HH
#define G4TWISTTUBSFLATSIDE_HH

#include "G4VTwistSurface.hh"

/**
 * @brief G4TwistTubsFlatSide describes a flat boundary surface for a cylinder.
 */

class G4TwistTubsFlatSide : public G4VTwistSurface
{
  public:

    /**
     * Constructs a cylinder flat boundary surface, given its parameters.
     *  @param[in] name The surface name.
     *  @param[in] rot Rotation.
     *  @param[in] tlate Translation.
     *  @param[in] n Normal vector.
     *  @param[in] axis0 Rho axis.
     *  @param[in] axis1 Phi axis.
     *  @param[in] axis0min Minimum in Rho.
     *  @param[in] axis1min Minimum in Phi.
     *  @param[in] axis0max Maximum in Rho.
     *  @param[in] axis1max Maximum in Phi.
     */
    G4TwistTubsFlatSide(const G4String& name,
                        const G4RotationMatrix& rot,
                        const G4ThreeVector& tlate,
                        const G4ThreeVector& n,
                        const EAxis axis0 = kRho, // RHO axis !
                        const EAxis axis1 = kPhi, // PHI axis !
                              G4double axis0min = -kInfinity,
                              G4double axis1min = -kInfinity,
                              G4double axis0max = kInfinity,
                              G4double axis1max = kInfinity);

    /**
     * Alternative Construct for a cylinder flat boundary surface.
     *  @param[in] name The surface name.
     *  @param[in] EndInnerRadius Inner-hype radius at z=0.
     *  @param[in] EndOuterRadius Outer-hype radius at z=0.
     *  @param[in] DPhi Phi angle.
     *  @param[in] EndPhi Total Phi.
     *  @param[in] EndZ Z length.
     *  @param[in] handedness Orientation: +z = +ve, -z = -ve.
     */
    G4TwistTubsFlatSide(const G4String& name,
                              G4double EndInnerRadius[2],
                              G4double EndOuterRadius[2],
                              G4double DPhi,
                              G4double EndPhi[2],
                              G4double EndZ[2], 
                              G4int handedness);

    /**
     * Default destructor.
     */
    ~G4TwistTubsFlatSide() override = default;

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
    G4TwistTubsFlatSide(__void__&);

  private:

    /**
     * Returns point on surface given 'phi' and 'u'.
     */
    inline G4ThreeVector SurfacePoint(G4double, G4double,
                                      G4bool isGlobal = false) override;  

    /**
     * Internal accessors.
     */
    inline G4double GetBoundaryMin(G4double phi) override;
    inline G4double GetBoundaryMax(G4double phi) override;
    inline G4double GetSurfaceArea() override { return fSurfaceArea ; }
    void GetFacets( G4int m, G4int n, G4double xyz[][3],
                    G4int faces[][4], G4int iside ) override;

    /**
     * Returns the area code for point 'xx' using or not surface tolerance.
     */
    G4int GetAreaCode(const G4ThreeVector& xx, 
                            G4bool withTol = true) override ;

    /**
     * Setters.
     */
    void SetCorners() override;
    void SetBoundaries() override;

  private:

    G4double fSurfaceArea = 0.0;
};

//========================================================
// inline functions
//========================================================

inline G4ThreeVector G4TwistTubsFlatSide::
SurfacePoint(G4double phi , G4double rho , G4bool isGlobal )
{
  G4ThreeVector SurfPoint (rho*std::cos(phi) , rho*std::sin(phi) , 0);

  if (isGlobal) { return (fRot * SurfPoint + fTrans); }
  return SurfPoint;
}

inline
G4double G4TwistTubsFlatSide::GetBoundaryMin(G4double)
{
  G4ThreeVector dphimin = GetCorner(sC0Max1Min);
  return  std::atan2( dphimin.y(), dphimin.x() );  
}

inline
G4double G4TwistTubsFlatSide::GetBoundaryMax(G4double)
{
  G4ThreeVector dphimax = GetCorner(sC0Max1Max);   
  return  std::atan2( dphimax.y(), dphimax.x() );  
}

#endif
