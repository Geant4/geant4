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
// G4TwistTubsHypeSide
//
// Class description:
//
// Class describing a hyperbolic boundary surface for a cylinder.

// Author: Kotoyo Hoshina (Chiba University), 01.08.2002 - Created.
//         Oliver Link (CERN), 13.11.2003 - Integration in Geant4
//               from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------
#ifndef G4TWISTTUBSHYPESIDE_HH
#define G4TWISTTUBSHYPESIDE_HH

#include "G4VTwistSurface.hh"
#include "G4Integrator.hh"
#include "G4SimpleIntegration.hh"

/**
 * @brief G4TwistTubsHypeSide describes hyperbolic boundary surface
 * for a cylinder.
 */

class G4TwistTubsHypeSide : public G4VTwistSurface
{
  public:
                       
    /**
     * Constructs a cylinder hyperbolic boundary surface, given its parameters.
     *  @param[in] name The surface name.
     *  @param[in] rot Rotation: 0.5*(phi-width segment).
     *  @param[in] tlate Translation.
     *  @param[in] handedness Orientation: R-hand = 1, L-hand = -1.
     *  @param[in] kappa Kappa=tan(TwistAngle/2)/fZHalfLen.
     *  @param[in] tanstereo Tangent of the stereo angle.
     *  @param[in] r0 Radius at z = 0.
     *  @param[in] axis0 Phi axis.
     *  @param[in] axis1 Z axis.
     *  @param[in] axis0min Minimum in Phi.
     *  @param[in] axis1min Minimum in Z.
     *  @param[in] axis0max Maximum in Phi.
     *  @param[in] axis1max Maximum in Z.
     */
    G4TwistTubsHypeSide(const G4String& name,
                        const G4RotationMatrix& rot,  // 0.5*(phi-width segment)
                        const G4ThreeVector& tlate,
                        const G4int handedness,   // R-hand = 1, L-hand = -1
                        const G4double kappa,     // tan(TwistAngle/2)/fZHalfLen
                        const G4double tanstereo, // tan(stereo angle)
                        const G4double r0,        // radius at z = 0
                        const EAxis axis0 = kPhi,
                        const EAxis axis1 = kZAxis,
                              G4double axis0min = -kInfinity,
                              G4double axis1min = -kInfinity,
                              G4double axis0max = kInfinity,
                              G4double axis1max = kInfinity); 
                             
    /**
     * Alternative Construct for a cylinder hyperbolic boundary surface.
     *  @param[in] name The surface name.
     *  @param[in] EndInnerRadius Inner-hype radius at z=0.
     *  @param[in] EndOuterRadius Outer-hype radius at z=0.
     *  @param[in] DPhi Phi angle.
     *  @param[in] EndPhi Total Phi.
     *  @param[in] EndZ Z length.
     *  @param[in] InnerRadius Inner radius.
     *  @param[in] OuterRadius Outer radius.
     *  @param[in] Kappa Kappa=tan(TwistAngle/2)/fZHalfLen.
     *  @param[in] TanInnerStereo Tangent inner stereo angle.
     *  @param[in] TanOuterStereo Tangent outer stereo angle.
     *  @param[in] handedness Orientation: R-hand = 1, L-hand = -1.
     */
   G4TwistTubsHypeSide(const G4String&  name,
                             G4double   EndInnerRadius[2],
                             G4double   EndOuterRadius[2],
                             G4double   DPhi,
                             G4double   EndPhi[2],
                             G4double   EndZ[2], 
                             G4double   InnerRadius,
                             G4double   OuterRadius,
                             G4double   Kappa,
                             G4double   TanInnerStereo,
                             G4double   TanOuterStereo,
                             G4int      handedness) ;

    /**
     * Default destructor.
     */
    ~G4TwistTubsHypeSide() override = default;

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
     * Returns a normal vector at a surface (or very close to the surface)
     * point at 'p'.
     *  @param[in] p The point where computing the normal.
     *  @param[in] isGlobal If true, it returns the normal in global coordinates.
     *  @returns The normal vector.
     */
    G4ThreeVector GetNormal(const G4ThreeVector& p,
                                  G4bool isGlobal = false) override ;

    /**
     * Returns if point at 'gp' is inside surface.
     */
    EInside Inside(const G4ThreeVector& gp) ;

    /**
     * Gets Rho at p.z() on Hyperbolic Surface.
     */
    inline G4double GetRhoAtPZ(const G4ThreeVector& p,
                                     G4bool isglobal = false) const ;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4TwistTubsHypeSide(__void__&);

  private:

    /**
     * Returns point on surface given 'phi' and 'z'.
     */
    inline G4ThreeVector SurfacePoint(G4double phi, G4double z,
                                      G4bool isGlobal = false) override ;  

    /**
     * Internal accessors.
     */
    inline G4double GetBoundaryMin(G4double phi) override ;
    inline G4double GetBoundaryMax(G4double phi) override ;
    inline G4double GetSurfaceArea() override ;
    void GetFacets( G4int m, G4int n, G4double xyz[][3],
                    G4int faces[][4], G4int iside ) override ;
    
    /**
     * Returns the area code for point 'xx' using or not surface tolerance.
     */
    G4int GetAreaCode(const G4ThreeVector& xx, 
                            G4bool withTol = true) override;
    G4int GetAreaCodeInPhi(const G4ThreeVector& xx, 
                                 G4bool withTol = true);

    /**
     * Setters.
     */
    void SetCorners() override;
    void SetCorners(G4double EndInnerRadius[2],
                    G4double EndOuterRadius[2],
                    G4double DPhi,
                    G4double EndPhi[2],
                    G4double EndZ[2]);
    void SetBoundaries() override;

  private:
   
    G4double fKappa;        // std::tan(TwistedAngle/2)/HalfLenZ;
    G4double fTanStereo;    // std::tan(StereoAngle)
    G4double fTan2Stereo;   // std::tan(StereoAngle)**2
    G4double fR0;           // radius at z = 0
    G4double fR02;          // radius**2 at z = 0
    G4double fDPhi ;        // segment

    class Insidetype
    {
      public:

        G4ThreeVector gp;
        EInside       inside;
    };
    Insidetype fInside;
};

//========================================================
// inline functions
//========================================================

inline
G4double G4TwistTubsHypeSide::GetRhoAtPZ(const G4ThreeVector& p,
                                               G4bool isglobal) const 
{
  // Get Rho at p.z() on Hyperbolic Surface.
  G4ThreeVector tmpp;
  if (isglobal) { tmpp = fRot.inverse()*p - fTrans; }
  else          { tmpp = p; }

  return std::sqrt(fR02 + tmpp.z() * tmpp.z() * fTan2Stereo); 
}

inline
G4ThreeVector G4TwistTubsHypeSide::
SurfacePoint(G4double phi , G4double z , G4bool isGlobal)
{
  G4double rho = std::sqrt(fR02 + z * z * fTan2Stereo) ;

  G4ThreeVector SurfPoint (rho*std::cos(phi), rho*std::sin(phi), z) ;

  if (isGlobal) { return (fRot * SurfPoint + fTrans); }
  return SurfPoint;
}

inline
G4double G4TwistTubsHypeSide::GetBoundaryMin(G4double z)
{
  G4ThreeVector ptmp(0,0,z) ;  // temporary point with z Komponent only
  G4ThreeVector lowerlimit;    // lower phi-boundary limit at z = ptmp.z()
  lowerlimit = GetBoundaryAtPZ(sAxis0 & sAxisMin, ptmp);
  return  std::atan2( lowerlimit.y(), lowerlimit.x() ) ;  
}

inline
G4double G4TwistTubsHypeSide::GetBoundaryMax(G4double z )
{
  G4ThreeVector ptmp(0,0,z) ;  // temporary point with z Komponent only
  G4ThreeVector upperlimit;    // upper phi-boundary limit at z = ptmp.z()
  upperlimit = GetBoundaryAtPZ(sAxis0 & sAxisMax, ptmp);
  return   std::atan2( upperlimit.y(), upperlimit.x() ) ;
}

inline
G4double G4TwistTubsHypeSide::GetSurfaceArea()
{
  // approximation with tube surface

  return ( fAxisMax[1] - fAxisMin[1] ) * fR0 * fDPhi ;
}

#endif
