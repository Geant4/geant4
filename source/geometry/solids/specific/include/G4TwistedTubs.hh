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
// G4TwistedTubs
//
// Class description:
//
// G4TwistedTubs is a sector of a twisted hollow cylinder.
// A twisted cylinder which is placed along with z-axis and is
// separated into phi-segments should become a hyperboloid, and
// its each segmented piece should be tilted with a stereo angle.
// G4TwistedTubs is a G4VSolid.
//
// Details of the implementation: "Development of a Geant4 solid
// for stereo mini-jet cells in a cylindrical drift chamber",
// Computer Physics Communications 153 (2003) pp.373-391

// Author: Kotoyo Hoshina (Chiba University), 01.08.2002 - created.
//         Oliver Link (CERN), 13.11.2003 - Integration in Geant4
//                from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------
#ifndef G4TWISTEDTUBS_HH
#define G4TWISTEDTUBS_HH

#include "G4VSolid.hh"
#include "G4TwistTubsFlatSide.hh"
#include "G4TwistTubsSide.hh"
#include "G4TwistTubsHypeSide.hh"

class G4SolidExtentList;
class G4ClippablePolygon;

/**
 * @brief G4TwistedTubs is a sector of a twisted hollow cylinder.
 * A twisted cylinder which is placed along with Z axis and is separated into
 * phi-segments should become a hyperboloid, and its each segmented piece
 * should be tilted with a stereo angle.
 */

class G4TwistedTubs : public G4VSolid
{
  public:
 
    /**
     * Constructs a twisted tube, given radii and twist angle.
     *  @param[in] pName The solid name.
     *  @param[in] twistedangle The twisted angle.
     *  @param[in] endinnerrad Inner radius at endcap.
     *  @param[in] endouterrad Outer radius at endcap.
     *  @param[in] halfzlen Half Z length.
     *  @param[in] dphi Phi angle of a segment.
     */
    G4TwistedTubs(const G4String& pname,        // Name of instance
                        G4double  twistedangle, // Twisted angle
                        G4double  endinnerrad,  // Inner radius at endcap
                        G4double  endouterrad,  // Outer radius at endcap
                        G4double  halfzlen,     // half z length
                        G4double  dphi);        // Phi angle of a segment
                      
    /**
     * Constructs a twisted tube, given radii, twist angle and number of segments.
     *  @param[in] pName The solid name.
     *  @param[in] twistedangle The twisted angle.
     *  @param[in] endinnerrad Inner radius at endcap.
     *  @param[in] endouterrad Outer radius at endcap.
     *  @param[in] halfzlen Half Z length.
     *  @param[in] nseg Number of segments in Phi.
     *  @param[in] totphi Total angle of all segments.
     */
    G4TwistedTubs(const G4String& pname,        // Name of instance
                        G4double  twistedangle, // Stereo angle
                        G4double  endinnerrad,  // Inner radius at endcap
                        G4double  endouterrad,  // Outer radius at endcap
                        G4double  halfzlen,     // half z length 
                        G4int     nseg,         // Number of segments in totalPhi
                        G4double  totphi);      // Total angle of all segments
                      
    /**
     * Constructs a twisted tube, given radii, twist angle and EndZ values.
     *  @param[in] pName The solid name.
     *  @param[in] twistedangle The twisted angle.
     *  @param[in] innerrad Inner radius at z=0.
     *  @param[in] outerrad Outer radius at z=0.
     *  @param[in] negativeEndz -ve Z endplate.
     *  @param[in] positiveEndz +ve Z endplate.
     *  @param[in] dphi Phi angle of a segment.
     */
    G4TwistedTubs(const G4String& pname,        // Name of instance
                        G4double  twistedangle, // Twisted angle
                        G4double  innerrad,     // Inner radius at z=0
                        G4double  outerrad,     // Outer radius at z=0
                        G4double  negativeEndz, // -ve z endplate
                        G4double  positiveEndz, // +ve z endplate
                        G4double  dphi);        // Phi angle of a segment

    /**
     * Constructs a twisted tube, given radii, twist angle, EndZ values
     * and number of segments.
     *  @param[in] pName The solid name.
     *  @param[in] twistedangle The twisted angle.
     *  @param[in] innerrad Inner radius at z=0.
     *  @param[in] outerrad Outer radius at z=0.
     *  @param[in] negativeEndz -ve Z endplate.
     *  @param[in] positiveEndz +ve Z endplate.
     *  @param[in] nseg Number of segments in Phi.
     *  @param[in] totphi Total angle of all segments.
     */
    G4TwistedTubs(const G4String& pname,        // Name of instance
                        G4double  twistedangle, // Stereo angle
                        G4double  innerrad,     // Inner radius at z=0 
                        G4double  outerrad,     // Outer radius at z=0 
                        G4double  negativeEndz, // -ve z endplate
                        G4double  positiveEndz, // +ve z endplate
                        G4int     nseg,         // Number of segments in totalPhi
                        G4double  totphi);      // Total angle of all segments

    /**
     * Destructor.
     */
    ~G4TwistedTubs() override;
             
    /**
     * Thows a fatal exception. Parameterisations are not allowed on this shape.
     */
    void ComputeDimensions(G4VPVParameterisation*    /* p  */ ,
                           const G4int               /* n  */ ,
                           const G4VPhysicalVolume*  /* prep */ ) override;

    /**
     * Computes the bounding limits of the solid.
     *  @param[out] pMin The minimum bounding limit point.
     *  @param[out] pMax The maximum bounding limit point.
     */
    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;

    /**
     * Calculates the minimum and maximum extent of the solid, when under the
     * specified transform, and within the specified limits.
     *  @param[in] pAxis The axis along which compute the extent.
     *  @param[in] pVoxelLimit The limiting space dictated by voxels.
     *  @param[in] pTransform The internal transformation applied to the solid.
     *  @param[out] pMin The minimum extent value.
     *  @param[out] pMax The maximum extent value.
     *  @returns True if the solid is intersected by the extent region.
     */
    G4bool CalculateExtent(const EAxis               pAxis,
                           const G4VoxelLimits&     pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double&          pMin,
                                 G4double&          pMax ) const override;

    /**
     * Concrete implementations of the expected query interfaces for
     * solids, as defined in the base class G4VSolid.
     */
    G4double DistanceToIn (const G4ThreeVector& p,
                           const G4ThreeVector& v ) const override;
    G4double DistanceToIn (const G4ThreeVector& p ) const override;
    G4double DistanceToOut(const G4ThreeVector& p, 
                           const G4ThreeVector& v,
                           const G4bool calcnorm = false,
                                 G4bool* validnorm = nullptr, 
                                 G4ThreeVector* n = nullptr ) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;
    EInside Inside (const G4ThreeVector& p) const override;
    G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo (G4VGraphicsScene& scene) const override;
    G4Polyhedron* CreatePolyhedron () const override;
    G4Polyhedron* GetPolyhedron () const override;
    G4VisExtent GetExtent () const override;

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream &StreamInfo(std::ostream& os) const override;

    /**
     * Accessors.
     */
    inline G4double GetDPhi        () const { return fDPhi       ; }
    inline G4double GetPhiTwist    () const { return fPhiTwist   ; }
    inline G4double GetInnerRadius () const { return fInnerRadius; }
    inline G4double GetOuterRadius () const { return fOuterRadius; }
    inline G4double GetInnerStereo () const { return fInnerStereo; }
    inline G4double GetOuterStereo () const { return fOuterStereo; }
    inline G4double GetZHalfLength () const { return fZHalfLength; }
    inline G4double GetKappa       () const { return fKappa      ; }
    //
    inline G4double GetTanInnerStereo () const { return fTanInnerStereo  ; }
    inline G4double GetTanInnerStereo2() const { return fTanInnerStereo2 ; }
    inline G4double GetTanOuterStereo () const { return fTanOuterStereo  ; }
    inline G4double GetTanOuterStereo2() const { return fTanOuterStereo2 ; }
    //
    inline G4double GetEndZ           (G4int i) const { return fEndZ[i]  ; }
    inline G4double GetEndPhi         (G4int i) const { return fEndPhi[i]; }
    inline G4double GetEndInnerRadius (G4int i) const
                    { return fEndInnerRadius[i]; }
    inline G4double GetEndOuterRadius (G4int i) const
                    { return fEndOuterRadius[i]; }
    inline G4double GetEndInnerRadius () const 
                    { return (fEndInnerRadius[0] > fEndInnerRadius[1] ?
                      fEndInnerRadius[0] : fEndInnerRadius[1]); }
    inline G4double GetEndOuterRadius () const
                    { return (fEndOuterRadius[0] > fEndOuterRadius[1] ?
                      fEndOuterRadius[0] : fEndOuterRadius[1]); }
  
    /**
     * Returns the type ID, "G4TwistedTubs" of the solid.
     */
    G4GeometryType  GetEntityType() const override;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Returning an estimation of the solid volume (capacity) and surface area,
     * in internal units. Caches the computed value once computed the first time.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    /**
     * Returns a random point located and uniformly distributed on the
     * surface of the solid.
     */
    G4ThreeVector GetPointOnSurface() const override ;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4TwistedTubs(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4TwistedTubs(const G4TwistedTubs& rhs);
    G4TwistedTubs& operator=(const G4TwistedTubs& rhs); 

#ifdef G4TWISTDEBUG
    inline G4VTwistSurface* GetOuterHype() const { return fOuterHype; }
#endif
  
  private:
 
    /**
     * Initialiser of data.
     */
    inline void SetFields(G4double phitwist, G4double innerrad,
                          G4double outerrad,
                          G4double negativeEndz, G4double positiveEndz);

    /**
     * Generates the surfaces. Used in constructors.
     */
    void CreateSurfaces();

    /**
     * Utilities for area computation.
     */
    G4double GetLateralArea(G4double a, G4double r, G4double z) const;
    G4double GetPhiCutArea(G4double a, G4double r, G4double z) const;

  private:
 
    G4double fPhiTwist;       // Twist angle from -fZHalfLength to fZHalfLength
    G4double fInnerRadius;    // Inner-hype radius at z=0
    G4double fOuterRadius;    // Outer-hype radius at z=0
    G4double fEndZ[2];        // z at endcaps, [0] = -ve z, [1] = +ve z
    G4double fDPhi;           // Phi-width of a segment fDPhi > 0
    G4double fZHalfLength;    // Half length along z-axis
 
    G4double fInnerStereo;       // Inner-hype stereo angle
    G4double fOuterStereo;       // Outer-hype stereo angle
    G4double fTanInnerStereo;    // std::tan(innerStereoAngle)
    G4double fTanOuterStereo;    // std::tan(outerStereoAngle)
    G4double fKappa;             // std::tan(fPhiTwist/2)/fZHalfLen;
    G4double fEndInnerRadius[2]; // Inner-hype radii endcaps [0] -ve z, [1] +ve z
    G4double fEndOuterRadius[2]; // Outer-hype radii endcaps [0] -ve z, [1] +ve z
    G4double fEndPhi[2];         // Phi endcaps, [0] = -ve z, [1] = +ve z
  
    G4double fInnerRadius2;      // fInnerRadius * fInnerRadius
    G4double fOuterRadius2;      // fOuterRadius * fOuterRadius
    G4double fTanInnerStereo2;   // fInnerRadius * fInnerRadius
    G4double fTanOuterStereo2;   // fInnerRadius * fInnerRadius
    G4double fEndZ2[2];          // fEndZ * fEndZ
  
    G4VTwistSurface* fLowerEndcap;    // Surface of -ve z
    G4VTwistSurface* fUpperEndcap;    // Surface of +ve z
    G4VTwistSurface* fLatterTwisted;  // Surface of -ve phi
    G4VTwistSurface* fFormerTwisted;  // Surface of +ve phi
    G4VTwistSurface* fInnerHype;      // Surface of -ve r
    G4VTwistSurface* fOuterHype;      // Surface of +ve r

    G4double fCubicVolume = 0.0;      // Cached value for cubic volume
    G4double fSurfaceArea = 0.0;      // Cached value for surface area

    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr; // polyhedron for vis

};

//=====================================================================

//---------------------
// inline functions
//---------------------

inline
void G4TwistedTubs::SetFields(G4double phitwist, G4double innerrad, 
                              G4double outerrad, G4double negativeEndz, 
                              G4double positiveEndz)
{
   fCubicVolume  = 0.;
   fPhiTwist     = phitwist;
   fEndZ[0]      = negativeEndz;
   fEndZ[1]      = positiveEndz;
   fEndZ2[0]     = fEndZ[0] * fEndZ[0];
   fEndZ2[1]     = fEndZ[1] * fEndZ[1];
   fInnerRadius  = innerrad;
   fOuterRadius  = outerrad;
   fInnerRadius2 = fInnerRadius * fInnerRadius;
   fOuterRadius2 = fOuterRadius * fOuterRadius;
   
   if (std::fabs(fEndZ[0]) >= std::fabs(fEndZ[1]))
   {
      fZHalfLength = std::fabs(fEndZ[0]);
   }
   else
   {
      fZHalfLength = std::fabs(fEndZ[1]);
   }

   G4double parity         = (fPhiTwist > 0 ? 1 : -1); 
   G4double tanHalfTwist   = std::tan(0.5 * fPhiTwist);
   G4double innerNumerator = std::fabs(fInnerRadius * tanHalfTwist) * parity;
   G4double outerNumerator = std::fabs(fOuterRadius * tanHalfTwist) * parity;

   fTanInnerStereo    = innerNumerator / fZHalfLength; 
   fTanOuterStereo    = outerNumerator / fZHalfLength; 
   fTanInnerStereo2   = fTanInnerStereo * fTanInnerStereo;
   fTanOuterStereo2   = fTanOuterStereo * fTanOuterStereo;
   fInnerStereo       = std::atan2(innerNumerator,  fZHalfLength); 
   fOuterStereo       = std::atan2(outerNumerator,  fZHalfLength); 
   fEndInnerRadius[0] = std::sqrt(fInnerRadius2 + fEndZ2[0] * fTanInnerStereo2);
   fEndInnerRadius[1] = std::sqrt(fInnerRadius2 + fEndZ2[1] * fTanInnerStereo2);
   fEndOuterRadius[0] = std::sqrt(fOuterRadius2 + fEndZ2[0] * fTanOuterStereo2);
   fEndOuterRadius[1] = std::sqrt(fOuterRadius2 + fEndZ2[1] * fTanOuterStereo2);

   fKappa          = tanHalfTwist / fZHalfLength;
   fEndPhi[0]      = std::atan2(fEndZ[0] * tanHalfTwist, fZHalfLength);
   fEndPhi[1]      = std::atan2(fEndZ[1] * tanHalfTwist, fZHalfLength);

#ifdef G4TWISTDEBUG
   G4cout << "/********* G4TwistedTubs::SetFields() Field Parameters ***************** " << G4endl;
   G4cout << "/*   fPhiTwist                  : " << fPhiTwist << G4endl;
   G4cout << "/*   fEndZ(0, 1)                : " << fEndZ[0] << " , " << fEndZ[1] << G4endl; 
   G4cout << "/*   fEndPhi(0, 1)              : " << fEndPhi[0] << " , " << fEndPhi[1] << G4endl; 
   G4cout << "/*   fInnerRadius, fOuterRadius : " << fInnerRadius << " , " << fOuterRadius << G4endl; 
   G4cout << "/*   fEndInnerRadius(0, 1)      : " << fEndInnerRadius[0] << " , " 
          << fEndInnerRadius[1] << G4endl; 
   G4cout << "/*   fEndOuterRadius(0, 1)      : " << fEndOuterRadius[0] << " , " 
          << fEndOuterRadius[1] << G4endl; 
   G4cout << "/*   fInnerStereo, fOuterStereo : " << fInnerStereo << " , " << fOuterStereo << G4endl; 
   G4cout << "/*   tanHalfTwist, fKappa       : " << tanHalfTwist << " , " << fKappa << G4endl; 
   G4cout << "/*********************************************************************** " << G4endl;
#endif
}

#endif
