//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Orb.hh,v 1.4 2003/11/05 10:56:58 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4Orb
//
// Class description:
//
//   A G4Orb is a simple case of G4Sphere. It has only:
//   fRmax  outer radius

//  History:
// 20.08.03 V.Grichine - created
// --------------------------------------------------------------------

#ifndef G4Orb_HH
#define G4Orb_HH

#include "G4CSGSolid.hh"

class G4Orb : public G4CSGSolid
{
  public:  // with description

    G4Orb(const G4String& pName, G4double pRmax);
       
    virtual ~G4Orb() ;
    
    // Accessors
       
    inline G4double GetRadius() const { return fRmax;};

    // Modifiers

    inline void SetRadius   (G4double newRmax) {fRmax=newRmax;};

    // Methods for solid

    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pmin, G4double& pmax) const;
         
    EInside Inside(const G4ThreeVector& p) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const;
    
    G4double DistanceToIn(const G4ThreeVector& p) const;
    
    G4double DistanceToOut(const G4ThreeVector& p,
                           const G4ThreeVector& v,
                           const G4bool calcNorm=G4bool(false),
                                 G4bool *validNorm=0,
                                 G4ThreeVector *n=0) const;
         
    G4double DistanceToOut(const G4ThreeVector& p) const;

    G4GeometryType GetEntityType() const;
 
    std::ostream& StreamInfo(std::ostream& os) const;

    // Visualisation functions
  
    void          DescribeYourselfTo(G4VGraphicsScene& scene) const;
    G4Polyhedron* CreatePolyhedron() const;
    G4NURBS*      CreateNURBS() const;

  protected:
  
    // Used by distanceToOut
  
    enum ESide {kNull,kRMax};
  
    // used by normal
  
    enum ENorm {kNRMax};

  private:

    static const G4double fEpsilon;
    G4double fRmax;
    G4double fRmaxTolerance;
};

#endif
