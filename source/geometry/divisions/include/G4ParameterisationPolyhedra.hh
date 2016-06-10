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
//
// $Id: G4ParameterisationPolyhedra.hh 73433 2013-08-27 11:05:39Z gcosmo $
// 
// classes G4ParameterisationPolyhedraRho,
//         G4ParameterisationPolyhedraPhi,
//         G4ParameterisationPolyhedraZ
//
// Class description:
//
// These classes represent the parameterised positioning equivalent to 
// dividing a G4Polyhedra along one of each axis Rho, Phi, Z.

// History:
// -------
// 09.05.01 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
//---------------------------------------------------------------------
#ifndef G4ParameterisationPolyhedra_H
#define G4ParameterisationPolyhedra_H 1

#include "G4VDivisionParameterisation.hh"
#include "G4Polyhedra.hh"

class G4VPhysicalVolume;

// Dummy declarations to get rid of warnings ...
//
class G4Trd;
class G4Trap;
class G4Cons;
class G4Sphere;
class G4Orb;
class G4Ellipsoid;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;

class G4VParameterisationPolyhedra : public G4VDivisionParameterisation
{ 
  public:  // with description
  
    G4VParameterisationPolyhedra( EAxis axis, G4int nCopies,
                            G4double offset, G4double step,
                            G4VSolid* msolid, DivisionType divType );
  
    virtual ~G4VParameterisationPolyhedra();

  private:

    G4double ConvertRadiusFactor(const G4Polyhedra& phedra) const;
      // Converts radius of the sides to radius of the corners:
      // r_corners = r_sides/factor
};

//---------------------------------------------------------------------
// Class G4ParameterisationPolyhedraRho
//---------------------------------------------------------------------

class G4ParameterisationPolyhedraRho : public G4VParameterisationPolyhedra
{ 
  public:  // with description

    G4ParameterisationPolyhedraRho( EAxis axis, G4int nCopies,
                                    G4double offset, G4double step,
                                    G4VSolid* motherSolid,
                                    DivisionType divType );
   ~G4ParameterisationPolyhedraRho();

    void CheckParametersValidity();

    G4double GetMaxParameter() const;

    void ComputeTransformation( const G4int copyNo,
                                      G4VPhysicalVolume* physVol ) const;
    void ComputeDimensions( G4Polyhedra& phedra, const G4int copyNo,
                            const G4VPhysicalVolume* physVol ) const;

  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Box&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const {}
};

//---------------------------------------------------------------------
// Class G4ParameterisationPolyhedraPhi
//---------------------------------------------------------------------

class G4ParameterisationPolyhedraPhi : public G4VParameterisationPolyhedra
{ 
  public:  // with description

    G4ParameterisationPolyhedraPhi( EAxis axis, G4int nCopies,
                                   G4double offset, G4double step,
                                   G4VSolid* motherSolid,
                                   DivisionType divType );
   ~G4ParameterisationPolyhedraPhi();

    void CheckParametersValidity();

    G4double GetMaxParameter() const;

    void ComputeTransformation( const G4int copyNo,
                                      G4VPhysicalVolume* physVol ) const;
    void ComputeDimensions( G4Polyhedra& phedra, const G4int copyNo,
                            const G4VPhysicalVolume* physVol ) const;

  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Box&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const {}
};

//---------------------------------------------------------------------
// Class G4ParameterisationPolyhedraZ
//---------------------------------------------------------------------

class G4ParameterisationPolyhedraZ : public G4VParameterisationPolyhedra
{ 
  public:  // with description

    G4ParameterisationPolyhedraZ( EAxis axis, G4int nCopies,
                                  G4double offset, G4double step,
                                  G4VSolid* motherSolid,
                                  DivisionType divType );
   ~G4ParameterisationPolyhedraZ();

    void CheckParametersValidity();

    G4double GetMaxParameter() const;

    void ComputeTransformation( const G4int copyNo,
                                      G4VPhysicalVolume* physVol ) const;
    void ComputeDimensions( G4Polyhedra& phedra, const G4int copyNo,
                            const G4VPhysicalVolume* physVol ) const;

  private:

    G4double GetR(G4double z, G4double z1, G4double r1,
                  G4double z2, G4double r2) const;
    G4double GetRmin(G4double z, G4int nsegment) const;
    G4double GetRmax(G4double z, G4int nsegment) const;

    // Dummy declarations to get rid of warnings ...
    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Box&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const {}
  private:

    G4int fNSegment;
    G4PolyhedraHistorical* fOrigParamMother;
};

#endif
