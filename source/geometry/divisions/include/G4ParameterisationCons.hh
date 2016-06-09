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
// $Id: G4ParameterisationCons.hh,v 1.5 2004/05/13 14:57:12 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// classes G4ParameterisationConsRho,
//         G4ParameterisationConsPhi,
//         G4ParameterisationConsZ
//
// Class description:
//
// These classes represent the parameterised positioning equivalent to 
// dividing a G4Cons along one of each axis Rho, Phi, Z.

// History:
// -------
// 09.05.01 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// --------------------------------------------------------------------
#ifndef G4ParameterisationCons_H
#define G4ParameterisationCons_H 1

#include "G4VDivisionParameterisation.hh"

class G4VSolid;
class G4VPhysicalVolume;

// Dummy declarations to get rid of warnings ...
class G4Trd;
class G4Trap;
class G4Box;
class G4Sphere;
class G4Orb;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;

class G4VParameterisationCons : public G4VDivisionParameterisation
{ 
  public:  // with description
  
    G4VParameterisationCons( EAxis axis, G4int nCopies,
                            G4double offset, G4double step,
                            G4VSolid* msolid, DivisionType divType );
  
    virtual ~G4VParameterisationCons();
};

class G4ParameterisationConsRho : public G4VParameterisationCons
{ 
  public:  // with description

    G4ParameterisationConsRho( EAxis axis, G4int nCopies,
                               G4double offset, G4double step,
                               G4VSolid* motherSolid, DivisionType divType );
    virtual ~G4ParameterisationConsRho();

    virtual G4double GetMaxParameter() const;

    virtual void ComputeTransformation( const G4int copyNo,
                                        G4VPhysicalVolume* physVol ) const;
    void ComputeDimensions( G4Cons& tubs, const G4int copyNo,
                            const G4VPhysicalVolume* physVol) const;

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
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const {}
};

class G4ParameterisationConsPhi : public G4VParameterisationCons
{ 
  public:  // with description

    G4ParameterisationConsPhi( EAxis axis, G4int nCopies,
                               G4double offset, G4double step,
                               G4VSolid* motherSolid, DivisionType divType );
    virtual ~G4ParameterisationConsPhi();

    virtual G4double GetMaxParameter() const;

    virtual void ComputeTransformation( const G4int copyNo,
                                        G4VPhysicalVolume* physVol ) const;
    void ComputeDimensions( G4Cons& tubs, const G4int copyNo,
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
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const {}
};

class G4ParameterisationConsZ : public G4VParameterisationCons
{ 
  public:  // with description

    G4ParameterisationConsZ( EAxis axis, G4int nCopies,
                             G4double offset, G4double step,
                             G4VSolid* motherSolid, DivisionType divType );
    virtual ~G4ParameterisationConsZ();

    virtual G4double GetMaxParameter() const;

    virtual void ComputeTransformation( const G4int copyNo,
                                        G4VPhysicalVolume* physVol ) const;
    void ComputeDimensions( G4Cons& tubs, const G4int copyNo,
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
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const {}
};

#endif
