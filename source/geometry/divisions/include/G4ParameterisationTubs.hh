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
// $Id: G4ParameterisationTubs.hh,v 1.5 2004/05/13 14:57:13 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// classes G4ParameterisationTubsRho
//         G4ParameterisationTubsPhi
//         G4ParameterisationTubsZ
//
// Class description:
//
// This class represents the parameterised positioning equivalent to 
// dividing a G4Tubs along one of each axis Rho, Phi, Z.

// History:
// -------
// 09.05.01 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// --------------------------------------------------------------------
#ifndef G4ParameterisationTubsRho_H
#define G4ParameterisationTubsRho_H 1

#include "G4VDivisionParameterisation.hh"

class G4VPhysicalVolume;

// Dummy declarations to get rid of warnings ...
class G4Trd;
class G4Trap;
class G4Cons;
class G4Sphere;
class G4Orb;
class G4Torus;
class G4Para;
class G4Hype;
class G4Polycone;
class G4Polyhedra;

class G4VParameterisationTubs : public G4VDivisionParameterisation
{ 
  public:  // with description
  
    G4VParameterisationTubs( EAxis axis, G4int nCopies,
                            G4double offset, G4double step,
                            G4VSolid* msolid, DivisionType divType );
  
    virtual ~G4VParameterisationTubs();
};

class G4ParameterisationTubsRho : public G4VParameterisationTubs
{ 
  public:  // with description

    G4ParameterisationTubsRho( EAxis axis, G4int nCopies,
                               G4double offset, G4double step,
                               G4VSolid* motherSolid, DivisionType divType );
    virtual ~G4ParameterisationTubsRho();

    virtual G4double GetMaxParameter() const;

    virtual void ComputeTransformation(const G4int copyNo,
                                       G4VPhysicalVolume* physVol) const;
    void ComputeDimensions(G4Tubs& tubs, const G4int copyNo,
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
    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const {}
};


class G4ParameterisationTubsPhi : public G4VParameterisationTubs
{ 
  public:  // with description

    G4ParameterisationTubsPhi( EAxis axis, G4int nCopies,
                               G4double offset, G4double step,
                               G4VSolid* motherSolid, DivisionType divType );
    virtual ~G4ParameterisationTubsPhi();

    virtual G4double GetMaxParameter() const;

    virtual void ComputeTransformation(const G4int copyNo,
                                       G4VPhysicalVolume* physVol) const;
    void ComputeDimensions(G4Tubs& tubs, const G4int copyNo,
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
    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const {}
};


class G4ParameterisationTubsZ : public G4VParameterisationTubs
{ 
  public:  // with description

    G4ParameterisationTubsZ( EAxis axis, G4int nCopies,
                             G4double offset, G4double step,
                             G4VSolid* motherSolid, DivisionType divType );
    virtual ~G4ParameterisationTubsZ();

    virtual G4double GetMaxParameter() const;

    virtual void ComputeTransformation(const G4int copyNo,
                                       G4VPhysicalVolume* physVol) const;
    void ComputeDimensions(G4Tubs& tubs, const G4int copyNo,
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
    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const {}
};

#endif
