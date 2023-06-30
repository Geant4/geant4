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
// G4ParameterisationPara[X,Y,Z]
//
// Class description:
//
// These classes represent the parameterised positioning equivalent to 
// dividing a G4Para along one of each axis X, Y, Z.

// 09.05.01 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// --------------------------------------------------------------------
#ifndef G4PARAMETERISATIONPARA_HH
#define G4PARAMETERISATIONPARA_HH 1

#include "G4VDivisionParameterisation.hh"

class G4VSolid;
class G4VPhysicalVolume;

// Dummy declarations to get rid of warnings ...
//
class G4Cons;
class G4Cons;
class G4Sphere;
class G4Orb;
class G4Ellipsoid;
class G4Torus;
class G4Trd;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;

class G4VParameterisationPara : public G4VDivisionParameterisation
{ 
  public:  // with description
  
    G4VParameterisationPara( EAxis axis, G4int nCopies,
                            G4double offset, G4double step,
                            G4VSolid* msolid, DivisionType divType );
  
    ~G4VParameterisationPara() override;
};

class G4ParameterisationParaX : public G4VParameterisationPara
{ 
  public:  // with description

    G4ParameterisationParaX( EAxis axis, G4int nCopies,
                             G4double offset, G4double step,
                             G4VSolid* msolid, DivisionType divType );
   ~G4ParameterisationParaX() override;

    G4double GetMaxParameter() const override;

    void ComputeTransformation( const G4int copyNo,
                                      G4VPhysicalVolume* physVol ) const override;
    void ComputeDimensions(G4Para& para, const G4int copyNo,
                           const G4VPhysicalVolume* pv) const override;

  
  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Box&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const override {}
};


class G4ParameterisationParaY : public G4VParameterisationPara
{ 
  public:  // with description

    G4ParameterisationParaY( EAxis axis, G4int nCopies,
                             G4double offset, G4double step,
                             G4VSolid* msolid, DivisionType divType );
   ~G4ParameterisationParaY() override;
  
    G4double GetMaxParameter() const override;

    void ComputeTransformation( const G4int copyNo,
                                      G4VPhysicalVolume* physVol ) const override;
    void ComputeDimensions(G4Para& para, const G4int copyNo,
                           const G4VPhysicalVolume* pv) const override;

  
  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Box&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const override {}
};


class G4ParameterisationParaZ : public G4VParameterisationPara
{ 
  public:  // with description

    G4ParameterisationParaZ( EAxis axis, G4int nCopies,
                             G4double offset, G4double step,
                             G4VSolid* msolid, DivisionType divType );
   ~G4ParameterisationParaZ() override;

    G4double GetMaxParameter() const override;

    void ComputeTransformation( const G4int copyNo,
                                      G4VPhysicalVolume* physVol ) const override;
    void ComputeDimensions(G4Para& para, const G4int copyNo,
                           const G4VPhysicalVolume* pv) const override;

  
  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Box&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const override {}
};

#endif
