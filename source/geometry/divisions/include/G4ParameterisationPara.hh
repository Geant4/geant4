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
// $Id: G4ParameterisationPara.hh,v 1.6 2004/05/13 14:57:12 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-02 $
// 
// classes G4ParameterisationParaX,
//         G4ParameterisationParaY,
//         G4ParameterisationParaZ
//
// Class description:
//
// These classes represent the parameterised positioning equivalent to 
// dividing a G4Para along one of each axis X, Y, Z.

// History:
// -------
// 09.05.01 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// --------------------------------------------------------------------
#ifndef G4ParameterisationPara_H
#define G4ParameterisationPara_H 1

#include "G4VDivisionParameterisation.hh"

class G4VSolid;
class G4VPhysicalVolume;
// Dummy declarations to get rid of warnings ...
class G4Cons;
class G4Cons;
class G4Sphere;
class G4Orb;
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
  
    virtual ~G4VParameterisationPara();
};

class G4ParameterisationParaX : public G4VParameterisationPara
{ 
  public:  // with description

    G4ParameterisationParaX( EAxis axis, G4int nCopies,
                             G4double offset, G4double step,
                             G4VSolid* msolid, DivisionType divType );
    virtual ~G4ParameterisationParaX();

    virtual G4double GetMaxParameter() const;

    virtual void ComputeTransformation( const G4int copyNo,
                                        G4VPhysicalVolume* physVol ) const;
    void ComputeDimensions(G4Para& para, const G4int copyNo,
                           const G4VPhysicalVolume* pv) const;

  
  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Box&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&,const G4int,
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


class G4ParameterisationParaY : public G4VParameterisationPara
{ 
  public:  // with description

    G4ParameterisationParaY( EAxis axis, G4int nCopies,
                             G4double offset, G4double step,
                             G4VSolid* msolid, DivisionType divType );
    virtual ~G4ParameterisationParaY();
  
    virtual G4double GetMaxParameter() const;

    virtual void ComputeTransformation( const G4int copyNo,
                                        G4VPhysicalVolume* physVol ) const;
    void ComputeDimensions(G4Para& para, const G4int copyNo,
                           const G4VPhysicalVolume* pv) const;

  
  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Box&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&,const G4int,
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


class G4ParameterisationParaZ : public G4VParameterisationPara
{ 
  public:  // with description

    G4ParameterisationParaZ( EAxis axis, G4int nCopies,
                             G4double offset, G4double step,
                             G4VSolid* msolid, DivisionType divType );
    virtual ~G4ParameterisationParaZ();

    virtual G4double GetMaxParameter() const;

    virtual void ComputeTransformation( const G4int copyNo,
                                        G4VPhysicalVolume* physVol ) const;
    void ComputeDimensions(G4Para& para, const G4int copyNo,
                           const G4VPhysicalVolume* pv) const;

  
  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Box&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&,const G4int,
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
