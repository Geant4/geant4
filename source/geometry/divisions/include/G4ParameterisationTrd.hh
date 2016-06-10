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
// $Id: G4ParameterisationTrd.hh 73433 2013-08-27 11:05:39Z gcosmo $
//
// classes G4ParameterisationTrdX
//         G4ParameterisationTrdY
//         G4ParameterisationTrdZ
//
// Class description:
//
// This class represents the parameterised positioning equivalent to 
// dividing a trapezoid along one of each axis X, Y, Z.

// History:
// -------
// 09.05.01 - P.Arce, Initial version
// 08.04.04 - I.Hrivnacova, Implemented reflection
// --------------------------------------------------------------------
#ifndef G4ParameterisationTrd_H
#define G4ParameterisationTrd_H 1

#include <vector>

#include "G4VDivisionParameterisation.hh"
#include "G4VSolid.hh" 

class G4VPhysicalVolume;

// Dummy declarations to get rid of warnings ...
//
class G4Cons;
class G4Box;
class G4Sphere;
class G4Orb;
class G4Ellipsoid;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;

class G4VParameterisationTrd : public G4VDivisionParameterisation
{ 
  public:  // with description
  
    G4VParameterisationTrd( EAxis axis, G4int nCopies,
                            G4double offset, G4double step,
                            G4VSolid* msolid, DivisionType divType );
  
    virtual ~G4VParameterisationTrd();

  protected:

    G4bool bDivInTrap;
};

class G4ParameterisationTrdX : public G4VParameterisationTrd
{ 
  public:  // with description

    G4ParameterisationTrdX( EAxis axis, G4int nCopies,
                            G4double width, G4double offset,
                            G4VSolid* motherSolid, DivisionType divType );
   ~G4ParameterisationTrdX();

    void CheckParametersValidity();

    G4double GetMaxParameter() const;

    void ComputeTransformation(const G4int copyNo,
                                     G4VPhysicalVolume* physVol) const;

    void ComputeDimensions(G4Trd& trd, const G4int copyNo,
                           const G4VPhysicalVolume* pv) const;

    void ComputeDimensions(G4Trap& trd, const G4int copyNo,
                           const G4VPhysicalVolume* pv) const;
  
    G4VSolid*   ComputeSolid(const G4int, G4VPhysicalVolume *);


  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Cons&,const G4int,
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
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const {}

    void ComputeTrapParams();
};


class G4ParameterisationTrdY : public G4VParameterisationTrd
{ 
  public:  // with description

    G4ParameterisationTrdY( EAxis axis, G4int nCopies,
                            G4double width, G4double offset,
                            G4VSolid* motherSolid, DivisionType divType );
   ~G4ParameterisationTrdY();

    void CheckParametersValidity();

    G4double GetMaxParameter() const;

    void ComputeTransformation(const G4int copyNo,
                                     G4VPhysicalVolume *physVol) const;
 
    void ComputeDimensions(G4Trd& trd, const G4int copyNo,
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
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,
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


class G4ParameterisationTrdZ : public G4VParameterisationTrd
{ 
  public:  // with description

    G4ParameterisationTrdZ( EAxis axis, G4int nCopies,
                            G4double width, G4double offset,
                            G4VSolid* motherSolid, DivisionType divType );
   ~G4ParameterisationTrdZ();

    G4double GetMaxParameter() const;

    void ComputeTransformation(const G4int copyNo,
                                     G4VPhysicalVolume* physVol) const;
    void ComputeDimensions(G4Trd& trd, const G4int copyNo,
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
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,
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
