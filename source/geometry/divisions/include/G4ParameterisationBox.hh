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
// $Id: G4ParameterisationBox.hh,v 1.1 2003-06-16 15:11:40 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// classes G4ParameterisationBoxX,
//         G4ParameterisationBoxY,
//         G4ParameterisationBoxZ
//
// Class description:
//
// These classes represent the parameterised positioning equivalent to 
// dividing a G4Box along one of each axis X, Y, Z.

// History:
// 09.05.01 - P.Arce Initial version
// ********************************************************************

#ifndef G4ParameterisationBox_H
#define G4ParameterisationBox_H 1

#include "G4VDivisionParameterisation.hh"

class G4VSolid;
class G4VPhysicalVolume;
class HepTransform3D;

class G4ParameterisationBoxX : public G4VDivisionParameterisation
{ 
  public:  // with description

    G4ParameterisationBoxX( EAxis axis, G4int nCopies,
                            G4double offset, G4double step,
                            G4VSolid* msolid, DivisionType divType );

    virtual ~G4ParameterisationBoxX();

    virtual void ComputeTransformation( const G4int copyNo,
                                        G4VPhysicalVolume* physVol ) const;
};

class G4ParameterisationBoxY : public G4VDivisionParameterisation
{ 
  public:  // with description

    G4ParameterisationBoxY( EAxis axis, G4int nCopies,
                            G4double offset, G4double step,
                            G4VSolid* msolid, DivisionType divType );
    virtual ~G4ParameterisationBoxY();

    virtual void ComputeTransformation( const G4int copyNo,
                                        G4VPhysicalVolume* physVol ) const;
};

class G4ParameterisationBoxZ : public G4VDivisionParameterisation
{ 
  public:  // with description

  G4ParameterisationBoxZ( EAxis axis, G4int nCopies,
                          G4double offset, G4double step,
                          G4VSolid* msolid, DivisionType divType );
  virtual ~G4ParameterisationBoxZ();

  virtual void ComputeTransformation( const G4int copyNo,
                                      G4VPhysicalVolume* physVol ) const;
};

#endif
