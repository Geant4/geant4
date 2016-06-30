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
// $Id: G4UniformMagField.hh 96751 2016-05-04 09:39:38Z gcosmo $
//
// 
// class G4UniformMagField
//
// Class description:
//
// Class for creation of Uniform Magnetic Field.

// History:
// - 30.01.97 V.Grichine, Created.
// - 01.08.97 J.Apostolakis, cleanup, new 3-vector constructor, 
//            and removal of helix-stepper (to separate file).
// - 05.11.97 G.Cosmo, added copy constructor and assignment operator.

#ifndef G4UNIFORMMAGFIELD_HH
#define G4UNIFORMMAGFIELD_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4MagneticField.hh"

class G4UniformMagField : public G4MagneticField
{
  public:  // with description
  
    G4UniformMagField(const G4ThreeVector& FieldVector );
      // A field with value equal to FieldVector.

    G4UniformMagField(G4double vField,
                      G4double vTheta,
                      G4double vPhi     ) ;

    virtual ~G4UniformMagField() ;

    G4UniformMagField(const G4UniformMagField &p);
    G4UniformMagField& operator = (const G4UniformMagField &p);
      // Copy constructor and assignment operator.

    virtual void GetFieldValue(const G4double yTrack[4],
                                     G4double *MagField) const ;

    void SetFieldValue(const G4ThreeVector& newFieldValue);

    G4ThreeVector GetConstantFieldValue() const;
      // Return the field value
    
    virtual G4Field* Clone() const;

  private:

    G4double fFieldComponents[3] ;
};

#endif
