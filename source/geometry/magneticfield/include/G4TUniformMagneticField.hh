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

#ifndef G4UniformMagneticField_HH
#define G4UniformMagneticField_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4MagneticField.hh"

class G4TUniformMagneticField : public G4MagneticField
{
  public:  // with description

    G4TUniformMagneticField(const G4ThreeVector& FieldVector ) 
    // A field with value equal to FieldVector.
    {
      fFieldComponents[0] = FieldVector.x();
      fFieldComponents[1] = FieldVector.y();
      fFieldComponents[2] = FieldVector.z();
    }
   

    G4TUniformMagneticField(G4double vField,
                            G4double vTheta,
                            G4double vPhi     ) 
    {
      if ( (vField<0) || (vTheta<0) || (vTheta>pi) || (vPhi<0) || (vPhi>twopi) )
      {
         G4Exception("G4TUniformMagneticField::G4TUniformMagneticField()",
                     "GeomField0002", FatalException, "Invalid parameters.") ;
      }
      fFieldComponents[0] = vField*std::sin(vTheta)*std::cos(vPhi) ;
      fFieldComponents[1] = vField*std::sin(vTheta)*std::sin(vPhi) ;
      fFieldComponents[2] = vField*std::cos(vTheta) ;
    }

    virtual ~G4TUniformMagneticField() {;}

    G4TUniformMagneticField(const G4TUniformMagneticField &p)
       : G4MagneticField(p)
    {
       for (G4int i=0; i<3; ++i)
          fFieldComponents[i] = p.fFieldComponents[i];
    }
   
    G4TUniformMagneticField& operator = (const G4TUniformMagneticField &p)
            // Copy constructor and assignment operator.
    {
      if (&p == this) return *this;
      for (G4int i=0; i<3; ++i)
         fFieldComponents[i] = p.fFieldComponents[i];
      return *this;
    }

    inline void GetFieldValue(const G4double yTrack[4],
                              G4double *B) const 
    {
       B[0]= fFieldComponents[0] ;
       B[1]= fFieldComponents[1] ;
       B[2]= fFieldComponents[2] ;
    }
   
    void SetFieldValue(const G4ThreeVector& newFieldVector)
    {
       fFieldComponents[0] = newFieldVector.x();
       fFieldComponents[1] = newFieldVector.y();
       fFieldComponents[2] = newFieldVector.z();
    }

    G4ThreeVector GetConstantFieldValue() const
    {
       G4ThreeVector B(fFieldComponents[0],
                       fFieldComponents[1],
                       fFieldComponents[2]);
       return B;
    }
    // Return the field value
   
    virtual G4TUniformMagneticField* Clone() const
    { 
       return new G4TUniformMagneticField( G4ThreeVector(this->fFieldComponents[0],
                                                         this->fFieldComponents[1],
                                                         this->fFieldComponents[2]) );
    }

    private:
        G4double fFieldComponents[3] ;
};

#endif
