// $Id: A01MagneticField.hh,v 1.1 2002-11-13 07:18:38 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef A01MagneticField_H
#define A01MagneticField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"
class A01MagneticFieldMessenger;

class A01MagneticField : public G4MagneticField
{
  public:
    A01MagneticField();
    ~A01MagneticField();

    virtual void GetFieldValue( const  double Point[3],
                               double *Bfield ) const;

  private:
    A01MagneticFieldMessenger* messenger;
    G4double By;
    G4double rmax_sq;
    G4double ymax;

  public:
    inline void SetField(G4double val) { By = val; }
    inline G4double GetField() const { return By; }
};

#endif

