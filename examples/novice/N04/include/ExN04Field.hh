
#ifndef ExN04Field_H
#define ExN04Field_H 1

#include "globals.hh"
#include "G4MagneticField.hh"

class ExN04Field : public G4MagneticField
{
  public:
    ExN04Field();
    ~ExN04Field();

    void GetFieldValue( const  double Point[3],
                               double *Bfield ) const;

  private:
    G4double Bz;
    G4double rmax_sq;
    G4double zmax;
};

#endif

