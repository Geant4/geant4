
#ifndef field_h
#define field_h 1

#include "G4UniformMagField.hh"
//#include "G4Field.hh"
//#include "G4ThreeVector.hh"
//#include "G4Types.hh"

class FieldMessenger;

class Field : public G4UniformMagField
{
  public:
    Field();
    Field(G4double);
    virtual ~Field();

    //virtual void GetFieldValue(const double*, double*);
    void SetMagFieldValue(G4double field);

  private:
    FieldMessenger* fFieldMessenger;
    //G4double fFieldVal;

};

#endif 
