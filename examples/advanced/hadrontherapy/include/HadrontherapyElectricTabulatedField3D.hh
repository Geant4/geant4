#include "globals.hh"
#include "G4ElectricField.hh"
#include "G4ElectroMagneticField.hh"
#include "G4ios.hh"

#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

class HadrontherapyElectricTabulatedField3D

 : public G4ElectricField

{
  
  // Storage space for the table
  vector< vector< vector< G4double > > > xEField;
  vector< vector< vector< G4double > > > yEField;
  vector< vector< vector< G4double > > > zEField;
  // The dimensions of the table
  G4int Enx,Eny,Enz; 
  // The physical limits of the defined region
  G4double Eminx, Emaxx, Eminy, Emaxy, Eminz, Emaxz;
  // The physical extent of the defined region
  G4double dx1, dy1, dz1;
  G4double feXoffset;
  G4double feYoffset;
  G4double feZoffset;
  G4bool einvertX, einvertY, einvertZ;

public:
  HadrontherapyElectricTabulatedField3D(const char* filename, G4double exOffset, G4double eyOffset, G4double ezOffset );
  void  GetFieldValue( const  G4double Epoint[4],
		       G4double *Efield) const;
};

