///////////////////////////////////////////////////////////////////////////////
// File: CCalMagneticField.cc
// Description: User Field class implementation.
///////////////////////////////////////////////////////////////////////////////
#include "CCalMagneticField.hh"
#include "CCalutils.hh"
#include "G4FieldManager.hh"
#include <fstream.h>

//#define ddebug
//#define debug

//Constructor and destructor:

CCalMagneticField::CCalMagneticField(const G4String &filename) :
  fval(0), pos(0), slope(0), intercept(0) {

  //Let's open the file
  cout << " ==> Opening file " << filename << " to read magnetic field..."
       << endl;
  G4String pathName = getenv("CCAL_GLOBALPATH");
  ifstream is;
  bool ok = openGeomFile(is, pathName, filename);

  if (ok) {
    findDO(is, G4String("FLDM"));
    is >> fval >> npts >> xoff;
#ifdef debug
    cout << "Field value " << fval << " # points " << npts << " offset in x "
	 << xoff*mm << endl;
#endif

    if (npts > 0) {
      pos       = new G4double[npts];
      slope     = new G4double[npts];
      intercept = new G4double[npts];

      for (G4int i = 0; i < npts; i++) {
	is >> pos[i] >> slope[i] >> intercept[i];
#ifdef debug
	cout << tab << "Position " << i << " " << pos[i] << " Slope "
	     << slope[i] << " Intercept " << intercept[i] << endl;
#endif
      }
    }

    ///////////////////////////////////////////////////////////////
    // Close the file
    cout << " ==> Closing file " << filename << endl;
    is.close();
  }
}


CCalMagneticField::~CCalMagneticField() {
  if (pos)
    delete[] pos;
  if (slope)
    delete[] slope;
  if (intercept)
    delete[] intercept;
}


// Member functions

void CCalMagneticField::MagneticField(const double x[3], double B[3]) const {

  G4int i=0;
  for (i=0; i<2; i++) {
    B[i]   = 0*kilogauss;
  }

  G4double m=0, c=1;
  G4double xnew = x[0]/mm + xoff;
  if (npts > 0) {
    m = slope[npts-1];
    c = intercept[i-1];
    for (i=npts-2; i>=0; i--) {
      if (xnew < pos[i]*mm) {
	m = slope[i];
        c = intercept[i];
      }
    }
  }
  G4double scor = c + m*xnew;
  if (scor < 0.) scor = 0.;
  //  cout << "X: " << xnew << " Slope " << m << " " << c << " " << scor << endl;
  B[2] = scor*fval*kilogauss;
#ifdef ddebug
  cout << "Field at x: " << x[0]/mm << "mm = " << B[2]/tesla << "T" << endl;
#endif
}


Hep3Vector CCalMagneticField::MagneticField(const Hep3Vector point) const {

  G4double x[3],B[3];
  Hep3Vector v;
  
  x[0] = point.x();
  x[1] = point.y();
  x[2] = point.z();
  CCalMagneticField::MagneticField(x, B);
  v.setX(B[0]);   
  v.setY(B[1]);   
  v.setZ(B[2]);   
  return v;
}


void CCalMagneticField::GetFieldValue(const double x[3], double* B) const {
  CCalMagneticField::MagneticField(x, B);
}

