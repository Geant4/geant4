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
///////////////////////////////////////////////////////////////////////////////
// File: CCalMagneticField.cc
// Description: User Field class implementation.
///////////////////////////////////////////////////////////////////////////////
#include "CCalMagneticField.hh"
#include "CCalutils.hh"
#include "G4FieldManager.hh"
#include "g4std/fstream"

//#define ddebug
//#define debug

//Constructor and destructor:

CCalMagneticField::CCalMagneticField(const G4String &filename) :
  fval(0), pos(0), slope(0), intercept(0) {

  //Let's open the file
  G4cout << " ==> Opening file " << filename << " to read magnetic field..."
       << G4endl;
  G4String pathName = getenv("CCAL_GLOBALPATH");
  G4std::ifstream is;
  bool ok = openGeomFile(is, pathName, filename);

  if (ok) {
    findDO(is, G4String("FLDM"));
    is >> fval >> npts >> xoff;
#ifdef debug
    G4cout << "Field value " << fval << " # points " << npts << " offset in x "
	 << xoff*mm << G4endl;
#endif

    if (npts > 0) {
      pos       = new G4double[npts];
      slope     = new G4double[npts];
      intercept = new G4double[npts];

      for (G4int i = 0; i < npts; i++) {
	is >> pos[i] >> slope[i] >> intercept[i];
#ifdef debug
	G4cout << tab << "Position " << i << " " << pos[i] << " Slope "
	     << slope[i] << " Intercept " << intercept[i] << G4endl;
#endif
      }
    }

    ///////////////////////////////////////////////////////////////
    // Close the file
    G4cout << " ==> Closing file " << filename << G4endl;
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
  //  G4cout << "X: " << xnew << " Slope " << m << " " << c << " " << scor << G4endl;
  B[2] = scor*fval*kilogauss;
#ifdef ddebug
  G4cout << "Field at x: " << x[0]/mm << "mm = " << B[2]/tesla << "T" << G4endl;
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

