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
///////////////////////////////////////////////////////////////////////////////
// File: CCalMagneticField.cc
// Description: User Field class implementation.
///////////////////////////////////////////////////////////////////////////////
#include <fstream>

#include "CCalMagneticField.hh"
#include "CCalutils.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldManager.hh"

//#define ddebug
//#define debug

//Constructor and destructor:

CCalMagneticField::CCalMagneticField(const G4String &filename) :
  fval(0), pos(0), slope(0), intercept(0) 
{
#ifdef debug
  fVerbosity = 1;
#else
  fVerbosity = 0;
#endif

  //Let's open the file
  G4cout << " ==> Opening file " << filename << " to read magnetic field..."
         << G4endl;
  G4String pathName = std::getenv("CCAL_GLOBALPATH");
  std::ifstream is;
  G4bool ok = openGeomFile(is, pathName, filename);
  
  if (ok) {
    findDO(is, G4String("FLDM"));
    is >> fval >> npts >> xoff;

    if (fVerbosity)
      G4cout << "Field value " << fval << " # points " << npts << 
        " offset in x "
             << xoff*mm << G4endl;

    if (npts > 0) {
      pos       = new G4double[npts];
      slope     = new G4double[npts];
      intercept = new G4double[npts];

      for (G4int i = 0; i < npts; i++) {
        is >> pos[i] >> slope[i] >> intercept[i];
        if (fVerbosity)
          G4cout << tab << "Position " << i << " " << pos[i] << " Slope "
                 << slope[i] << " Intercept " << intercept[i] << G4endl;
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

void CCalMagneticField::MagneticField(const G4double x[3], G4double B[3]) const 
{
  G4int i=0;
  for (i=0; i<2; i++) {
    B[i]   = 0*kilogauss;
  }

  G4double m1=0;
  G4double c1=0;
  G4double xnew = x[0]/mm + xoff;
  if (npts > 0) {
    for (i=0; i<npts; i++) {
      if (xnew > pos[i]*mm) {
        m1 = slope[i];
        c1 = intercept[i];
      }
    }
  }
  G4double scor = c1 + m*xnew;
  if (scor < 0.) scor = 0.;
  if (scor > 1.) scor = 1.0;

  B[2] = scor*fval*kilogauss;
  if (fVerbosity)
    {

      G4cout << "Field at x: " << x[0]/mm << "mm (" << xnew << ") = " << 
        B[2]/tesla
             << "T (m = " << m1 << ", c = " << 
        c1 << ", scale = " << scor << ")"
             << G4endl;
    }
}


CLHEP::Hep3Vector CCalMagneticField::
MagneticField(const CLHEP::Hep3Vector point) const {

  G4double x[3],B[3];
  CLHEP::Hep3Vector v;
  
  x[0] = point.x();
  x[1] = point.y();
  x[2] = point.z();
  CCalMagneticField::MagneticField(x, B);
  v.setX(B[0]);   
  v.setY(B[1]);   
  v.setZ(B[2]);   
  return v;
}


void CCalMagneticField::GetFieldValue(const G4double x[3], G4double* B) const {
  CCalMagneticField::MagneticField(x, B);
}

