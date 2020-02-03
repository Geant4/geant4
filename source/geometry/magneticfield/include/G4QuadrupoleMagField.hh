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
// G4QuadrupoleMagField
//
// Class description:
//
// Class for creation of quadrupole magnetic field
//   fGradient - is the gradient value for quadrupole magnetic lense.
// Then the magnetic field components are:
//   Bx = B[0] = fGradient*X ,
//   By = B[1] = fGradient*Y ,
//   Bz = B[2] = 0 .
// Here X,Y,Z are the coordinates of a space point of interest.

// 03.02.1997, V.Grichine - Created
// 11.05.2012, B.Riese  - Allow displaced origin and rotation
// -------------------------------------------------------------------
#ifndef G4QUADRUPOLEMAGFIELD_HH
#define G4QUADRUPOLEMAGFIELD_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class G4QuadrupoleMagField : public G4MagneticField
{
  public: // with description

    G4QuadrupoleMagField(G4double          pGradient);
    G4QuadrupoleMagField(G4double          pGradient, 
                         G4ThreeVector     pOrigin, 
                         G4RotationMatrix* pMatrix);
   ~G4QuadrupoleMagField();

    void GetFieldValue(const G4double yTrack[],
                             G4double B[]     ) const;
    G4Field* Clone() const;

  private:

    G4double          fGradient = 0.0;
    G4ThreeVector     fOrigin = G4ThreeVector(0.0, 0.0, 0.0);
    G4RotationMatrix* fpMatrix = nullptr;
};

#endif

