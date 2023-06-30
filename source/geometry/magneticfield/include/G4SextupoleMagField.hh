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
// Class description:
//
// Sextupole magnetic field
//   fGradient - is the gradient value for a sextupole magnet
//   the magnetic field components are:
//   Bx = B[0] = fGradient*X*Y,
//   By = B[1] = fGradient*(X*X-Y*Y)/2,
//   Bz = B[2] = 0

// Author: H. Burkhardt - 23/10/2019
// -------------------------------------------------------------------
#ifndef G4SEXTUPOLEMAGFIELD_HH
#define G4SEXTUPOLEMAGFIELD_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class G4SextupoleMagField : public G4MagneticField
{
  public:

    G4SextupoleMagField(G4double pGradient);
    G4SextupoleMagField(      G4double pGradient,
                        const G4ThreeVector& pOrigin,
                        G4RotationMatrix* pMatrix);
   ~G4SextupoleMagField() override;

    void GetFieldValue(const G4double yTrack[], G4double B[] ) const override;
    G4Field* Clone() const override;

  private:

    G4double          fGradient = 0.0;
    G4ThreeVector     fOrigin   = G4ThreeVector(0.0, 0.0, 0.0);
    G4RotationMatrix* fpMatrix  = nullptr;
};

#endif
