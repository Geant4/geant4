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
//
/// file: EllipticalChromosome.hh
/// brief: declaration file for a rod-shaped chromosome model.

#ifndef MOLECULAR_ELLIPTICAL_CHROMOSOME_HH
#define MOLECULAR_ELLIPTICAL_CHROMOSOME_HH

#include "VirtualChromosome.hh"

#include "globals.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EllipticalChromosome : public VirtualChromosome
{
 public:
  EllipticalChromosome(const G4String&, G4ThreeVector, const G4double&,
                       const G4double&, const G4double&);

  EllipticalChromosome(const G4String&, G4ThreeVector, const G4double&,
                       const G4double&, const G4double&, G4RotationMatrix);

  ~EllipticalChromosome() override;

  G4bool PointInChromosome(G4ThreeVector const& position) override;

  G4ThreeVector RandomPointInChromosome() override;

  G4String GetShape() override { return fShape; };
  static const G4String fShape;

  void Print() override
  {
    G4cout << "type: " << fShape << G4endl;
    G4cout << "semi-x: " << fSemiX << G4endl;
    G4cout << "semi-y: " << fSemiY << G4endl;
    G4cout << "semi-z: " << fSemiZ << G4endl;
    G4cout << "center: " << fCenter << G4endl;
    G4cout << "rotation: " << fRotation.getPhi() << " " << fRotation.getTheta()
           << " " << fRotation.getPhi() << G4endl;
  }

  void SetRotation(const G4RotationMatrix& rot)
  {
    fRotation        = rot;
    fInverseRotation = fRotation.inverse();
  };

  G4RotationMatrix GetRotation() const { return fRotation; };

  G4ThreeVector GetPosition() const { return fCenter; };

  G4double GetSemiX() const { return fSemiX; };

  G4double GetSemiY() const { return fSemiY; };

  G4double GetSemiZ() const { return fSemiZ; };

 private:
  G4ThreeVector fCenter;
  G4double fSemiX, fSemiY, fSemiZ;
  G4RotationMatrix fRotation, fInverseRotation;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_ELLIPTICAL_CHROMOSOME_HH
