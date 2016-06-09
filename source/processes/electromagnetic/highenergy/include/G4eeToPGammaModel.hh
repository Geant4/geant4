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
// $Id: G4eeToPGammaModel.hh,v 1.1 2008/07/10 18:07:26 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eeToPGammaModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 10.07.2008
//
// Modifications:
//

//
// Class Description:
//

// -------------------------------------------------------------------
//

#ifndef G4eeToPGammaModel_h
#define G4eeToPGammaModel_h 1

#include "G4Vee2hadrons.hh"
#include "globals.hh"
#include "G4eeCrossSections.hh"
#include "G4ParticleDefinition.hh"

class G4DynamicParticle;
class G4PhysicsVector;

class G4eeToPGammaModel : public G4Vee2hadrons
{

public:

  G4eeToPGammaModel(G4eeCrossSections*, const G4String&);

  virtual ~G4eeToPGammaModel();

  G4double ThresholdEnergy() const;

  G4double PeakEnergy() const;

  G4double ComputeCrossSection(G4double) const;

  G4PhysicsVector* PhysicsVector(G4double, G4double) const;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
              G4double, const G4ThreeVector&);

private:

  // hide assignment operator
  G4eeToPGammaModel & operator=(const  G4eeToPGammaModel &right);
  G4eeToPGammaModel(const  G4eeToPGammaModel&);

  G4eeCrossSections* cross;

  G4ParticleDefinition* particle;
  G4ParticleDefinition* pi0;

  G4double massP;
  G4double massR;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToPGammaModel::ThresholdEnergy() const
{
  return LowEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToPGammaModel::PeakEnergy() const
{
  return massR;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToPGammaModel::ComputeCrossSection(G4double e) const
{
  G4double ee = std::min(HighEnergy(),e);
  G4double xs;
  if(particle == pi0) xs = cross->CrossSectionPi0G(ee);
  else                xs = cross->CrossSectionEtaG(ee);
  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
