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
// $Id: G4eeToTwoPiModel.hh,v 1.1 2004/11/19 18:44:04 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eeToTwoPiModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 25.10.2003
//
// Modifications:
//

//
// Class Description:
//

// -------------------------------------------------------------------
//

#ifndef G4eeToTwoPiModel_h
#define G4eeToTwoPiModel_h 1

#include "G4Vee2hadrons.hh"
#include "globals.hh"
#include "G4eeCrossSections.hh"

class G4DynamicParticle;
class G4PhysicsVector;

class G4eeToTwoPiModel : public G4Vee2hadrons
{

public:

  G4eeToTwoPiModel(G4eeCrossSections*);

  virtual ~G4eeToTwoPiModel();

  G4double ThresholdEnergy() const;

  G4double PeakEnergy() const;

  G4double ComputeCrossSection(G4double) const;

  G4PhysicsVector* PhysicsVector(G4double, G4double) const;

  virtual std::vector<G4DynamicParticle*>* SampleSecondaries(
              G4double, const G4ThreeVector&) const;

private:

  void Initialise();

  // hide assignment operator
  G4eeToTwoPiModel & operator=(const  G4eeToTwoPiModel &right);
  G4eeToTwoPiModel(const  G4eeToTwoPiModel&);

  G4eeCrossSections* cross;

  G4double massPi;
  G4double massRho;
  G4double highEnergy;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToTwoPiModel::ThresholdEnergy() const
{
  return 2.0*massPi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToTwoPiModel::PeakEnergy() const
{
  return massRho;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToTwoPiModel::ComputeCrossSection(G4double e) const
{
  G4double ee = std::min(GeV,e);
  return cross->CrossSection2pi(ee);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
