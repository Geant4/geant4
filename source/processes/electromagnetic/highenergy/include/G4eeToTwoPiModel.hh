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
// $Id: G4eeToTwoPiModel.hh,v 1.2 2006-06-29 19:32:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
