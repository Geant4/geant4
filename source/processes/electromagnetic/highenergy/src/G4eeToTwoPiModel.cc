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
// $Id: G4eeToTwoPiModel.cc 91869 2015-08-07 15:21:02Z gcosmo $
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
// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eeToTwoPiModel.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4eeCrossSections.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eeToTwoPiModel::G4eeToTwoPiModel(G4eeCrossSections* cr,
				   G4double maxkinEnergy,
				   G4double binWidth)
:   G4Vee2hadrons(cr,
		  2.0*G4PionPlus::PionPlus()->GetPDGMass(),
                  maxkinEnergy,
		  binWidth)
{
  G4cout << "#####G4eeToTwoPiModel####" << G4endl;  

  massPi = G4PionPlus::PionPlus()->GetPDGMass();
  massRho = 775.5*MeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToTwoPiModel::~G4eeToTwoPiModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToTwoPiModel::PeakEnergy() const
{
  return massRho;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToTwoPiModel::ComputeCrossSection(G4double e) const
{
  return cross->CrossSection2pi(e);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToTwoPiModel::SampleSecondaries(std::vector<G4DynamicParticle*>* newp,
	    G4double e, const G4ThreeVector& direction)
{

  G4double tkin = 0.5*e - massPi;
  if(tkin < 0.0) tkin = 0.0;
  G4double cost;
  do {
    cost = 2.0*G4UniformRand() - 1.0;
    // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  } while( G4UniformRand() > 1.0 - cost*cost );

  G4double sint = sqrt(1.0 - cost*cost);
  G4double phi  = twopi * G4UniformRand();

  G4ThreeVector dir(sint*cos(phi),sint*sin(phi), cost);
  dir.rotateUz(direction);

  // create G4DynamicParticle objects
  G4DynamicParticle* pip = 
     new G4DynamicParticle(G4PionPlus::PionPlus(),dir,tkin);
  G4DynamicParticle* pin = 
     new G4DynamicParticle(G4PionMinus::PionMinus(),-dir,tkin);
  newp->push_back(pip);
  newp->push_back(pin);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

