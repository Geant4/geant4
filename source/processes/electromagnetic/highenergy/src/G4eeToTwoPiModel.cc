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
// $Id: G4eeToTwoPiModel.cc,v 1.4 2006/06/29 19:32:50 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
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
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4eeCrossSections.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eeToTwoPiModel::G4eeToTwoPiModel(G4eeCrossSections* cr):
  cross(cr)
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToTwoPiModel::~G4eeToTwoPiModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToTwoPiModel::Initialise()
{
  massPi = G4PionPlus::PionPlus()->GetPDGMass();
  massRho = 770.*MeV;
  highEnergy = 1.*GeV;
  cross = new G4eeCrossSections();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4eeToTwoPiModel::PhysicsVector(G4double emin, 
                                                 G4double emax) const
{
  G4double tmin = max(emin, 2.0*massPi);
  G4double tmax = max(tmin, emax);
  G4int nbins = (G4int)((tmax - tmin)/(5.*MeV));
  G4PhysicsVector* v = new G4PhysicsLinearVector(emin,emax,nbins);
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

vector<G4DynamicParticle*>* G4eeToTwoPiModel::SampleSecondaries(
	    G4double e, const G4ThreeVector& direction) const
{
  vector<G4DynamicParticle*>* newp = new vector<G4DynamicParticle*>;
  G4double tkin = 0.5*e - massPi;
  if(tkin < 0.0) tkin = 0.0;
  G4double cost;
  do {
    cost = 2.0*G4UniformRand() - 1.0;
  } while( G4UniformRand() > 1.0 - cost*cost );

  G4double sint = sqrt(1.0 - cost*cost);
  G4double phi  = twopi * G4UniformRand();

  G4ThreeVector dir(sint*cos(phi),sint*sin(phi), cost);
  dir.rotateUz(direction);

  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* pip = 
     new G4DynamicParticle(G4PionPlus::PionPlus(),dir,tkin);
  G4DynamicParticle* pin = 
     new G4DynamicParticle(G4PionMinus::PionMinus(),-dir,tkin);
  newp->push_back(pip);
  newp->push_back(pin);
  return newp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

