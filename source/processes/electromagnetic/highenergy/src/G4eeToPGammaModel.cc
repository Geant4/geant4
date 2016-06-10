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
// $Id: G4eeToPGammaModel.cc 91869 2015-08-07 15:21:02Z gcosmo $
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
// Creation date: 25.10.2003
//
// Modifications:
//
//
// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eeToPGammaModel.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PionZero.hh"
#include "G4Eta.hh"
#include "G4Gamma.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4eeCrossSections.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eeToPGammaModel::G4eeToPGammaModel(G4eeCrossSections* cr,
	 			     const G4String& npart,
				     G4double maxkinEnergy,
				     G4double binWidth)
: G4Vee2hadrons(cr,
		npart=="pi0" ? 782.62*MeV:1019.46*MeV,	
		maxkinEnergy,
		binWidth)
{
  G4cout << "####G4eeToPGammaModel & particle:" << npart 
	 << "####" << G4endl; 

  pi0 = G4PionZero::PionZero();
  if(npart == "pi0") {
    massR = 782.62*MeV;
    particle = pi0;
  } else {
    massR = 1019.46*MeV;
    particle = G4Eta::Eta();
  }
  massP = particle->GetPDGMass();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToPGammaModel::~G4eeToPGammaModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToPGammaModel::PeakEnergy() const
{
  return massR;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToPGammaModel::ComputeCrossSection(G4double e) const
{
  G4double xs;
  if(particle == pi0) xs = cross->CrossSectionPi0G(e);
  else                xs = cross->CrossSectionEtaG(e);
  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToPGammaModel::SampleSecondaries(std::vector<G4DynamicParticle*>* newp,
	    G4double e, const G4ThreeVector& direction)
{
  G4double egam = 0.5*e*(1.0 - massP*massP/(massR*massR));
  G4double tkin = e - egam - massP;
  if(tkin < 0.0) tkin = 0.0;
  G4double cost;
  do {
    cost = 2.0*G4UniformRand() - 1.0;
    // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  } while( 2.0*G4UniformRand() > 1.0 + cost*cost );

  G4double sint = sqrt(1.0 - cost*cost);
  G4double phi  = twopi * G4UniformRand();

  G4ThreeVector dir(sint*cos(phi),sint*sin(phi), cost);
  dir.rotateUz(direction);

  // create G4DynamicParticle objects
  G4DynamicParticle* p1 = 
     new G4DynamicParticle(particle,dir,tkin);
  G4DynamicParticle* p2 = 
     new G4DynamicParticle(G4Gamma::Gamma(),-dir,egam);
  newp->push_back(p1);
  newp->push_back(p2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

