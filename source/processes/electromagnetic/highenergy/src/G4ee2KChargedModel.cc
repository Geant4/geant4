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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4ee2KChargedModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 09.07.2008
//
// Modifications:
//
//
// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ee2KChargedModel.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4eeCrossSections.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4ee2KChargedModel::G4ee2KChargedModel(G4eeCrossSections* cr,
				       G4double maxkinEnergy,
				       G4double binWidth)
:  G4Vee2hadrons(cr,				      
		 2.0*G4KaonPlus::KaonPlus()->GetPDGMass(),
		 maxkinEnergy,
		 binWidth)
{
  G4cout << "####G4ee2KChargedModel####" << G4endl;

  massK = G4KaonPlus::KaonPlus()->GetPDGMass();
  massPhi = 1019.46*MeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ee2KChargedModel::~G4ee2KChargedModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ee2KChargedModel::PeakEnergy() const
{
  return massPhi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ee2KChargedModel::ComputeCrossSection(G4double e) const
{
  return cross->CrossSection2Kcharged(e); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ee2KChargedModel::SampleSecondaries(std::vector<G4DynamicParticle*>* newp,
	    G4double e, const G4ThreeVector& direction)
{

  G4double tkin = 0.5*e - massK;
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
  G4DynamicParticle* p1 = 
     new G4DynamicParticle(G4KaonPlus::KaonPlus(),dir,tkin);
  G4DynamicParticle* p2 = 
     new G4DynamicParticle(G4KaonMinus::KaonMinus(),-dir,tkin);
  newp->push_back(p1);
  newp->push_back(p2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

