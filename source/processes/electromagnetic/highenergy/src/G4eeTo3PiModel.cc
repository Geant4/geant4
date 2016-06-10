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
// $Id: G4eeTo3PiModel.cc 91869 2015-08-07 15:21:02Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eeTo3PiModel
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

#include "G4eeTo3PiModel.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4eeCrossSections.hh"
#include "G4RandomDirection.hh"
#include <complex>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eeTo3PiModel::G4eeTo3PiModel(G4eeCrossSections* cr,
			       G4double maxkinEnergy,
			       G4double binWidth)
:  G4Vee2hadrons(cr,
	         0.41612*GeV,	 //threshold
		 maxkinEnergy,
		 binWidth)
{
  G4cout << "####G4eeTo3PiModel####" << G4endl;

  massPi  = G4PionPlus::PionPlus()->GetPDGMass();
  massPi0 = G4PionZero::PionZero()->GetPDGMass();
  massOm  = 782.62*MeV;
  massPhi = 1019.46*MeV;
  gmax    = 3.0e-8;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeTo3PiModel::~G4eeTo3PiModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeTo3PiModel::PeakEnergy() const
{
  G4double e = massOm;
  if(HighEnergy() > massPhi) { e = massPhi; } 
  return e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeTo3PiModel::ComputeCrossSection(G4double e) const
{
  return cross->CrossSection3pi(e);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeTo3PiModel::SampleSecondaries(std::vector<G4DynamicParticle*>* newp,
	    G4double e, const G4ThreeVector& direction) 
{

  G4double x0 = massPi0/e;
  G4double x1 = massPi/e;

  G4LorentzVector w0(0.,0.,0.,0.), w1(0.,0.,0.,0.), w2(0.,0.,0.,0.);
  G4ThreeVector dir0, dir1;
  G4double e0, p0, e2, p, gg, m01, m02, m12;

  // max pi0 energy
  G4double edel  = 0.5*e*(1.0 + x0*x0 - 4.0*x1*x1) - massPi0;

  const G4int nmax = 200;
  G4int nn = 0;
  do {
    ++nn;
    // pi0 sample
    e0 = edel*G4UniformRand() + massPi0;
    p0 = sqrt(e0*e0 - massPi0*massPi0);
    dir0 = G4RandomDirection();
    w0 = G4LorentzVector(p0*dir0.x(),p0*dir0.y(),p0*dir0.z(),e0);

    // pi+pi- pair
    w1 = G4LorentzVector(-p0*dir0.x(),-p0*dir0.y(),-p0*dir0.z(),e-e0);
    G4ThreeVector bst = w1.boostVector();
    e2 = 0.25*w1.m2();

    // pi+ 
    p = sqrt(e2 - massPi*massPi);
    dir1 = G4RandomDirection();
    w2 = G4LorentzVector(p*dir1.x(),p*dir1.y(),p*dir1.z(),sqrt(e2));
    // pi- 
    w1.set(-w2.px(), -w2.py(), -w2.pz(), w2.e());

    w1.boost(bst);
    w2.boost(bst);

    G4double px2 = w2.x();
    G4double py2 = w2.y();
    G4double pz2 = w2.z();

    G4double px1 = w1.x();
    G4double py1 = w1.y();
    G4double pz1 = w1.z();

    m01 = w0*w1;
    m02 = w0*w2;
    m12 = w1*w2;

    G4double px = py1*pz2 - py2*pz1;
    G4double py = pz1*px2 - pz2*px1;
    G4double pz = px1*py2 - px2*py1;

    gg = (px*px + py*py + pz*pz)*
      norm( 1.0/cross->DpRho(m01) +  1.0/cross->DpRho(m02)
	    + 1.0/cross->DpRho(m12) );

    if(gg > gmax) {
      G4cout << "G4eeTo3PiModel::SampleSecondaries WARNING matrix element g= "
	     << gg << " > " << gmax << " (majoranta)" << G4endl;
      gmax = gg;
    }
    // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  } while( gmax*G4UniformRand() > gg || nn < nmax);

  w0.rotateUz(direction);
  w1.rotateUz(direction);
  w2.rotateUz(direction);

  // create G4DynamicParticle objects 
  G4DynamicParticle* dp0 = 
     new G4DynamicParticle(G4PionZero::PionZero(), w0);
  G4DynamicParticle* dp1 = 
     new G4DynamicParticle(G4PionPlus::PionPlus(), w1);
  G4DynamicParticle* dp2 = 
     new G4DynamicParticle(G4PionMinus::PionMinus(), w2);
  newp->push_back(dp0);
  newp->push_back(dp1);
  newp->push_back(dp2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

