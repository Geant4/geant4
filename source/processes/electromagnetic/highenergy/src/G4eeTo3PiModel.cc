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
// $Id: G4eeTo3PiModel.cc,v 1.1 2008/07/10 18:07:27 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
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

G4eeTo3PiModel::G4eeTo3PiModel(G4eeCrossSections* cr):
  cross(cr)
{
  massPi  = G4PionPlus::PionPlus()->GetPDGMass();
  massPi0 = G4PionZero::PionZero()->GetPDGMass();
  massOm  = 782.62*MeV;
  massPhi = 1019.46*MeV;
  gcash   = 0.0;
  gmax    = 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeTo3PiModel::~G4eeTo3PiModel()
{
  G4cout << "### G4eeTo3PiModel::~G4eeTo3PiModel: gmax= "
	 << gmax << " gcash= " << gcash << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4eeTo3PiModel::PhysicsVector(G4double emin, 
					       G4double emax) const
{
  G4double tmin = std::max(emin, ThresholdEnergy());
  G4double tmax = std::max(tmin, emax);
  G4int nbins = (G4int)((tmax - tmin)/(1.*MeV));
  G4PhysicsVector* v = new G4PhysicsLinearVector(emin,emax,nbins);
  v->SetSpline(true);
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeTo3PiModel::SampleSecondaries(std::vector<G4DynamicParticle*>* newp,
	    G4double e, const G4ThreeVector& direction) 
{
  if(e < ThresholdEnergy()) return;

  G4double x0 = massPi0/e;
  G4double x1 = massPi/e;

  G4LorentzVector w0, w1, w2;
  G4ThreeVector dir0, dir1, mom, mom1, mom2;
  G4double e0, p0, e2, p, g, m01, m02, m12;

  // max pi0 energy
  G4double edel  = 0.5*e*(1.0 + x0*x0 - 4.0*x1*x1) - massPi0;

  do {
    // pi0 sample
    e0 = edel*G4UniformRand() + massPi0;
    p0 = sqrt(e0 - massPi0*massPi0);
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
    w2.boost(bst);
    mom2 = w2.vect();

    // pi- 
    w1 -= w2;
    mom1 = w2.vect();

    m01 = w0*w1;
    m02 = w0*w2;
    m12 = w1*w2;

    mom = mom1*mom2;
    g = mom.mag2()*norm( 1.0/cross->DpRho(m01) +  1.0/cross->DpRho(m02)
			 + 1.0/cross->DpRho(m12) );
    if(g > gmax) {
      G4cout << "G4eeTo3PiModel::SampleSecondaries WARNING matrix element g= "
	     << g << " > " << gmax << " (majoranta)" << G4endl;
    }
    if(g > gcash) gcash = g;
    
  } while( gmax*G4UniformRand() > g );

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

