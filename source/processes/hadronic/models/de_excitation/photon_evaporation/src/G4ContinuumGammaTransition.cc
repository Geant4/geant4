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
// $Id: G4ContinuumGammaTransition.cc 87443 2014-12-04 12:26:31Z gunter $
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4ContinuumGammaTransition
//
//      Authors:       Carlo Dallapiccola (dallapiccola@umdhep.umd.edu)
//                     Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
//        02 May 2003,   V. Ivanchenko change interface to G4NuclearlevelManager
//        06 Oct 2010,   M. Kelsey -- follow changes to G4NuclearLevelManager
//        17 Nov 2010,   V. Ivanchenko use exponential law for sampling of time
//                                     and extra cleanup
// ----------------------------------------------------------------------------
//
//  Class G4ContinuumGammaTransition.cc
//

#include "G4ContinuumGammaTransition.hh"
#include "G4NuclearLevelManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Pow.hh"
#include "G4Log.hh"

//static const G4double normC = CLHEP::millibarn/
//  (CLHEP::pi*CLHEP::hbarc*CLHEP::pi*CLHEP::hbarc);
static const G4double factor = 5;
static const G4double tolerance = 2*CLHEP::keV;

G4ContinuumGammaTransition::G4ContinuumGammaTransition(
                            const G4NuclearLevelManager* manager,
			    G4int Z, G4int A, G4double exc, G4int verb)
{
  nBins = 100;
  verbose = verb;
  eMin = keV;
  eMax = minLevelE = eGamma = gammaCreationTime = 0.0;
  sampleArray.resize(nBins+1,0.0);
  g4pow = G4Pow::GetInstance();
  Update(manager, Z, A, exc);
}

G4ContinuumGammaTransition::~G4ContinuumGammaTransition() 
{}

void G4ContinuumGammaTransition::Update(const G4NuclearLevelManager* manager,
					G4int Z, G4int A, G4double exc) 
{
  levelManager = manager;
  nucleusZ = Z;
  nucleusA = A;
  excitation = exc;
  eGamma = 0.;
  gammaCreationTime = 0.;

  minLevelE = DBL_MAX;
  if(levelManager) { 
    minLevelE = std::max(0.5*levelManager->MinLevelEnergy(), tolerance);
  }
  // Energy range for photon generation; 
  // upper limit is defined 5*Gamma(GDR) from GDR peak
  // Giant Dipole Resonance energy
  G4double energyGDR = (40.3 / g4pow->powZ(nucleusA,0.2) ) * MeV;
  energyGDR2 = energyGDR*energyGDR;
  // Giant Dipole Resonance width
  widthGDR = 0.30 * energyGDR;
  widthGDR2 = widthGDR*widthGDR;
  // Extend 
  eMax = energyGDR + factor * widthGDR;
  if (eMax > excitation) { eMax = excitation; }
}

void G4ContinuumGammaTransition::SelectGamma()
{
  eGamma = 0.;
  sampleArray[0] = 0.0;
  G4int i;
  G4double del = (eMax - eMin) / G4double(nBins);
  G4double sum = 0;
  G4double w1 = E1Pdf(eMin);
  G4double w2;
  //G4cout << eMin << "  " << eMax << "  " << del << G4endl;
  for (i=1; i<=nBins; i++) {
    G4double e = eMin + del * i;
    w2 = E1Pdf(e);
    sum += 0.5*(w1 + w2);
    w1 = w2;
    sampleArray[i] = sum;
    if(verbose > 1) {
      G4cout << "*---* G4ContinuumTransition: e = " << e 
	     << " pdf = " << sampleArray[i] << G4endl;
    }
  }
  sum *= G4UniformRand();
  eGamma = eMax;
  for (i=1; i<=nBins; i++) {
    if(sum <= sampleArray[i]) {
      eGamma = eMin + del * i;
      G4double w = sampleArray[i] - sampleArray[i-1];
      //G4cout << eGamma << "  " << w << G4endl;
      if(w != 0.0) {
	eGamma -= (sampleArray[i] - sum)*del/w;
      }
      break;
    }
  }   
  
  G4double finalExcitation = excitation - eGamma;
  
  if(verbose > 1) {
    G4cout << "*---*---* G4ContinuumTransition: eGamma = " << eGamma
	   << "   finalExcitation = " << finalExcitation << G4endl;
  }

  if(finalExcitation <= minLevelE) {
    eGamma = excitation;

  } else {
    finalExcitation = levelManager->NearestLevel(finalExcitation)->Energy();
    eGamma = excitation - finalExcitation;
  }  

  gammaCreationTime = GammaTime();

  if(verbose > 1) {
    G4cout << "*---*---* G4ContinuumTransition: gammaCreationTime = "
	   << gammaCreationTime/second << G4endl;
  }
}

G4double G4ContinuumGammaTransition::GetGammaEnergy()
{
  return eGamma;
}

G4double G4ContinuumGammaTransition::GetGammaCreationTime()
{
  return gammaCreationTime;
}

void G4ContinuumGammaTransition::SetEnergyFrom(G4double energy)
{
  excitation = energy;
}

G4double G4ContinuumGammaTransition::E1Pdf(G4double e)
{
  G4double theProb = 0.0;
  G4double U = excitation - e;

  if(U < 0.0) { return theProb; }

  G4double aLevelDensityParam = 
    ldPar.LevelDensityParameter(nucleusA,nucleusZ,excitation);

  //G4double levelDensBef = G4Exp(2.0*std::sqrt(aLevelDensityParam*excitation));
  //G4double levelDensAft = G4Exp(2.0*std::sqrt(aLevelDensityParam*(excitation - e)));
  G4double coeff = G4Exp(2.0*(std::sqrt(aLevelDensityParam*U) 
			      - std::sqrt(aLevelDensityParam*excitation)));

  //if(verbose > 20)
  //  G4cout << nucleusA << " LevelDensityParameter = " <<  aLevelDensityParam
  //	   << " Bef Aft " << levelDensBef << " " << levelDensAft << G4endl;
  
  // Now form the probability density
  // Define constants for the photoabsorption cross-section (the reverse
  // process of our de-excitation)

  G4double sigma0 = 2.5 * nucleusA;  
 
  G4double e2 = e*e;
  G4double numerator = sigma0 * e2 * widthGDR2;
  G4double denominator = (e2 - energyGDR2)* (e2 - energyGDR2) + widthGDR2*e2;
  //  if (denominator < 1.0e-9) denominator = 1.0e-9;

  G4double sigmaAbs = numerator/denominator; 

  if(verbose > 2) {
    G4cout << "E_GDR(MeV)= " << std::sqrt(energyGDR2) << " W_GDR(MeV)= " << widthGDR 
	   << " sigAbs= " << sigmaAbs  
	   << " E(MeV)= " << e << " coeff= " << coeff
	   << G4endl;
  }

  theProb = sigmaAbs * e2 * coeff;

  return theProb;
}

G4double G4ContinuumGammaTransition::GammaTime()
{
  G4double tau = hbar_Planck/widthGDR;
  G4double creationTime = -tau*G4Log(G4UniformRand());

  return creationTime;
}












