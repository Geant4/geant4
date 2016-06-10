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
// $Id: G4ContinuumGammaTransition.cc 74869 2013-10-23 09:26:17Z gcosmo $
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
#include "G4VLevelDensityParameter.hh"
#include "G4ConstantLevelDensityParameter.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Pow.hh"
#include "G4Log.hh"

//
// Constructor
//

G4ContinuumGammaTransition::G4ContinuumGammaTransition(
                            const G4NuclearLevelManager* levelManager,
			    G4int Z, G4int A,
			    G4double excitation,
			    G4int verbose):
  _nucleusA(A), _nucleusZ(Z), _excitation(excitation), _levelManager(levelManager) 
{
  G4double eTolerance = 0.;
  G4int lastButOne = _levelManager->NumberOfLevels() - 2;
  if (lastButOne >= 0)
    {
      eTolerance = (_levelManager->MaxLevelEnergy() -
		    _levelManager->GetLevel(lastButOne)->Energy());
      if (eTolerance < 0.) eTolerance = 0.;
    }
  

  _verbose = verbose;
  _eGamma = 0.;
  _gammaCreationTime = 0.;

  _maxLevelE = _levelManager->MaxLevelEnergy() + eTolerance;
  _minLevelE = _levelManager->MinLevelEnergy();

  // Energy range for photon generation; upper limit is defined 5*Gamma(GDR) from GDR peak
  _eMin = 0.001 * MeV;
  // Giant Dipole Resonance energy
  G4double energyGDR = (40.3 / G4Pow::GetInstance()->powZ(_nucleusA,0.2) ) * MeV;
  // Giant Dipole Resonance width
  G4double widthGDR = 0.30 * energyGDR;
  // Extend 
  G4double factor = 5;
  _eMax = energyGDR + factor * widthGDR;
  if (_eMax > excitation) _eMax = _excitation;

}

//
// Destructor
//

G4ContinuumGammaTransition::~G4ContinuumGammaTransition() 
{}

void G4ContinuumGammaTransition::SelectGamma()
{

  _eGamma = 0.;

  G4int nBins = 100;
  G4double sampleArray[101];
  sampleArray[0] = 0.0;
  G4int i;
  G4double del = (_eMax - _eMin) / G4double(nBins);
  G4double sum = 0;
  G4double w1 = E1Pdf(_eMin);
  G4double w2;
  //G4cout << _eMin << "  " << _eMax << "  " << del << G4endl;
  for (i=1; i<=nBins; i++) {
    G4double e = _eMin + del * i;
    w2 = E1Pdf(e);
    sum += 0.5*(w1 + w2);
    w1 = w2;
    sampleArray[i] = sum;
    if(_verbose > 10) {
      G4cout << "*---* G4ContinuumTransition: e = " << e 
	     << " pdf = " << sampleArray[i] << G4endl;
    }
  }
  sum *= G4UniformRand();
  _eGamma = _eMax;
  for (i=1; i<=nBins; i++) {
    if(sum <= sampleArray[i]) {
      _eGamma = _eMin + del * i;
      G4double w = sampleArray[i] - sampleArray[i-1];
      //G4cout << _eGamma << "  " << w << G4endl;
      if(w != 0.0) {
	_eGamma -= (sampleArray[i] - sum)*del/w;
      }
      break;
    }
  }   
  
  G4double finalExcitation = _excitation - _eGamma;
  
  if(_verbose > 10) {
    G4cout << "*---*---* G4ContinuumTransition: eGamma = " << _eGamma
	   << "   finalExcitation = " << finalExcitation << G4endl;
  }
  //  if (finalExcitation < 0)
  if(finalExcitation < _minLevelE/2.)
    {
      _eGamma = _excitation;
      finalExcitation = 0.;
    }
  
  if (finalExcitation < _maxLevelE && finalExcitation > 0.) 
    {
      G4double levelE = _levelManager->NearestLevel(finalExcitation)->Energy();
      G4double diff = finalExcitation - levelE;
      _eGamma = _eGamma + diff;
    }  

  _gammaCreationTime = GammaTime();

  if(_verbose > 10) {
    G4cout << "*---*---* G4ContinuumTransition: _gammaCreationTime = "
	   << _gammaCreationTime/second << G4endl;
  }
  return;  
}

G4double G4ContinuumGammaTransition::GetGammaEnergy()
{
  return _eGamma;
}

G4double G4ContinuumGammaTransition::GetGammaCreationTime()
{
  return _gammaCreationTime;
}


void G4ContinuumGammaTransition::SetEnergyFrom(G4double energy)
{
  if (energy > 0.) _excitation = energy;
}


G4double G4ContinuumGammaTransition::E1Pdf(G4double e)
{
  G4double theProb = 0.0;
  G4double U = std::max(0.0, _excitation - e);

  if(e < 0.0 || _excitation < 0.0) { return theProb; }

  G4ConstantLevelDensityParameter ldPar;
  G4double aLevelDensityParam = 
    ldPar.LevelDensityParameter(_nucleusA,_nucleusZ,_excitation);

  //G4double levelDensBef = std::exp(2.0*std::sqrt(aLevelDensityParam*_excitation));
  //G4double levelDensAft = std::exp(2.0*std::sqrt(aLevelDensityParam*(_excitation - e)));
  G4double coeff = std::exp(2.0*(std::sqrt(aLevelDensityParam*U) 
				 - std::sqrt(aLevelDensityParam*_excitation)));

  //if(_verbose > 20)
  //  G4cout << _nucleusA << " LevelDensityParameter = " <<  aLevelDensityParam
  //	   << " Bef Aft " << levelDensBef << " " << levelDensAft << G4endl;
  
  // Now form the probability density

  // Define constants for the photoabsorption cross-section (the reverse
  // process of our de-excitation)

  //  G4double sigma0 = 2.5 * _nucleusA * millibarn;  
  G4double sigma0 = 2.5 * _nucleusA;  

  G4double Egdp = (40.3 /G4Pow::GetInstance()->powZ(_nucleusA,0.2) )*MeV;
  G4double GammaR = 0.30 * Egdp;
 
  static const G4double normC = 1.0 / (pi * hbarc)*(pi * hbarc);

  G4double numerator = sigma0 * e*e * GammaR*GammaR;
  G4double denominator = (e*e - Egdp*Egdp)* (e*e - Egdp*Egdp) + GammaR*GammaR*e*e;
  //  if (denominator < 1.0e-9) denominator = 1.0e-9;

  G4double sigmaAbs = numerator/denominator ; 

  if(_verbose > 20) {
    G4cout << ".. " << Egdp << " .. " << GammaR 
	   << " .. " << normC << " .. " << sigmaAbs  
	   << " .. " << e*e << " .. " << coeff
	   << G4endl;
  }

  //  theProb = normC * sigmaAbs * e*e * levelDensAft/levelDensBef;
  theProb =  sigmaAbs * e*e * coeff;

  return theProb;
}


G4double G4ContinuumGammaTransition::GammaTime()
{

  G4double GammaR = 0.30 * (40.3 /G4Pow::GetInstance()->powZ(_nucleusA,0.2) )*MeV;
  G4double tau = hbar_Planck/GammaR;
  G4double creationTime = -tau*G4Log(G4UniformRand());

  return creationTime;
}












