// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
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
//      
// -------------------------------------------------------------------
//
//  Class G4ContinuumGammaTransition.cc
//

#include "G4ContinuumGammaTransition.hh"
#include "G4VLevelDensityParameter.hh"
#include "G4ConstantLevelDensityParameter.hh"
#include "G4RandGeneralTmp.hh"
//
// Constructor
//

G4ContinuumGammaTransition::G4ContinuumGammaTransition(
                            const G4NuclearLevelManager& levelManager,
			    G4int Z, G4int A,
			    G4double excitation,
			    G4int verbose):
  _nucleusZ(Z), _nucleusA(A), _excitation(excitation), _levelManager(levelManager) 
{
  const G4PtrLevelVector* levels = levelManager.GetLevels();
  G4double eTolerance = 0.;
  if (levels != 0)
  {
    G4int lastButOne = levelManager.NumberOfLevels() - 2;
    if (lastButOne >= 0)
    {
      eTolerance = levelManager.MaxLevelEnergy() - levels->at(lastButOne)->Energy();
      if (eTolerance < 0.) eTolerance = 0.;
    }
  }

  _verbose = verbose;
  _eGamma = 0.;
  _gammaCreationTime = 0.;

  _maxLevelE = levelManager.MaxLevelEnergy() + eTolerance;
  _minLevelE = levelManager.MinLevelEnergy();

  // Energy range for photon generation; upper limit is defined 5*Gamma(GDR) from GDR peak
  _eMin = 0.001 * MeV;
  // Giant Dipole Resonance energy
  G4double energyGDR = (40.3 / pow(_nucleusA,0.2) ) * MeV;
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

G4ContinuumGammaTransition::~G4ContinuumGammaTransition() {}

//
// Override GammaEnergy function from G4VGammaTransition
//

void G4ContinuumGammaTransition::SelectGamma()
{

  _eGamma = 0.;

  G4int nBins = 200;
  G4double sampleArray[200];
  G4int i;
  for (i=0; i<nBins; i++)
    {
      G4double e = _eMin + ( (_eMax - _eMin) / nBins) * i;
      sampleArray[i] = E1Pdf(e);

      if(_verbose > 10)
	G4cout << "*---* G4ContinuumTransition: e = " << e 
	       << " pdf = " << sampleArray[i] << endl;
    }
  G4RandGeneralTmp randGeneral(sampleArray, nBins);
  G4double random = randGeneral.shoot();
  
  _eGamma = _eMin + (_eMax - _eMin) * random;
  
  G4double finalExcitation = _excitation - _eGamma;
  
  if(_verbose > 10)
    G4cout << "*---*---* G4ContinuumTransition: eGamma = " << _eGamma
	   << "   finalExcitation = " << finalExcitation 
	   << " random = " << random << endl;

//  if (finalExcitation < 0)
  if(finalExcitation < _minLevelE/2.)
    {
      _eGamma = _excitation;
      finalExcitation = 0.;
    }
  
  if (finalExcitation < _maxLevelE && finalExcitation > 0.) 
    {
      G4double levelE = _levelManager.NearestLevel(finalExcitation)->Energy();
      G4double diff = finalExcitation - levelE;
      _eGamma = _eGamma + diff;
    }  

  _gammaCreationTime = GammaTime();

  if(_verbose > 10)
    G4cout << "*---*---* G4ContinuumTransition: _gammaCreationTime = "
	   << _gammaCreationTime/second << endl;

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


void G4ContinuumGammaTransition::SetEnergyFrom(const G4double energy)
{

  if (energy > 0.) _excitation = energy;
  return;  

}


G4double G4ContinuumGammaTransition::E1Pdf(G4double e)
{
  G4double theProb = 0.0;

  if( (_excitation - e) < 0.0 || e < 0 || _excitation < 0) return theProb;

  G4ConstantLevelDensityParameter ldPar;
  G4double aLevelDensityParam = ldPar.LevelDensityParameter(_nucleusA,_nucleusZ,_excitation);

  G4double levelDensBef = exp(2.0*sqrt(aLevelDensityParam*_excitation));
  G4double levelDensAft = exp(2.0*sqrt(aLevelDensityParam*(_excitation - e)));

  if(_verbose > 20)
    G4cout << _nucleusA << " LevelDensityParameter = " <<  aLevelDensityParam
	   << " Bef Aft " << levelDensBef << " " << levelDensAft << endl;
  
  // Now form the probability density

  // Define constants for the photoabsorption cross-section (the reverse
  // process of our de-excitation)

  //  G4double sigma0 = 2.5 * _nucleusA * millibarn;  
  G4double sigma0 = 2.5 * _nucleusA;  

  G4double Egdp = (40.3 / pow(_nucleusA,0.2) )*MeV;
  G4double GammaR = 0.30 * Egdp;
 
  G4double normC = 1.0 / (pi * hbarc)*(pi * hbarc);

  G4double numerator = sigma0 * e*e * GammaR*GammaR;
  G4double denominator = (e*e - Egdp*Egdp)* (e*e - Egdp*Egdp) + GammaR*GammaR*e*e;
  //  if (denominator < 1.0e-9) denominator = 1.0e-9;

  G4double sigmaAbs = numerator/denominator ; 

  if(_verbose > 20)
    G4cout << ".. " << Egdp << " .. " << GammaR 
	   << " .. " << normC << " .. " << sigmaAbs  
	   << " .. " << e*e << " .. " << levelDensAft/levelDensBef
	   << endl;

  //  theProb = normC * sigmaAbs * e*e * levelDensAft/levelDensBef;
  theProb =  sigmaAbs * e*e * levelDensAft/levelDensBef;

  return theProb;
}


G4double G4ContinuumGammaTransition::GammaTime()
{

  G4double GammaR = 0.30 * (40.3 / pow(_nucleusA,0.2) )*MeV;
  G4double tau = hbar_Planck/GammaR;

  G4double tMin = 0;
  G4double tMax = 10.0 * tau;
  G4int nBins = 200;
  G4double sampleArray[200];

  for(G4int i = 0;i<nBins;i++)
  {
    G4double t = tMin + ((tMax-tMin)/nBins)*i;
    sampleArray[i] = (exp(-t/tau))/tau;
  }

  G4RandGeneralTmp randGeneral(sampleArray, nBins);
  G4double random = randGeneral.shoot();
  
  G4double creationTime = tMin + (tMax - tMin) * random;

  return creationTime;
}












