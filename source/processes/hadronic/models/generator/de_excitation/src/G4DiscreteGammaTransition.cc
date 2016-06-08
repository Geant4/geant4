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
//      File name:     G4ContDiscrGammaTransition
//
//      Author:        Maria Grazia Pia   (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
//      
// -------------------------------------------------------------------

#include "G4DiscreteGammaTransition.hh"
#include "Randomize.hh"
#include "G4RandGeneralTmp.hh"


G4DiscreteGammaTransition::G4DiscreteGammaTransition(const G4NuclearLevel& level): 
  _level(level), _excitation(0.), _gammaEnergy(0.), _gammaCreationTime(0.)
{ }


G4DiscreteGammaTransition::~G4DiscreteGammaTransition() 
{ }


void G4DiscreteGammaTransition::SelectGamma()
{

  _gammaEnergy = 0.;

  G4int nGammas = _level.NumberOfGammas();
  if (nGammas > 0)
  {
    G4double random = G4UniformRand();

    G4int iGamma = 0;
    for(iGamma=0;iGamma < nGammas;iGamma++)
    {
      if(random <= _level.GammaCumulativeProbabilities().at(iGamma))
	break;
    }


    // Small correction due to the fact that there are mismatches between 
    // nominal level energies and emitted gamma energies

    G4double eCorrection = _level.Energy() - _excitation;

    _gammaEnergy = _level.GammaEnergies().at(iGamma) - eCorrection;

//  Warning: the following check is needed to avoid loops:
//  Due essentially to missing nuclear levels in data files, it is
//  possible that _gammaEnergy is so low as the nucleus doesn't change
//  its level after the transition.
//  When such case is found, force the full deexcitation of the nucleus.
//
//    NOTE: you should force the transition to the next lower level,
//          but this change needs a more complex revision of actual design.
//          I leave this for a later revision.

    if (_gammaEnergy < _level.Energy()*10e-5) _gammaEnergy = _excitation;

  }

  G4double tau = _level.HalfLife() / log(2.0);

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
  
  _gammaCreationTime = tMin + (tMax - tMin) * random;

//  if(_verbose > 10)
//    G4cout << "*---*---* G4DiscreteTransition: _gammaCreationTime = "
//	   << _gammaCreationTime/second << endl;
  return;
}

G4double G4DiscreteGammaTransition::GetGammaEnergy()
{
  return _gammaEnergy;
}

G4double G4DiscreteGammaTransition::GetGammaCreationTime()
{
  return _gammaCreationTime;
}

void G4DiscreteGammaTransition::SetEnergyFrom(const G4double energy)
{
  _excitation = energy;
  return;
}





