// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
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
// -------------------------------------------------------------------

#include "G4DiscreteGammaTransition.hh"
#include "Randomize.hh"


G4DiscreteGammaTransition::G4DiscreteGammaTransition(const G4NuclearLevel& level): 
  _level(level), _excitation(0.), _gammaEnergy(0.)
{ }


G4DiscreteGammaTransition::~G4DiscreteGammaTransition() 
{ }


G4double G4DiscreteGammaTransition::GammaEnergy()
{

  _gammaEnergy = 0.;

      G4int nGammas = _level.NumberOfGammas();
      if (nGammas > 0)
	{
	  G4double random = G4UniformRand();
	  G4int iGamma = 0;
	  if (random <= _level.GammaCumulativeProbabilities().at(0)) iGamma = 0;
	  else
	    {
	      G4int i;
	      
	      for (i=1; i<nGammas; i++)
		{
		  if (random > _level.GammaCumulativeProbabilities().at(i-1) &&
		      random <= _level.GammaCumulativeProbabilities().at(i)) 
		    { iGamma = i; }
		}
	    }

	  // Small correction due to the fact that there are mismatches between 
	  // nominal level energies and emitted gamma energies
	  G4double eCorrection = _level.Energy() - _excitation;

	  _gammaEnergy = _level.GammaEnergies().at(iGamma) - eCorrection;
          if (_gammaEnergy < 0.) _gammaEnergy = 0.;
	}

  return _gammaEnergy;
}


G4double G4DiscreteGammaTransition::GetEnergyTo() const
{
  G4double energyTo = _excitation - _gammaEnergy;
  if (energyTo < 0.) energyTo = 0.;
  
  return energyTo;
}


void G4DiscreteGammaTransition::SetEnergyFrom(const G4double energy)
{
  _excitation = energy;
  return;
}
