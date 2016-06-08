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
//      File name:     G4NuclearLevel
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 24 October 1998
//
//      Modifications: 
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data. 
//      
// -------------------------------------------------------------------

#include "G4NuclearLevel.hh"

#include "globals.hh"

G4NuclearLevel::G4NuclearLevel(const G4double energy, const G4double halfLife,
			       const G4double angularMomentum,
			       const G4DataVector& eGamma,
			       const G4DataVector& wGamma,
			       const G4DataVector& polarities)
{
  _energy = energy;
  _halfLife = halfLife;
  _angularMomentum = angularMomentum;
  G4int i;
  for (i=0; i<eGamma.entries(); i++)
    {
      _energies.insert(eGamma.at(i));
      _weights.insert(wGamma.at(i));
      _polarities.insert(polarities.at(i));
    }
  _nGammas = _energies.entries();
  MakeProbabilities();
  MakeCumProb();
}

G4NuclearLevel::~G4NuclearLevel()
{ }


G4bool G4NuclearLevel::operator==(const G4NuclearLevel &right) const
{
  return (this == (G4NuclearLevel *) &right);
}


G4bool G4NuclearLevel::operator!=(const G4NuclearLevel &right) const
{
  return (this != (G4NuclearLevel *) &right);
}


G4bool G4NuclearLevel::operator<(const G4NuclearLevel &right) const  
{
  if (_energy < right.Energy()) return true;
  else return false;
}


const G4DataVector& G4NuclearLevel::GammaEnergies() const
{
  return _energies;
}
 
const G4DataVector& G4NuclearLevel::GammaWeights() const
{
  return _weights;
}
 

const G4DataVector& G4NuclearLevel::GammaProbabilities() const
{
  return _prob;
}
 

const G4DataVector& G4NuclearLevel::GammaCumulativeProbabilities() const
{
  return _cumProb;
}
 

const G4DataVector& G4NuclearLevel::GammaPolarities() const
{
  return _polarities;
}
 
G4double G4NuclearLevel::Energy() const
{
  return _energy;
}
 
G4double G4NuclearLevel::AngularMomentum() const
{
  return _angularMomentum;
}
 
G4double G4NuclearLevel::HalfLife() const
{
  return _halfLife;
}
 

G4int G4NuclearLevel::NumberOfGammas() const
{
  return _nGammas;
}
 

void G4NuclearLevel::PrintAll() const 
{
  G4cout << "---- Level energy = " << _energy << ", angular momentum = "
	 << _angularMomentum << ", half life " << _halfLife
	 << ", " << _nGammas << " photons" << G4endl;
  G4int i;
  G4cout << "     Gammas: ";
  for (i=0; i<_nGammas; i++) { G4cout << _energies.at(i) << " "; }
  G4cout << G4endl << "     Weights: ";
  for (i=0; i<_nGammas; i++) { G4cout << _weights.at(i) << " "; }
  G4cout << G4endl << "     Relative transition probabilities ";
  for (i=0; i<_nGammas; i++) { G4cout << _prob.at(i) << " "; }
  G4cout << G4endl << "     Cumulative probabilities: ";
  for (i=0; i<_nGammas; i++) { G4cout << _cumProb.at(i) << " "; }
  G4cout << G4endl << "     Polarities: ";
  for (i=0; i<_nGammas; i++) { G4cout << _polarities.at(i) << " "; }
  G4cout << G4endl;      

  return;
}
 

void G4NuclearLevel::MakeProbabilities()
{
  G4double sum = 0.;
  G4int i = 0;
  for (i=0; i<_nGammas; i++)
    {
      sum += _weights.at(i);
    }

  for (i=0; i<_nGammas; i++)
    {
      if (sum > 0.) { _prob.insert(_weights.at(i) / sum); }
      else { _prob.insert(1./_nGammas); }
    }
  return;
}


void G4NuclearLevel::MakeCumProb()
{
  if (_nGammas > 0)
    {
      G4double sum = _prob.at(0);
      _cumProb.insert(sum);
      
      G4int i = 0;
      for (i=1; i<_nGammas; i++)
	{
	  sum += _prob.at(i);
          _cumProb.insert(sum);
	}
    }
  return;
}
