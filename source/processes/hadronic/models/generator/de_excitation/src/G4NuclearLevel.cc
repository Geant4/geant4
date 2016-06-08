//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
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
  for (i=0; i<eGamma.size(); i++)
    {
      _energies.push_back(eGamma[i]);
      _weights.push_back(wGamma[i]);
      _polarities.push_back(polarities[i]);
    }
  _nGammas = _energies.size();
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
  for (i=0; i<_nGammas; i++) { G4cout << _energies[i] << " "; }
  G4cout << G4endl << "     Weights: ";
  for (i=0; i<_nGammas; i++) { G4cout << _weights[i] << " "; }
  G4cout << G4endl << "     Relative transition probabilities ";
  for (i=0; i<_nGammas; i++) { G4cout << _prob[i] << " "; }
  G4cout << G4endl << "     Cumulative probabilities: ";
  for (i=0; i<_nGammas; i++) { G4cout << _cumProb[i] << " "; }
  G4cout << G4endl << "     Polarities: ";
  for (i=0; i<_nGammas; i++) { G4cout << _polarities[i] << " "; }
  G4cout << G4endl;      

  return;
}
 

void G4NuclearLevel::MakeProbabilities()
{
  G4double sum = 0.;
  G4int i = 0;
  for (i=0; i<_nGammas; i++)
    {
      sum += _weights[i];
    }

  for (i=0; i<_nGammas; i++)
    {
      if (sum > 0.) { _prob.push_back(_weights[i] / sum); }
      else { _prob.push_back(1./_nGammas); }
    }
  return;
}


void G4NuclearLevel::MakeCumProb()
{
  if (_nGammas > 0)
    {
      G4double sum = _prob[0];
      _cumProb.push_back(sum);
      
      G4int i = 0;
      for (i=1; i<_nGammas; i++)
	{
	  sum += _prob[i];
          _cumProb.push_back(sum);
	}
    }
  return;
}
