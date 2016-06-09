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
//        09 Sep. 2002, Fan Lei  (flei@space.qinetiq.com)
//              Added IC probability when calculate the channel probabilities in 
//              MakeProbabilities().
//
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              Added K->N+ internal  conversion coefficiencies and their access
//              functions.
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data. 
//      
// -------------------------------------------------------------------

#include "G4NuclearLevel.hh"

#include "globals.hh"

G4int G4NuclearLevel::Increment(G4int aF)
{
  static G4int instanceCount = 0;
  instanceCount+=aF;
  return instanceCount;
}

G4NuclearLevel::G4NuclearLevel(const G4double energy, const G4double halfLife,
			       const G4double angularMomentum,
			       const std::vector<double>& eGamma,
			       const std::vector<double>& wGamma,
			       const std::vector<double>& polarities,
			       const std::vector<double>& kCC, const std::vector<double>& l1CC,
			       const std::vector<double>& l2CC, const std::vector<double>& l3CC,
			       const std::vector<double>& m1CC, const std::vector<double>& m2CC,
			       const std::vector<double>& m3CC, const std::vector<double>& m4CC,
			       const std::vector<double>& m5CC, const std::vector<double>& nPlusCC,
			       const std::vector<double>& totalCC)

{
  _energy = energy;
  _halfLife = halfLife;
  _angularMomentum = angularMomentum;
  unsigned int i;
  for (i=0; i<eGamma.size(); i++)
    {
      _energies.push_back(eGamma[i]);
      _weights.push_back(wGamma[i]);
      _polarities.push_back(polarities[i]);
      _kCC.push_back( kCC[i]);
      _l1CC.push_back( l1CC[i]);
      _l2CC.push_back( l2CC[i]);
      _l3CC.push_back( l3CC[i]);
      _m1CC.push_back( m1CC[i]);
      _m2CC.push_back( m2CC[i]);
      _m3CC.push_back( m3CC[i]);
      _m4CC.push_back( m4CC[i]);
      _m5CC.push_back( m5CC[i]);
      _nPlusCC.push_back( nPlusCC[i]);
      _totalCC.push_back( totalCC[i]);
    }
  _nGammas = _energies.size();
  MakeProbabilities();
  MakeCumProb();
 // G4cout << "####### Incrementing "<<Increment(1)<<G4endl;
}

G4NuclearLevel::~G4NuclearLevel()
{ 
 // G4cout << "####### Decrementing "<<Increment(-1)<<G4endl;
}


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


const std::vector<double>& G4NuclearLevel::GammaEnergies() const
{
  return _energies;
}
 
const std::vector<double>& G4NuclearLevel::GammaWeights() const
{
  return _weights;
}
 

const std::vector<double>& G4NuclearLevel::GammaProbabilities() const
{
  return _prob;
}
 

const std::vector<double>& G4NuclearLevel::GammaCumulativeProbabilities() const
{
  return _cumProb;
}
 

const std::vector<double>& G4NuclearLevel::GammaPolarities() const
{
  return _polarities;
}
 
const std::vector<double>& G4NuclearLevel::KConvertionProbabilities() const
{
  return _kCC;
}
 
const std::vector<double>& G4NuclearLevel::L1ConvertionProbabilities() const
{
  return _l1CC;
}
 
const std::vector<double>& G4NuclearLevel::L2ConvertionProbabilities() const
{
  return _l2CC;
}
 
const std::vector<double>& G4NuclearLevel::L3ConvertionProbabilities() const
{
  return _l3CC;
}
 
const std::vector<double>& G4NuclearLevel::M1ConvertionProbabilities() const
{
  return _m1CC;
}
 
const std::vector<double>& G4NuclearLevel::M2ConvertionProbabilities() const
{
  return _m2CC;
}
 
const std::vector<double>& G4NuclearLevel::M3ConvertionProbabilities() const
{
  return _m3CC;
}
 
const std::vector<double>& G4NuclearLevel::M4ConvertionProbabilities() const
{
  return _m4CC;
}
 
const std::vector<double>& G4NuclearLevel::M5ConvertionProbabilities() const
{
  return _m5CC;
}
 
const std::vector<double>& G4NuclearLevel::NPlusConvertionProbabilities() const
{
  return _nPlusCC;
}
 
const std::vector<double>& G4NuclearLevel::TotalConvertionProbabilities() const
{
  return _totalCC;
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
      sum += _weights[i]*(1+_totalCC[i]);
    }

  for (i=0; i<_nGammas; i++)
    {
      if (sum > 0.) { _prob.push_back(_weights[i]*(1+_totalCC[i])/ sum); }
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



