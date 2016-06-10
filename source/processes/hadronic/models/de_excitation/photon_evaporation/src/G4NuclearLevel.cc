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
// $Id: G4NuclearLevel.cc 86986 2014-11-21 13:00:05Z gcosmo $
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
//	  06 Oct 2010, M. Kelsey (kelsey@slac.stanford.edu)
//		Add friendship for G4NuclearLevelManager; define private
//		constructors without vectors.
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
//        28 October 2010, V.Ivanchenko moved copy constructor to source, cleanup
//
// -------------------------------------------------------------------

#include "G4NuclearLevel.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

G4int G4NuclearLevel::Increment(G4int aF)
{
  static G4ThreadLocal G4int instanceCount = 0;
  instanceCount+=aF;
  return instanceCount;
}

G4NuclearLevel::G4NuclearLevel()
  : _energy(0.), _halfLife(0.), _angularMomentum(0.), _nGammas(0) {
  // G4cout << "####### Incrementing "<<Increment(1)<<G4endl;
}

G4NuclearLevel::G4NuclearLevel(G4double energy, G4double halfLife,
			       G4double angularMomentum)
  : _energy(energy), _halfLife(halfLife), _angularMomentum(angularMomentum),
    _nGammas(0) {
  // G4cout << "####### Incrementing "<<Increment(1)<<G4endl;
}

G4NuclearLevel::G4NuclearLevel(G4double energy, G4double halfLife,
			       G4double angularMomentum,
			       const std::vector<G4double>& eGamma,
			       const std::vector<G4double>& wGamma,
			       const std::vector<G4double>& polarities,
			       const std::vector<G4double>& kCC, const std::vector<G4double>& l1CC,
			       const std::vector<G4double>& l2CC, const std::vector<G4double>& l3CC,
			       const std::vector<G4double>& m1CC, const std::vector<G4double>& m2CC,
			       const std::vector<G4double>& m3CC, const std::vector<G4double>& m4CC,
			       const std::vector<G4double>& m5CC, const std::vector<G4double>& nPlusCC,
			       const std::vector<G4double>& totalCC)

  : _energies(eGamma), _weights(wGamma), _polarities(polarities),
     _kCC(kCC), _l1CC(l1CC), _l2CC(l2CC), _l3CC(l3CC),
    _m1CC(m1CC), _m2CC(m2CC), _m3CC(m3CC), _m4CC(m4CC), _m5CC(m5CC),
    _nPlusCC(nPlusCC), _totalCC(totalCC),
    _energy(energy), _halfLife(halfLife), _angularMomentum(angularMomentum)
{
  Finalize();
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

const std::vector<G4double>& G4NuclearLevel::GammaEnergies() const
{
  return _energies;
}
 
const std::vector<G4double>& G4NuclearLevel::GammaWeights() const
{
  return _weights;
}
 

const std::vector<G4double>& G4NuclearLevel::GammaProbabilities() const
{
  return _prob;
}
 

const std::vector<G4double>& G4NuclearLevel::GammaCumulativeProbabilities() const
{
  return _cumProb;
}
 

const std::vector<G4double>& G4NuclearLevel::GammaPolarities() const
{
  return _polarities;
}
 
const std::vector<G4double>& G4NuclearLevel::KConvertionProbabilities() const
{
  return _kCC;
}
 
const std::vector<G4double>& G4NuclearLevel::L1ConvertionProbabilities() const
{
  return _l1CC;
}
 
const std::vector<G4double>& G4NuclearLevel::L2ConvertionProbabilities() const
{
  return _l2CC;
}
 
const std::vector<G4double>& G4NuclearLevel::L3ConvertionProbabilities() const
{
  return _l3CC;
}
 
const std::vector<G4double>& G4NuclearLevel::M1ConvertionProbabilities() const
{
  return _m1CC;
}
 
const std::vector<G4double>& G4NuclearLevel::M2ConvertionProbabilities() const
{
  return _m2CC;
}
 
const std::vector<G4double>& G4NuclearLevel::M3ConvertionProbabilities() const
{
  return _m3CC;
}
 
const std::vector<G4double>& G4NuclearLevel::M4ConvertionProbabilities() const
{
  return _m4CC;
}
 
const std::vector<G4double>& G4NuclearLevel::M5ConvertionProbabilities() const
{
  return _m5CC;
}
 
const std::vector<G4double>& G4NuclearLevel::NPlusConvertionProbabilities() const
{
  return _nPlusCC;
}
 
const std::vector<G4double>& G4NuclearLevel::TotalConvertionProbabilities() const
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
}

void G4NuclearLevel::PrintLevels() const 
{
  G4cout << "   Eexc(MeV)= " << _energy 
	 << " Time(ns)= " << _halfLife/ns << "  Ntrans= " << _nGammas
	 << G4endl;
}

void G4NuclearLevel::Finalize() {
  _nGammas = _energies.size();
  MakeProbabilities();
  MakeCumProb();
}

void G4NuclearLevel::MakeProbabilities()
{
  G4double sum = 0.;
  G4int i = 0;
  for (i=0; i<_nGammas; i++) {
    sum += _weights[i]*(1.+_totalCC[i]);
  }

  if (sum <= 0.) _prob.resize(_nGammas, 1./_nGammas);	// Fast fill
  else {
    _prob.reserve(_nGammas);
    for (i=0; i<_nGammas; i++) {
      _prob.push_back(_weights[i]*(1.+_totalCC[i])/sum);
    }
  }
}


void G4NuclearLevel::MakeCumProb()
{
  if (_nGammas <= 0) return;

  _cumProb.reserve(_nGammas);

  G4double sum = _prob[0];
  _cumProb.push_back(sum);
  
  for (G4int i=1; i<_nGammas; i++) {
    sum += _prob[i];
    _cumProb.push_back(sum);
  }
}

G4NuclearLevel& G4NuclearLevel::operator=(const G4NuclearLevel &right)
{
  if(this != &right)
    {
      _energies = right._energies;
      _weights =right._weights;
      _prob =right._prob;
      _cumProb =right._cumProb;
      _polarities =right._polarities;
      _kCC = right._kCC;
      _l1CC =right._l1CC;
      _l2CC =right._l2CC;
      _l3CC =right._l3CC;
      _m1CC = right._m1CC;
      _m2CC = right._m2CC;
      _m3CC = right._m3CC;
      _m4CC = right._m4CC;
      _m5CC = right._m5CC;
      _nPlusCC = right._nPlusCC;
      _totalCC = right._totalCC;
      _energy = right._energy;
      _halfLife = right._halfLife;
      _angularMomentum = right._angularMomentum;
      _nGammas = right._nGammas;
    }
  return *this;
}

G4NuclearLevel::G4NuclearLevel(const G4NuclearLevel &right)
{
  _energies = right._energies;
  _weights =right._weights;
  _prob =right._prob;
  _cumProb =right._cumProb;
  _polarities =right._polarities;
  _kCC = right._kCC;
  _l1CC =right._l1CC;
  _l2CC =right._l2CC;
  _l3CC =right._l3CC;
  _m1CC = right._m1CC;
  _m2CC = right._m2CC;
  _m3CC = right._m3CC;
  _m4CC = right._m4CC;
  _m5CC = right._m5CC;
  _nPlusCC = right._nPlusCC;
  _totalCC = right._totalCC;
  _energy = right._energy;
  _halfLife = right._halfLife;
  _angularMomentum = right._angularMomentum;
  _nGammas = right._nGammas;
}



