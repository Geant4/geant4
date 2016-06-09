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
//      Creation date: 25 October 1998
//
//      Modifications: 
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              Added K->N+ internal  conversion coefficiencies and their access
//              functions      
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data. 
//      
// -------------------------------------------------------------------

#ifndef G4NUCLEARLEVEL_HH
#define G4NUCLEARLEVEL_HH

#include "globals.hh"
#include "G4NuclearLevel.hh"
#include <vector>

class G4NuclearLevel 
{

public:

  G4NuclearLevel(const G4double energy, const G4double halfLife,
		 const G4double angularMomentum, const std::vector<double>& eGamma,
		 const std::vector<double>& wGamma, const std::vector<double>& polarities,
		 const std::vector<double>& kCC, const std::vector<double>& l1CC,
		 const std::vector<double>& l2CC, const std::vector<double>& l3CC,
		 const std::vector<double>& m1CC, const std::vector<double>& m2CC,
		 const std::vector<double>& m3CC, const std::vector<double>& m4CC,
		 const std::vector<double>& m5CC, const std::vector<double>& nPlusCC,
		 const std::vector<double>& totalCC);

  ~G4NuclearLevel();

  const std::vector<double>& GammaEnergies() const;
 
  const std::vector<double>& GammaWeights() const;

  const std::vector<double>& GammaProbabilities() const;

  const std::vector<double>& GammaCumulativeProbabilities() const;

  const std::vector<double>& GammaPolarities() const;

  const std::vector<double>& KConvertionProbabilities() const;

  const std::vector<double>& L1ConvertionProbabilities() const;

  const std::vector<double>& L2ConvertionProbabilities() const;

  const std::vector<double>& L3ConvertionProbabilities() const;

  const std::vector<double>& M1ConvertionProbabilities() const;

  const std::vector<double>& M2ConvertionProbabilities() const;

  const std::vector<double>& M3ConvertionProbabilities() const;

  const std::vector<double>& M4ConvertionProbabilities() const;

  const std::vector<double>& M5ConvertionProbabilities() const;

  const std::vector<double>& NPlusConvertionProbabilities() const;

  const std::vector<double>& TotalConvertionProbabilities() const;

  G4double Energy() const;

  G4double AngularMomentum() const;

  G4double HalfLife() const;

  G4int NumberOfGammas() const;

  void PrintAll() const;  

  G4bool operator==(const G4NuclearLevel &right) const;
  G4bool operator!=(const G4NuclearLevel &right) const;
  G4bool operator<(const G4NuclearLevel &right) const;

    const G4NuclearLevel& operator=(const G4NuclearLevel &right)
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

    G4NuclearLevel(const G4NuclearLevel &right)
    {
      if(this != &right) *this = right;
//      G4cout << "####### Incrementing "<<Increment(1)<<G4endl;
    }

protected:

private:  

  G4NuclearLevel() {G4cout << "Calling default constructor"<<G4endl;};
      
  void MakeProbabilities();
  void MakeCumProb();
  
  G4int Increment(G4int aF);
  
  std::vector<double> _energies;
  std::vector<double> _weights;
  std::vector<double> _prob;
  std::vector<double> _cumProb;
  std::vector<double> _polarities;
  std::vector<double> _kCC;
  std::vector<double> _l1CC;  
  std::vector<double> _l2CC;  
  std::vector<double> _l3CC;  
  std::vector<double> _m1CC;  
  std::vector<double> _m2CC;  
  std::vector<double> _m3CC;  
  std::vector<double> _m4CC;  
  std::vector<double> _m5CC;  
  std::vector<double> _nPlusCC;  
  std::vector<double> _totalCC;  

  G4double _energy;
  G4double _halfLife;
  G4double _angularMomentum;
  G4int _nGammas;
  

};

#endif




















