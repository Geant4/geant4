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
// $Id: G4NuclearLevel.hh,v 1.7 2010-11-17 19:17:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//	  06 Oct 2010, M. Kelsey (kelsey@slac.stanford.edu)
//		Add friendship for G4NuclearLevelManager; define private
//		constructors without vectors.
//
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              Added K->N+ internal  conversion coefficiencies and their access
//              functions      
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data. 
//
//        28 October 2010, V.Ivanchenko moved copy constructor to source, cleanup
//
//      
// -------------------------------------------------------------------

#ifndef G4NUCLEARLEVEL_HH
#define G4NUCLEARLEVEL_HH 1

#include "globals.hh"
#include <vector>

class G4NuclearLevel 
{

public:
  G4NuclearLevel(G4double energy, G4double halfLife,
		 G4double angularMomentum, const std::vector<G4double>& eGamma,
		 const std::vector<G4double>& wGamma, const std::vector<G4double>& polarities,
		 const std::vector<G4double>& kCC, const std::vector<G4double>& l1CC,
		 const std::vector<G4double>& l2CC, const std::vector<G4double>& l3CC,
		 const std::vector<G4double>& m1CC, const std::vector<G4double>& m2CC,
		 const std::vector<G4double>& m3CC, const std::vector<G4double>& m4CC,
		 const std::vector<G4double>& m5CC, const std::vector<G4double>& nPlusCC,
		 const std::vector<G4double>& totalCC);

  ~G4NuclearLevel();

  const std::vector<G4double>& GammaEnergies() const;
 
  const std::vector<G4double>& GammaWeights() const;

  const std::vector<G4double>& GammaProbabilities() const;

  const std::vector<G4double>& GammaCumulativeProbabilities() const;

  const std::vector<G4double>& GammaPolarities() const;

  const std::vector<G4double>& KConvertionProbabilities() const;

  const std::vector<G4double>& L1ConvertionProbabilities() const;

  const std::vector<G4double>& L2ConvertionProbabilities() const;

  const std::vector<G4double>& L3ConvertionProbabilities() const;

  const std::vector<G4double>& M1ConvertionProbabilities() const;

  const std::vector<G4double>& M2ConvertionProbabilities() const;

  const std::vector<G4double>& M3ConvertionProbabilities() const;

  const std::vector<G4double>& M4ConvertionProbabilities() const;

  const std::vector<G4double>& M5ConvertionProbabilities() const;

  const std::vector<G4double>& NPlusConvertionProbabilities() const;

  const std::vector<G4double>& TotalConvertionProbabilities() const;

  G4double Energy() const;

  G4double AngularMomentum() const;

  G4double HalfLife() const;

  G4int NumberOfGammas() const;

  void PrintAll() const;  

  G4bool operator==(const G4NuclearLevel &right) const;
  G4bool operator!=(const G4NuclearLevel &right) const;
  G4bool operator<(const G4NuclearLevel &right) const;
  const G4NuclearLevel& operator=(const G4NuclearLevel &right);

  G4NuclearLevel(const G4NuclearLevel &right);

private:  

  friend class G4NuclearLevelManager;

  G4NuclearLevel();

  G4NuclearLevel(G4double energy, G4double halfLife,
		 G4double angularMomentum);

  void Finalize();

  void MakeProbabilities();
  void MakeCumProb();
  
  G4int Increment(G4int aF);
  
  std::vector<G4double> _energies;
  std::vector<G4double> _weights;
  std::vector<G4double> _prob;
  std::vector<G4double> _cumProb;
  std::vector<G4double> _polarities;
  std::vector<G4double> _kCC;
  std::vector<G4double> _l1CC;  
  std::vector<G4double> _l2CC;  
  std::vector<G4double> _l3CC;  
  std::vector<G4double> _m1CC;  
  std::vector<G4double> _m2CC;  
  std::vector<G4double> _m3CC;  
  std::vector<G4double> _m4CC;  
  std::vector<G4double> _m5CC;  
  std::vector<G4double> _nPlusCC;  
  std::vector<G4double> _totalCC;  

  G4double _energy;
  G4double _halfLife;
  G4double _angularMomentum;
  G4int _nGammas;
};

#endif






