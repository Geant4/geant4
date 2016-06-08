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
//      Creation date: 25 October 1998
//
//      Modifications: 
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
#include "G4DataVector.hh"

class G4NuclearLevel 
{

public:

  G4NuclearLevel(const G4double energy, const G4double halfLife,
		 const G4double angularMomentum, const G4DataVector& eGamma,
		 const G4DataVector& wGamma, const G4DataVector& polarities);

  G4NuclearLevel() {};
  ~G4NuclearLevel();

  const G4DataVector& GammaEnergies() const;
 
  const G4DataVector& GammaWeights() const;

  const G4DataVector& GammaProbabilities() const;

  const G4DataVector& GammaCumulativeProbabilities() const;

  const G4DataVector& GammaPolarities() const;

  G4double Energy() const;

  G4double AngularMomentum() const;

  G4double HalfLife() const;

  G4int NumberOfGammas() const;

  void PrintAll() const;  

  G4bool operator==(const G4NuclearLevel &right) const;
  G4bool operator!=(const G4NuclearLevel &right) const;
  G4bool operator<(const G4NuclearLevel &right) const;

protected:

private:  

  //  G4NuclearLevel(const G4NuclearLevel &right);  
  //  const G4NuclearLevel& operator=(const G4NuclearLevel &right);

  void MakeProbabilities();
  void MakeCumProb();

  G4DataVector _energies;
  G4DataVector _weights;
  G4DataVector _prob;
  G4DataVector _cumProb;
  G4DataVector _polarities;
  G4double _energy;
  G4double _halfLife;
  G4double _angularMomentum;
  G4int _nGammas;

};

#endif




















