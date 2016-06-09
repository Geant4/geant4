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
#ifndef G4REGIONMODEL
#define G4REGIONMODEL

#include "G4ios.hh"
#include <vector>
#include <cmath>
#include "globals.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
//#include "G4NucleusModel.hh"


typedef std::vector<G4double>::const_iterator my_iterator; 

class G4RegionModel //:public G4VRegionModel
{
public:
  G4RegionModel(const G4int numberOfLayers, const G4int A, const G4int Z);
  ~G4RegionModel();

  //instead of A and Z outer radius of the nucleus?
  //void Init(const G4int numberOfLayers, const G4int A, const G4int Z);
  G4double GetDensity(G4double radius);
  G4double GetPotentialEnergy(G4double r, G4int particle);
  G4double GetMaximumNucleonMomentum(G4double radius, G4int nucleon);
  // G4double NumberOfRegions();   

private:

  G4int massNumber;
  G4int protonNumber;
  std::vector<G4double> radius; //contains the outer radiuses of the shells
  std::vector<G4double> density;
  std::vector<G4double> protonFermiEnergy;
  std::vector<G4double> neutronFermiEnergy;
  std::vector<G4double> protonFermiMomentum;
  std::vector<G4double> neutronFermiMomentum;
  
  std::vector<G4double> protonPotentialEnergy;
  std::vector<G4double> neutronPotentialEnergy;

  static const G4double radius0; 
  static const G4double BE;
  //static const G4double pi = 3.141592;

  G4double GetFermiMomentum(G4double density, G4double mass);
  G4double GetFermiEnergy(G4double density, G4double mass);
};  
#endif








