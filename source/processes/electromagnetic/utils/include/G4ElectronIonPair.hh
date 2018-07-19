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
// $Id: G4ElectronIonPair.hh 95657 2016-02-17 13:03:36Z gcosmo $
//
//
#ifndef G4ElectronIonPair_h
#define G4ElectronIonPair_h 1

// -------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4ElectronIonPair
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 08.07.2008
//
// Modifications:
//
//
// Class Description:
//   Compution on number of electon-ion or electorn-hole pairs
//   at the step of a particle and sampling ionisation points 
//   in space
//
// Based on ICRU Report 31, 1979
// "Average Energy Required to Produce an Ion Pair" 
//
// 06.04.2010 V. Grichine, substitute Gauss by Gamma for ionisation 
//                         distribution at fixed energy deposition
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4VProcess.hh"
#include "Randomize.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Material;

class G4ElectronIonPair
{
public: 

  explicit G4ElectronIonPair(G4int verb);

  virtual ~G4ElectronIonPair();

  // compute mean number of ionisation points at a step
  G4double MeanNumberOfIonsAlongStep(const G4ParticleDefinition*, 
				     const G4Material*,
				     G4double edepTotal,
				     G4double edepNIEL = 0.0);

  inline G4double MeanNumberOfIonsAlongStep(const G4Step*); 

  inline G4int SampleNumberOfIonsAlongStep(const G4Step*); 

  // returns pointer to the new vector of positions of
  // ionisation points in the World coordinate system 
  std::vector<G4ThreeVector>* SampleIonsAlongStep(const G4Step*);

  // compute number of holes in the atom after PostStep interaction
  G4int ResidualeChargePostStep(const G4ParticleDefinition*,
				const G4TrackVector* secondary = nullptr,
				G4int processSubType = -1) const;

  inline G4int ResidualeChargePostStep(const G4Step*) const;

  // find mean energies per ionisation 
  G4double FindG4MeanEnergyPerIonPair(const G4Material*) const;

  // dump mean energies per ionisation used in run time
  void DumpMeanEnergyPerIonPair() const;

  // dump G4 list
  void DumpG4MeanEnergyPerIonPair() const;

  inline void SetVerbose(G4int);

private:

  void Initialise();

  G4double FindMeanEnergyPerIonPair(const G4Material*) const;

  // hide assignment operator
  G4ElectronIonPair & operator=(const G4ElectronIonPair &right) = delete;
  G4ElectronIonPair(const G4ElectronIonPair&) = delete;

  // cache
  const G4Material*  curMaterial;
  G4double           curMeanEnergy;

  G4double invFanoFactor;
  
  G4int    verbose;             
  G4int    nMaterials;

  // list of G4 NIST materials with mean energy per ion defined 
  std::vector<G4double>  g4MatData;
  std::vector<G4String>  g4MatNames;
};

inline G4double 
G4ElectronIonPair::MeanNumberOfIonsAlongStep(const G4Step* step) 
{
  return MeanNumberOfIonsAlongStep(step->GetTrack()->GetParticleDefinition(),
				   step->GetPreStepPoint()->GetMaterial(),
				   step->GetTotalEnergyDeposit(),
				   step->GetNonIonizingEnergyDeposit());
}

inline G4int 
G4ElectronIonPair::SampleNumberOfIonsAlongStep(const G4Step* step) 
{
  // use gamma distribution with mean value n=meanion and 
  // dispersion D=meanion/invFanoFactor
  G4double meanion = MeanNumberOfIonsAlongStep(step);
  return G4lrint(G4RandGamma::shoot(meanion*invFanoFactor,invFanoFactor));
} 

inline G4int 
G4ElectronIonPair::ResidualeChargePostStep(const G4Step* step) const
{
  G4int subtype = -1;
  const G4VProcess* proc = step->GetPostStepPoint()->GetProcessDefinedStep();
  if(proc) { subtype = proc->GetProcessSubType(); }
  return ResidualeChargePostStep(step->GetTrack()->GetParticleDefinition(),
				 step->GetSecondary(),
				 subtype);
}

inline void G4ElectronIonPair::SetVerbose(G4int val)
{
  verbose = val;
}

#endif

