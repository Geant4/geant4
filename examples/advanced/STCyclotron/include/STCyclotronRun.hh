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
//  Author: F. Poignant, floriane.poignant@gmail.com
//
// file STCyclotronRun.hh
#ifndef STCyclotronRun_h
#define STCyclotronRun_h 1

#include "G4Run.hh"
#include "globals.hh"
#include <map>
#include <fstream>

/// Run class
///
/// In RecordEvent() there is collected information event per event 
/// from Hits Collections, and accumulated statistic for the run 

class STCyclotronRun : public G4Run
{
  public:

  STCyclotronRun();
  virtual ~STCyclotronRun();
  
  virtual void Merge(const G4Run*);
  virtual void EndOfRun(G4double);
    
  public:

  //Accumulation functions
  void EnergyDepositionTarget(G4double);
  void EnergyDepositionFoil(G4double);
  void CountParticlesTarget();

  //Setting functions
  //parameters for the geometry
  //---> Target
  void SetTargetVolume(G4double);
  void SetTargetDiameter(G4double);
  void SetTargetThickness(G4double);
  //--->Foil
  void SetFoilVolume(G4double);
  void SetFoilThickness(G4double);

  //parameters for the beam
  void SetIrradiationTime(G4double);
  void SetBeamName(G4String);
  void SetBeamEnergy(G4double);
  void SetBeamCurrent(G4double);
  
  //parameters of the run
  void SetPrimariesPerEvent(G4int);
  void SetTimePerEvent(G4double);
  
  void StoreIsotopeID(G4int, G4String);
  std::map<G4int, G4String> GetIsotopeID();
  void ParticleParent(G4String,G4String);

  //Acumulation functions for maps
  //---->Accumulation of isotopes
  void PrimaryIsotopeCountTarget(G4String, G4double);
  void CountStableIsotopes(G4String);
  void DecayIsotopeCountTarget(G4String, G4String, G4double);

  //---->Accumulation of other particles
  void ParticleCountTarget(G4String);

  
  private:

  //Accumulable variables
  G4double fTotalEnergyDepositTarget;
  G4double fTotalEnergyDepositFoil;
  G4int fParticleTarget;

  //Store Isotopes created inside maps during the run
  std::map<G4String,G4int>    fPrimaryIsotopeCountTarget;  
  std::map<G4String,G4double> fPrimaryIsotopeTimeTarget;
  std::map<G4String,G4int>    fParticleCountTarget;
  std::map<G4String,G4double> fDecayIsotopeTimeTarget;
  std::map<G4String,G4String> fDecayIsotopeCountTarget;
  std::map<G4String,G4String> fParticleParent;
  std::map<G4String,G4int>    fStableIsotopeCountTarget;
  std::map<G4String,G4String> fStableIsotopeMumTarget;

  //Stored and used during the run
  std::map<G4int,G4String>    fIsotopeIDTarget;


  //Parameters that may be modified via messenger classes
 
  //--> geometry
  G4double fTargetThickness;
  G4double fTargetDiameter;
  G4double fFoilThickness;
  G4double fTargetVolume;
  G4double fFoilVolume;
  //---> run
  G4int fPrimariesPerEvent;
  G4double fTimePerEvent;
  //--> beam
  G4String fBeamName;
  G4double fBeamCurrent;
  G4double fBeamEnergy;
  
  //Write output in ASCII
  std::ofstream fOutPut;
  std::ofstream fOutPut1;
  std::ofstream fOutPut2;
  std::ofstream fOutPut3;
  std::ofstream fOutPut4;
  

};

#endif

    
