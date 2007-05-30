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
///////////////////////////////////////////////////////////////////////////////
//
// MODULE:       MultipleSource.hh
//
// Version:      1.0
// Date:         1/12/07
// Author:       Jerome Maurice Verbeke 
// Organisation: LLNL
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
// 01/12/2007, Jerome Maurice Verbeke
//
///////////////////////////////////////////////////////////////////////////////
//
// Class Description:
//
// The spontaneous Fission Source is similar to the General Particle
// source. The only difference is that currentSource is of class
// SingleSource (that has a virtual method for GeneratePrimaryVertex,
// as opposed to G4SingleParticleSource which has a non-virtual method)
//
// Documentation on G4GeneralParticleSource is avaialable at 
// http://reat.space.qinetiq.com/gps.  These include:
//       User Requirement Document (URD)
//       Software Specification Documents (SSD)
//       Software User Manual (SUM): on-line version available
//       Technical Note (TN) on the physics and algorithms
//
///////////////////////////////////////////////////////////////////////////////
//
// MEMBER FUNCTIONS
// ----------------
//
//  MultipleSource()
//     Constructor:  Initializes variables and instantiates the 
//
//  MultipleSource(SingleSource*, G4double)
//     Constructor:  Initializes variables, creates a 
//                   particle gun with the given strength
//                   and instantiates the generator classes
//
//  ~MultipleSource()
//     Destructor: 
//
//  G4int GetNumberofSource() 
//     Return the number of particle gun defined
//
//  void ListSource()
//     List the particle guns defined
//
//  void SetCurrentSourceto(G4int) 
//     set the current gun to the specified one so its definition can be changed
//
//  void SetCurrentSourceIntensity(G4double)
//     change the current particle gun strength
//  
//  void SetMultipleVertex(G4bool )  
//     Set if multiple vertex per event.

//  SingleSource* GetCurrentSource()
//     return the pointer to current particle gun
// 
//  G4int GetCurrentSourceIndex() 
//     return the index of the current particle gun
//
//  G4double GetCurrentSourceIntensity() 
//     return the strength of the current gun
//
//  void ClearAll()
//     remove all defined aprticle gun
//
//  void AddaSource (G4double)
//     add a new particle gun with the specified strength
//
//  void AddaSource (SingleSource*, G4double)
//     add this particle gun with the specified strength
//
//  void DeleteaSource(G4int);
//     delete the specified particle gun
//
//  void SetParticleDefinition ();
//  G4ParticleDefinition * GetParticleDefinition () 
//     Get/Set the particle definition of the primary track
//
//  void SetParticleCharge(G4double aCharge) 
//     set the charge state of the primary track
//
//  void SetParticlePolarization (G4ThreeVector aVal) 
//  G4ThreeVector GetParticlePolarization ()
//     Set/Get the polarization state of the primary track
//
//  void SetParticleTime(G4double aTime)  { particle_time = aTime; };
//  G4double GetParticleTime()  { return particle_time; };
//     Set/Get the Time.
//
//  void SetNumberOfParticles(G4int i) 
//  G4int GetNumberOfParticles() 
//     set/get the number of particles to be generated in the primary track
//
//  G4ThreeVector GetParticlePosition()  
//  G4ThreeVector GetParticleMomentumDirection()  
//  G4double GetParticleEnergy()  
//     get the position, direction, and energy of the current particle 
///////////////////////////////////////////////////////////////////////////////
//
#ifndef MultipleSource_h
#define MultipleSource_h 1

#include "globals.hh"
#include <vector>

#include "SingleSource.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"

class MultipleSource : public G4VPrimaryGenerator
{
  //
public:

  MultipleSource();
  MultipleSource(SingleSource* src, G4double);
  ~MultipleSource();

  void GeneratePrimaryVertex(G4Event*);

  G4int GetNumberofSource() { return G4int(sourceVector.size()); };
  void ListSource();
  void SetCurrentSourceto(G4int) ;
  void SetCurrentSourceIntensity(G4double);
  SingleSource* GetCurrentSource() {return currentSource;};
  G4int GetCurrentSourceIndex() { return currentSourceIdx; };
  G4double GetCurrentSourceIntensity() { return sourceIntensity[currentSourceIdx]; };
  void ClearAll();
  void AddaSource (G4double);
  void AddaSource (SingleSource*, G4double);
  void DeleteaSource(G4int);

  // Set the verbosity level.
  void SetVerbosity(G4int i) {currentSource->SetVerbosity(i);} ;

  // Set if multiple vertex per event.
  void SetMultipleVertex(G4bool av) {multiple_vertex = av;} ;

  // set if flat_sampling is applied in multiple source case

  void SetFlatSampling(G4bool av) {flat_sampling = av; normalised = false;} ;

  // Set the particle species
  void SetParticleDefinition (G4ParticleDefinition * aParticleDefinition) 
    {currentSource->SetParticleDefinition(aParticleDefinition); } ;

  G4ParticleDefinition * GetParticleDefinition () { return currentSource->GetParticleDefinition();} ;

  void SetParticleCharge(G4double aCharge) { currentSource->SetParticleCharge(aCharge); } ;

  // Set polarization
  void SetParticlePolarization (G4ThreeVector aVal) {currentSource->SetParticlePolarization(aVal);};
  G4ThreeVector GetParticlePolarization ()  {return currentSource->GetParticlePolarization();};

  // Set Time.
  void SetParticleTime(G4double aTime)  { currentSource->SetParticleTime(aTime); };
  G4double GetParticleTime()  { return currentSource->GetParticleTime(); };

  void SetNumberOfParticles(G4int i)  { currentSource->SetNumberOfParticles(i); };
  //
  G4int GetNumberOfParticles() { return currentSource->GetNumberOfParticles(); };
  G4ThreeVector GetParticlePosition()  { return currentSource->GetParticlePosition();};
  G4ThreeVector GetParticleMomentumDirection()  { return currentSource->GetParticleMomentumDirection();};
  G4double GetParticleEnergy()  {return currentSource->GetParticleEnergy();};

private:

  void IntensityNormalization();

private:
  G4bool multiple_vertex;
  G4bool flat_sampling;
  G4double weight_change;
  G4bool normalised;
  G4int currentSourceIdx;
  SingleSource* currentSource;
  std::vector <SingleSource*> sourceVector;
  std::vector <G4double> sourceIntensity;
  std::vector <G4double>sourceProbability;
};

#endif
