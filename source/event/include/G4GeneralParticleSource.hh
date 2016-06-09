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
///////////////////////////////////////////////////////////////////////////////
//
// MODULE:       G4GeneralParticleSource.hh
//
// Version:      2.0
// Date:         5/02/04
// Author:       Fan Lei 
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
// Documentation avaialable at http://reat.space.qinetiq.com/gps
//   These include:
//       User Requirement Document (URD)
//       Software Specification Documents (SSD)
//       Software User Manual (SUM): on-line version available
//       Technical Note (TN) on the physics and algorithms
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
// 26/10/2004, F Lei
//    Added the Multiple_vertex capability.
//    Removed "inline" from all Set/Get methods.
//
// Version 2.0, 05/02/2004, Fan Lei, Created.
//    based on version 1.1 in Geant4 v6.0
//     - Mutilple particle source definition
//     - Re-structured commands
//     - Split the task into smaller classes
//
//     - old commonds have been retained for backward compatibility, but will 
//       be removed in the future. 
//
///////////////////////////////////////////////////////////////////////////////
//
// Class Description:
//
// The General Particle Source is designed to replace the  G4ParticleGun class. 
// It is designed to allow specification of mutiple particle sources, each with 
// independent definitions of particle type, position, direction (or angular) 
// and energy distributions.  
//
///////////////////////////////////////////////////////////////////////////////
//
// MEMBER FUNCTIONS
// ----------------
//
// G4GeneralParticleSource()
//     Constructor:  Initializes variables and instantiates the 
//                   Messenger and generator classes
//
// ~G4GeneralParticleSourceMessenger()
//     Destructor:  deletes Messenger and others
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

//  G4SingleParticleSource* GetCurrentSource()
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
//
///////////////////////////////////////////////////////////////////////////////
//
#ifndef G4GeneralParticleSource_H
#define G4GeneralParticleSource_H 1

#include "globals.hh"
#include <vector>

#include "G4Event.hh"
#include "G4SingleParticleSource.hh"
//
#include "G4GeneralParticleSourceMessenger.hh"

class G4GeneralParticleSource
{
  //
public:

  G4GeneralParticleSource();
  ~G4GeneralParticleSource();

  void GeneratePrimaryVertex(G4Event*);

  G4int GetNumberofSource() { return G4int(sourceVector.size()); };
  void ListSource();
  void SetCurrentSourceto(G4int) ;
  void SetCurrentSourceIntensity(G4double);
  G4SingleParticleSource* GetCurrentSource() {return currentSource;};
  G4int GetCurrentSourceIndex() { return currentSourceIdx; };
  G4double GetCurrentSourceIntensity() { return sourceIntensity[currentSourceIdx]; };
  void ClearAll();
  void AddaSource (G4double);
  void DeleteaSource(G4int);

  // Set the verbosity level.
  void SetVerbosity(G4int i) {currentSource->SetVerbosity(i);} ;

  // Set if multiple vertex per event.
  void SetMultipleVertex(G4bool av) {multiple_vertex = av;} ;

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
  G4bool normalised;
  G4int currentSourceIdx;
  G4SingleParticleSource* currentSource;
  std::vector <G4SingleParticleSource*> sourceVector;
  std::vector <G4double> sourceIntensity;
  std::vector <G4double>sourceProbability;

  G4GeneralParticleSourceMessenger* theMessenger;
  
};

#endif
