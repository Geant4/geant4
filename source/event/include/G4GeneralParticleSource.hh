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

private:

  void IntensityNormalization();

private:
  G4bool normalised;
  G4int currentSourceIdx;
  G4SingleParticleSource* currentSource;
  std::vector <G4SingleParticleSource*> sourceVector;
  std::vector <G4double> sourceIntensity;
  std::vector <G4double>sourceProbability;

  G4GeneralParticleSourceMessenger* theMessenger;
  
};

#endif
