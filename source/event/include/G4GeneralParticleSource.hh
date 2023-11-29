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
// G4GeneralParticleSource
//
// Class Description:
//
// The General Particle Source is designed to replace the  G4ParticleGun class. 
// It is designed to allow specification of mutiple particle sources, each with 
// independent definitions of particle type, position, direction (or angular) 
// and energy distributions.  

// Author: Fan Lei, QinetiQ ltd
// Customer: ESA/ESTEC
// Version: 2.0
// History:
// - 05/02/2004, F.Lei - Version 2.0. Created.
// - 26/10/2004, F.Lei - Added the Multiple_vertex capability.
// - 25/03/2014, A.Green - Various changes to use the new class
//               G4GeneralParticleSourceData, mostly just transparent wrappers
//               around the thread safe object.
// --------------------------------------------------------------------
#ifndef G4GeneralParticleSource_hh
#define G4GeneralParticleSource_hh 1

#include "globals.hh"
#include <vector>

#include "G4Event.hh"
#include "G4SingleParticleSource.hh"
#include "G4GeneralParticleSourceMessenger.hh"
#include "G4GeneralParticleSourceData.hh"

class G4SingleParticleSource;

class G4GeneralParticleSource : public G4VPrimaryGenerator
{
  public:

    G4GeneralParticleSource();
      // Initialize variables and instantiates the messenger and
      // generator classes

   ~G4GeneralParticleSource() override;
      // Delete messenger and others

    void GeneratePrimaryVertex(G4Event*) override;

    inline G4int GetNumberofSource() { return GPSData->GetSourceVectorSize(); }
      // Return the number of particle gun defined

    void ListSource();
      // List the particle guns defined

    void SetCurrentSourceto(G4int) ;
      // Set the current gun to the specified one so its definition
      // can be changed

    void SetCurrentSourceIntensity(G4double);
      // Change the current particle gun strength

    inline G4SingleParticleSource* GetCurrentSource() const
      { return GPSData->GetCurrentSource(); }
      // Return the pointer to current particle gun
    inline G4int GetCurrentSourceIndex() const
      { return GPSData->GetCurrentSourceIdx(); }
      // Return the index of the current particle gun
    inline G4double GetCurrentSourceIntensity() const
      { return GPSData->GetIntensity(GetCurrentSourceIndex()); }
      // Return the strength of the current gun

    void ClearAll();
      // Remove all defined particle gun
    void AddaSource (G4double);
      // Add a new particle gun with the specified strength
    void DeleteaSource(G4int);
      // Delete the specified particle gun

    inline void SetVerbosity(G4int i) { GPSData->SetVerbosityAllSources(i); }
      // Set the verbosity level.

    inline void SetMultipleVertex(G4bool av) { GPSData->SetMultipleVertex(av); }
      // Set if multiple vertex per event.

    inline void SetFlatSampling(G4bool av)
      { GPSData->SetFlatSampling(av); normalised = false;}
      // Set if flat_sampling is applied in multiple source case

    inline void SetParticleDefinition (G4ParticleDefinition * aPDef) 
      { GPSData->GetCurrentSource()->SetParticleDefinition(aPDef); }
    inline G4ParticleDefinition* GetParticleDefinition () const
      { return GPSData->GetCurrentSource()->GetParticleDefinition(); }
      // Set/Get the particle definition of the primary track

    inline void SetParticleCharge(G4double aCharge)
      { GPSData->GetCurrentSource()->SetParticleCharge(aCharge); }
      // Set the charge state of the primary track

    inline void SetParticlePolarization (G4ThreeVector aVal)
      { GPSData->GetCurrentSource()->SetParticlePolarization(aVal); }
    inline G4ThreeVector GetParticlePolarization () const
      { return GPSData->GetCurrentSource()->GetParticlePolarization(); }
      // Set/Get polarization state of the primary track

    inline void SetParticleTime(G4double aTime)
      { GPSData->GetCurrentSource()->SetParticleTime(aTime); }
    inline G4double GetParticleTime() const
      { return GPSData->GetCurrentSource()->GetParticleTime(); }
      // Set/Get the Time.

    inline void SetNumberOfParticles(G4int i)
      { GPSData->GetCurrentSource()->SetNumberOfParticles(i); }
    inline G4int GetNumberOfParticles() const
      { return GPSData->GetCurrentSource()->GetNumberOfParticles(); }
      // Set/Get the number of particles to be generated in the primary track

    inline G4ThreeVector GetParticlePosition() const
      { return GPSData->GetCurrentSource()->GetParticlePosition(); }
    inline G4ThreeVector GetParticleMomentumDirection() const
      { return GPSData->GetCurrentSource()->GetParticleMomentumDirection(); }
    inline G4double GetParticleEnergy() const
      { return GPSData->GetCurrentSource()->GetParticleEnergy(); }
      // Get the position, direction, and energy of the current particle 

  private:

    void IntensityNormalization();

  private:

    G4bool normalised = false;
      // Helper Boolean, used to reduce number of locks
      // at run time (see GeneratePrimaryVertex)

    G4GeneralParticleSourceMessenger* theMessenger = nullptr;
      // Note this is a shared resource among MT workers
    G4GeneralParticleSourceData* GPSData = nullptr;
      // Note this is a shared resource among MT workers

/** Andrea Dotti Feb 2015
 * GPS messenger design requires some explanation for what distributions
 * parameters are concerned : Each thread has its own GPS
 * since primary generation is a user action.
 * However to save memory the underlying structures that provide the
 * GPS functionalities ( the G4SPS*Distribution classes and the
 * G4SPSRandomGenerator class)
 * are shared among threads. This implies that modifying parameters of sources
 * requires some attention:
 * 1- Only one thread should change source parameters.
 * 2- Changing of parameters can happen only between runs, when is guaranteed
 *    that no thread is accessing them
 * 2- UI commands require that even if messenger is instantiated in a thread
 *    the commands are executed in the master (this is possible since V10.1)
 * The simplest solution is to use UI commands to change GPS parameters and
 * avoid C++ APIs. If this is inevitable a simple solution is to instantiate
 * an instance of G4GeneralParticleSource explicitly in the master thread
 * (for example in G4VUserActionInitialization::BuildForMaster() and set the
 * defaults parameter there).
 */

};

#endif
