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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// ParticleSource header
// --------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
// This particle source is a shortened version of G4GeneralParticleSource by
// C Ferguson, F Lei & P Truscott (University of Southampton / DERA), with
// some minor modifications.
//////////////////////////////////////////////////////////////////////////////

#ifndef DMXParticleSource_h
#define DMXParticleSource_h 1

#include "G4VPrimaryGenerator.hh"
#include "G4Navigator.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"

#include "DMXParticleSourceMessenger.hh"


class DMXParticleSource : public G4VPrimaryGenerator {

   public:
     DMXParticleSource (); 
     ~DMXParticleSource ();
     void GeneratePrimaryVertex(G4Event *evt);

   public:

     // position distribution  
     void SetPosDisType(G4String);
     void SetPosDisShape(G4String);
     void SetCentreCoords(G4ThreeVector);
     void SetHalfZ(G4double);
     void SetRadius(G4double);
     void GeneratePointSource();
     void GeneratePointsInVolume();
     G4bool IsSourceConfined();
     void ConfineSourceToVolume(G4String);
  
     // angular distribution
     void SetAngDistType(G4String);
     void SetParticleMomentumDirection(G4ParticleMomentum);
     void GenerateIsotropicFlux();

     // energy distribution 
     void SetEnergyDisType(G4String);
     void SetMonoEnergy(G4double);
     void GenerateMonoEnergetic();

     // verbosity
     void SetVerbosity(G4int);
  
     // particle properties
     void SetParticleDefinition(G4ParticleDefinition * aParticleDefinition);
     inline void SetParticleCharge(G4double aCharge)
        { particle_charge = aCharge; }
  
   private:

     // position distribution
     G4String SourcePosType;
     G4String Shape;
     G4double halfz;
     G4double Radius;
     G4ThreeVector CentreCoords;
     G4bool Confine;
     G4String VolName;
     G4String AngDistType;
     G4double MinTheta, MaxTheta, MinPhi, MaxPhi;
     G4double Theta, Phi;
     G4String EnergyDisType;
     G4double MonoEnergy;

     // particle properties 
     G4int                  NumberOfParticlesToBeGenerated;
     G4ParticleDefinition*  particle_definition;
     G4ParticleMomentum     particle_momentum_direction;
     G4double               particle_energy;
     G4double               particle_charge;
     G4ThreeVector          particle_position;
     G4double               particle_time;
     G4ThreeVector          particle_polarization;

     // Verbose
     G4int verbosityLevel;

   private:
  
     DMXParticleSourceMessenger *theMessenger;
     G4Navigator *gNavigator;

  
};


#endif

