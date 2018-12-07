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
//
/// \file Par02PrimaryParticleInformation.hh
/// \brief Definition of the Par02PrimaryParticleInformation class

#ifndef PAR02_PRIMARY_PARTICLE_INFORMATION_H
#define PAR02_PRIMARY_PARTICLE_INFORMATION_H

#include "G4VUserPrimaryParticleInformation.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

/// Primary particle information
///
/// Describes the information that can be associated with a G4PrimaryParticle
/// class object.
/// @author Anna Zaborowska

class Par02PrimaryParticleInformation : public G4VUserPrimaryParticleInformation {
  public:

    /// A constructor.
    /// @param aID A unique particle ID within event.
    /// @param aPDG A PDG code of the particle.
    /// @param aMomentum An initial particle momentum (at the primary vertex).
    Par02PrimaryParticleInformation( G4int aID, G4int aPDG, G4ThreeVector aMomentum );

    virtual ~Par02PrimaryParticleInformation();
    
    /// Prints the information about the particle.
    virtual void Print() const;
    
    /// Sets the initial particle momentum (from particle generator).
    /// @param aMomentum The particle momentum.
    inline void SetMCMomentum( G4ThreeVector aMomentum ) { fMomentumMC = aMomentum; };
    
    /// Gets the initial particle momentum (from particle generator).
    inline G4ThreeVector GetMCMomentum() { return fMomentumMC; };
    
    /// Sets the particle momentum at the entrance to the tracker detector.
    /// @param aMomentum The particle momentum.
    inline void SetTrackerMomentum( G4ThreeVector aMomentum ) 
      { fMomentumTracker = aMomentum; };
    
    /// Gets the particle momentum at the entrance to the tracker detector.
    inline G4ThreeVector GetTrackerMomentum() { return fMomentumTracker; }
    
    /// Sets the tracker detector resolution. 
    /// Currently equal to -1 if AtlFast type of smearing is used.
    /// @param aResolution The detector resolution 
    ///                    (particle type and momentum dependent).
    inline void SetTrackerResolution( G4double aResolution ) 
      { fResolutionTracker = aResolution; };
    
    /// Gets the tracking detector resolution. 
    /// Currently equal to -1 if AtlFast type of smearing is used.
    inline G4double GetTrackerResolution() { return fResolutionTracker; };
    
    /// Sets the tracking detector efficiency. 
    /// Currently not used (efficiency is 1).
    /// @param aEfficiency The detector efficiency.
    inline void SetTrackerEfficiency( G4double aEfficiency ) 
      { fEfficiencyTracker = aEfficiency; };
    
    /// Gets the tracker detector efficiency. 
    /// Currently not used (efficiency is 1).
    inline G4double GetTrackerEfficiency() { return fEfficiencyTracker; };
    
    /// Sets the position of the energy deposit in the electromagnetic calorimeter.
    /// @param aPosition The position of the energy deposit.
    inline void SetEMCalPosition( G4ThreeVector aPosition ) 
      { fPositionEMCal = aPosition; };
    
    /// Gets the position of the energy deposit in the electromagnetic calorimeter.
    inline G4ThreeVector GetEMCalPosition() { return fPositionEMCal; };
    
    /// Sets the energy deposit in the electromagnetic calorimeter.
    /// @param aEnergy The energy deposited in the detector.
    inline void SetEMCalEnergy( G4double aEnergy ) { fEnergyEMCal = aEnergy; };
    
    /// Sets the energy deposit in the electromagnetic calorimeter.
    inline G4double GetEMCalEnergy() { return fEnergyEMCal; };
    
    /// Sets the electromagnetic calorimeter resolution. 
    /// Currently equal to -1 if AtlFast type of smearing is used.
    /// @param aResolution The calorimeter resolution 
    ///                    (particle type and momentum dependent).
    inline void SetEMCalResolution( G4double aResolution ) 
      { fResolutionEMCal = aResolution; };
    
    /// Gets the electromagnetic calorimeter resolution. 
    /// Currently equal to -1 if AtlFast type of smearing is used.
    inline G4double GetEMCalResolution() { return fResolutionEMCal; };
    
    /// Sets the electromagnetic calorimeter efficiency. 
    /// Currently not used (efficiency is 1).
    /// @param aEfficiency The detector efficiency.
    inline void SetEMCalEfficiency( G4double aEfficiency ) 
      { fEfficiencyEMCal = aEfficiency; };
    
    /// Gets the electromagnetic calorimeter efficiency. 
    /// Currently not used (efficiency is 1).
    inline G4double GetEMCalEfficiency() { return fEfficiencyEMCal; };
    
    /// Sets the position of the energy deposit in the hadronic calorimeter.
    /// @param aPosition The position of the energy deposit.
    inline void SetHCalPosition( G4ThreeVector aPosition ) 
      { fPositionHCal = aPosition; };
   
    /// Gets the position of the energy deposit in the hadronic calorimeter.
    inline G4ThreeVector GetHCalPosition() { return fPositionHCal; };
   
    /// Sets the energy deposit in the hadronic calorimeter.
    /// @param aEnergy The energy deposited in the detector.
    inline void SetHCalEnergy( G4double aEnergy ) { fEnergyHCal = aEnergy; };
    
    /// Sets the energy deposit in the hadronic calorimeter.
    inline G4double GetHCalEnergy() { return fEnergyHCal; };
    
    /// Sets the hadronic calorimeter resolution. 
    /// Currently equal to -1 if AtlFast type of smearing is used.
    /// @param aResolution The calorimeter resolution 
    ///                    (particle type and momentum dependent).
    inline void SetHCalResolution( G4double aResolution ) 
      { fResolutionHCal = aResolution; };
    
    /// Gets the hadronic calorimeter resolution. 
    /// Currently equal to -1 if AtlFast type of smearing is used.
    inline G4double GetHCalResolution() { return fResolutionHCal; };
   
    /// Sets the hadronic calorimeter efficiency. 
    /// Currently not used (efficiency is 1).
    /// @param aEfficiency The detector efficiency.
    inline void SetHCalEfficiency( G4double aEfficiency ) 
      { fEfficiencyHCal = aEfficiency; };
    
    /// Gets the hadronic calorimeter efficiency. 
    /// Currently not used (efficiency is 1).
    inline G4double GetHCalEfficiency() { return fEfficiencyHCal; };
    
    /// Gets the particle unique ID (within event). Can be set only in the constructor.
    inline G4int GetPartID() const { return fPartID; };
    
    /// Gets the standard PDG code. Can be set only in the constructor.
    inline G4int GetPDG() const { return fPDG; };

  private:
    
    /// A particle unique ID.
    G4int fPartID;
    
    /// A particle type (PDG code).
    G4int fPDG;
    
    /// A particle initial momentum (from particle generator).
    G4ThreeVector fMomentumMC;
    
    /// A particle momentum at the entrance to the tracking detector.
    G4ThreeVector fMomentumTracker;
    
    /// A resolution of the tracking detector.
    G4double fResolutionTracker;
    
    /// An efficiency of the tracking detector.
    /// Currently not used (equal to 1).
    G4double fEfficiencyTracker;
    
    /// A position of the energy deposited in the electromagnetic calorimeter.
    G4ThreeVector fPositionEMCal;
    
    /// An energy deposited in the electromagnetic calorimeter.
    G4double fEnergyEMCal;
    
    /// The resolution of the electromagnetic calorimeter.
    G4double fResolutionEMCal;
    
    /// The efficiency of the electromagnetic calorimeter. 
    /// Currently not used (equal to 1).
    G4double fEfficiencyEMCal;
    
    /// A position of the energy deposited in the hadronic calorimeter.
    G4ThreeVector fPositionHCal;
    
    /// An energy deposited in the hadronic calorimeter.
    G4double fEnergyHCal;
    
    /// The resolution of the hadronic calorimeter.
    G4double fResolutionHCal;
    
    /// The efficiency of the hadronic calorimeter. 
    /// Currently not used (equal to 1).
    G4double fEfficiencyHCal;
};

#endif

