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
////////////////////////////////////////////////////////////////////////////////
//
//  GEANT4 Class header file
//
//  G4VRadioactiveDecay
//
//  Authors: F. Lei and P.R. Truscott
//  Date:   15 June 2000
//
//  9 August 2017 D.H. Wright (SLAC) introduced multi-threaded mode
//  30 May 2024   V. Ivanchenko rename the class 
//
//  Description: This class is the base for simulation of radioactive
//               decays. It performs alpha, beta, gamma, electron capture,
//               and isomeric transition decays of radioactive nuclei.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef G4VRadioactiveDecay_h
#define G4VRadioactiveDecay_h 1

#include <vector>
#include <map>
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4ios.hh"
#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4ParticleChangeForRadDecay.hh"

#include "G4NucleusLimits.hh"
#include "G4ThreeVector.hh"
#include "G4RadioactiveDecayMode.hh"

class G4Fragment;
class G4RadioactiveDecayMessenger;
class G4PhotonEvaporation;
class G4Ions;
class G4DecayTable;
class G4ITDecay;

typedef std::map<G4String, G4DecayTable*> DecayTableMap;

class G4VRadioactiveDecay : public G4VRestDiscreteProcess 
{
  // class description

  // Implementation of the radioactive decay process which simulates the
  // decays of radioactive nuclei.  These nuclei are submitted to RDM as
  // G4Ions.  The required half-lives and decay schemes are retrieved from
  // the Radioactivity database which was derived from ENSDF.
  // All decay products are submitted back to the particle tracking process
  // through the G4ParticleChangeForRadDecay object.
  // class description - end 

  public: // with description

    G4VRadioactiveDecay(const G4String& processName="RadioactiveDecay",
                        const G4double timeThreshold=-1.0);
    ~G4VRadioactiveDecay() override;

    G4bool IsApplicable(const G4ParticleDefinition&) override;
    // Return true if the specified isotope is
    //  1) defined as "nucleus" and
    //  2) it is within theNucleusLimit

    G4VParticleChange* AtRestDoIt(const G4Track& theTrack,
                                  const G4Step& theStep) override;

    G4VParticleChange* PostStepDoIt(const G4Track& theTrack,
                                    const G4Step& theStep) override;

    void BuildPhysicsTable(const G4ParticleDefinition &) override;

    void ProcessDescription(std::ostream& outFile) const override;

    virtual G4VParticleChange* DecayIt(const G4Track& theTrack,
                                       const G4Step& theStep);

    // Return decay table if it exists, if not, load it from file
    G4DecayTable* GetDecayTable(const G4ParticleDefinition*);

    // Select a logical volume in which RDM applies
    void SelectAVolume(const G4String& aVolume);

    // Remove a logical volume from the RDM applied list
    void DeselectAVolume(const G4String& aVolume);

    // Select all logical volumes for the application of RDM
    void SelectAllVolumes();

    // Remove all logical volumes from RDM applications
    void DeselectAllVolumes();

    // Enable/disable ARM
    inline void SetARM(G4bool arm) { applyARM = arm; }

    G4DecayTable* LoadDecayTable(const G4Ions*);
    // Load the decay data of isotope theParentNucleus

    void AddUserDecayDataFile(G4int Z, G4int A, const G4String& filename);
    // Allow the user to replace the radio-active decay data provided in Geant4
    // by its own data file for a given isotope

    inline void SetNucleusLimits(G4NucleusLimits theNucleusLimits1)
      { theNucleusLimits = theNucleusLimits1; }
    // Sets theNucleusLimits which specifies the range of isotopes
    // the G4VRadioactiveDecay applies.

    // Returns theNucleusLimits which specifies the range of isotopes used
    // by G4VRadioactiveDecay
    inline G4NucleusLimits GetNucleusLimits() const { return theNucleusLimits; }

    inline void SetDecayDirection(const G4ThreeVector& theDir) {
      forceDecayDirection = theDir.unit();
    }

    inline const G4ThreeVector& GetDecayDirection() const {
      return forceDecayDirection; 
    }

    inline void SetDecayHalfAngle(G4double halfAngle=0.*CLHEP::deg) {
      forceDecayHalfAngle = std::min(std::max(0.*CLHEP::deg,halfAngle),180.*CLHEP::deg);
    }

    inline G4double GetDecayHalfAngle() const { return forceDecayHalfAngle; }

    // Force direction (random within half-angle) for "visible" daughters
    // (applies to electrons, positrons, gammas, neutrons, protons or alphas)
    inline void SetDecayCollimation(const G4ThreeVector& theDir,
                                    G4double halfAngle = 0.*CLHEP::deg) {
      SetDecayDirection(theDir);
      SetDecayHalfAngle(halfAngle);
    }

    // Ignore radioactive decays at rest of nuclides happening
    // after this (very long) time threshold
    inline void SetThresholdForVeryLongDecayTime(const G4double inputThreshold) {
      fThresholdForVeryLongDecayTime = std::max( 0.0, inputThreshold );
    }
    inline G4double GetThresholdForVeryLongDecayTime() const {
      return fThresholdForVeryLongDecayTime;
    }

    void StreamInfo(std::ostream& os, const G4String& endline);

    G4VRadioactiveDecay(const G4VRadioactiveDecay& right) = delete;
    G4VRadioactiveDecay& operator=(const G4VRadioactiveDecay& right) = delete;

  protected:

    G4double GetMeanFreePath(const G4Track& theTrack, G4double previousStepSize,
                             G4ForceCondition* condition) override;

    G4double GetMeanLifeTime(const G4Track& theTrack,
                             G4ForceCondition* condition) override;

    // sampling of products 
    void DecayAnalog(const G4Track& theTrack, G4DecayTable*);

    // sampling products at rest
    G4DecayProducts* DoDecay(const G4ParticleDefinition&, G4DecayTable*);

    // Apply directional bias for "visible" daughters (e+-, gamma, n, p, alpha)
    void CollimateDecay(G4DecayProducts* products);
    void CollimateDecayProduct(G4DynamicParticle* product);
    G4ThreeVector ChooseCollimationDirection() const;

    // ParticleChange for decay process
    G4ParticleChangeForRadDecay fParticleChangeForRadDecay;

    G4RadioactiveDecayMessenger* theRadioactiveDecayMessenger;
    G4PhotonEvaporation* photonEvaporation;
    G4ITDecay* decayIT;

    std::vector<G4String> ValidVolumes;
    G4bool isAllVolumesMode{true};

    static const G4double levelTolerance;

    // Library of decay tables
    static DecayTableMap* master_dkmap;

  private:

    G4NucleusLimits theNucleusLimits;

    G4bool isInitialised{false};
    G4bool applyARM{true};

    // Parameters for pre-collimated (biased) decay products
    G4ThreeVector forceDecayDirection{G4ThreeVector(0., 0., 0.)};
    G4double forceDecayHalfAngle{0.0};
    static const G4ThreeVector origin;	// (0,0,0) for convenience

    // Radioactive decay database directory path 
    static G4String dirPath;

    // User define radioactive decay data files replacing some files in the G4RADECAY database
    static std::map<G4int, G4String>* theUserRDataFiles;

    // The last RadDecayMode
    G4RadioactiveDecayMode theRadDecayMode{G4RadioactiveDecayMode::IT};
 
    // Ignore radioactive decays at rest of nuclides happening after this (very long) time threshold
    G4double fThresholdForVeryLongDecayTime;
};

#endif

