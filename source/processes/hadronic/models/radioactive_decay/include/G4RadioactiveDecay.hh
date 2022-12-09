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
//                                                                            //
//  File:   G4RadioactiveDecay.hh                                             //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   9 August 2017                                                     //
//  Description: version the G4RadioactiveDecay process by F. Lei and         //
//               P.R. Truscott with biasing and activation calculations       //
//               removed to a derived class.  It performs alpha, beta,        //
//               electron capture and isomeric transition decays of           //
//               radioactive nuclei.                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef G4RadioactiveDecay_h
#define G4RadioactiveDecay_h 1

#include <vector>
#include <map>
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4ios.hh"
#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4ParticleChangeForRadDecay.hh"

#include "G4NucleusLimits.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"
#include "G4RadioactiveDecayMode.hh"

class G4Fragment;
class G4RadioactiveDecayMessenger;
class G4PhotonEvaporation;

typedef std::map<G4String, G4DecayTable*> DecayTableMap;


class G4RadioactiveDecay : public G4VRestDiscreteProcess 
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

    G4RadioactiveDecay(const G4String& processName="RadioactiveDecay");
    ~G4RadioactiveDecay();

    virtual void ProcessDescription(std::ostream& outFile) const;

    G4bool IsApplicable(const G4ParticleDefinition&);
    // Return true if the specified isotope is
    //  1) defined as "nucleus" and
    //  2) it is within theNucleusLimit

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
    void SetARM(G4bool arm) {applyARM = arm;}

    G4DecayTable* LoadDecayTable(const G4ParticleDefinition& theParentNucleus);
    // Load the decay data of isotope theParentNucleus

    void AddUserDecayDataFile(G4int Z, G4int A,G4String filename);
    // Allow the user to replace the radio-active decay data provided in Geant4
    // by its own data file for a given isotope

    inline void  SetVerboseLevel(G4int value) {verboseLevel = value;}
    // Sets the VerboseLevel which controls duggering display

    inline G4int GetVerboseLevel() const {return verboseLevel;}
    // Returns the VerboseLevel which controls level of debugging output

    inline void SetNucleusLimits(G4NucleusLimits theNucleusLimits1)
      {theNucleusLimits = theNucleusLimits1 ;}
    // Sets theNucleusLimits which specifies the range of isotopes
    // the G4RadioactiveDecay applies.

    // Returns theNucleusLimits which specifies the range of isotopes used
    // by G4RadioactiveDecay
    inline G4NucleusLimits GetNucleusLimits() const {return theNucleusLimits;}

    inline void SetDecayDirection(const G4ThreeVector& theDir) {
      forceDecayDirection = theDir.unit();
    }

    inline const G4ThreeVector& GetDecayDirection() const {
      return forceDecayDirection; 
    }

    inline void SetDecayHalfAngle(G4double halfAngle=0.*CLHEP::deg) {
      forceDecayHalfAngle = std::min(std::max(0.*CLHEP::deg,halfAngle),180.*CLHEP::deg);
    }

    inline G4double GetDecayHalfAngle() const {return forceDecayHalfAngle;}

    // Force direction (random within half-angle) for "visible" daughters
    // (applies to electrons, positrons, gammas, neutrons, protons or alphas)
    inline void SetDecayCollimation(const G4ThreeVector& theDir,
                                    G4double halfAngle = 0.*CLHEP::deg) {
      SetDecayDirection(theDir);
      SetDecayHalfAngle(halfAngle);
    }

    // Ignore radioactive decays at rest of nuclides happening after this (very long) time threshold
    inline void SetThresholdForVeryLongDecayTime(const G4double inputThreshold) {
      fThresholdForVeryLongDecayTime = std::max( 0.0, inputThreshold );
    }
    inline G4double GetThresholdForVeryLongDecayTime() const {return fThresholdForVeryLongDecayTime;}

    void BuildPhysicsTable(const G4ParticleDefinition &);

    G4VParticleChange* DecayIt(const G4Track& theTrack,
                               const G4Step&  theStep);

  protected:

    void DecayAnalog(const G4Track& theTrack);

    G4DecayProducts* DoDecay(const G4ParticleDefinition& theParticleDef);

    // Apply directional bias for "visible" daughters (e+-, gamma, n, p, alpha)
    void CollimateDecay(G4DecayProducts* products);
    void CollimateDecayProduct(G4DynamicParticle* product);
    G4ThreeVector ChooseCollimationDirection() const;

    G4double GetMeanFreePath(const G4Track& theTrack, G4double previousStepSize,
                             G4ForceCondition* condition);

    G4double GetMeanLifeTime(const G4Track& theTrack,
                             G4ForceCondition* condition);

    // ParticleChange for decay process
    G4ParticleChangeForRadDecay fParticleChangeForRadDecay;

    G4RadioactiveDecayMessenger* theRadioactiveDecayMessenger;
    G4PhotonEvaporation* photonEvaporation;

    std::vector<G4String> ValidVolumes;
    bool isAllVolumesMode;

    static const G4double levelTolerance;

    // Library of decay tables
    DecayTableMap* dkmap;
#ifdef G4MULTITHREADED
    static DecayTableMap* master_dkmap;
#endif

  private:

    void StreamInfo(std::ostream& os, const G4String& endline);

    G4RadioactiveDecay(const G4RadioactiveDecay& right);
    G4RadioactiveDecay& operator=(const G4RadioactiveDecay& right);

    G4NucleusLimits theNucleusLimits;

    G4bool isInitialised;

    G4bool applyARM;

    // Parameters for pre-collimated (biased) decay products
    G4ThreeVector forceDecayDirection;
    G4double      forceDecayHalfAngle;
    static const G4ThreeVector origin;	// (0,0,0) for convenience

    // Radioactive decay database directory path 
    G4String dirPath;

    //User define radioactive decay data files replacing some files in the G4RADECAY database
    std::map<G4int, G4String> theUserRadioactiveDataFiles;

    //The last RadDecayMode
    G4RadioactiveDecayMode theRadDecayMode;

//    // Library of decay tables
//    DecayTableMap* dkmap;
// #ifdef G4MULTITHREADED
//     static DecayTableMap* master_dkmap;
// #endif

    // Remainder of life time at rest
    G4double fRemainderLifeTime;
    G4int verboseLevel;

    // Ignore radioactive decays at rest of nuclides happening after this (very long) time threshold
    G4double fThresholdForVeryLongDecayTime;

    // inline implementations 
    inline
    G4double AtRestGetPhysicalInteractionLength(const G4Track& track,
                                                G4ForceCondition* condition)
    {
      fRemainderLifeTime =
        G4VRestDiscreteProcess::AtRestGetPhysicalInteractionLength(track, condition);
      return fRemainderLifeTime;
    }

    inline
    G4VParticleChange* AtRestDoIt(const G4Track& theTrack,
                                  const G4Step& theStep)
      {return DecayIt(theTrack, theStep);}

    inline
    G4VParticleChange* PostStepDoIt(const G4Track& theTrack,
                                    const G4Step& theStep)
      {return DecayIt(theTrack, theStep);}

#ifdef G4MULTITHREADED
  public:
    static G4Mutex radioactiveDecayMutex;
  protected:
    G4int& NumberOfInstances();
#endif
};

#endif

