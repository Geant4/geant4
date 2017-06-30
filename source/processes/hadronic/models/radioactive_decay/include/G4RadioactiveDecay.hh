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
#ifndef G4RadioactiveDecay_h
#define G4RadioactiveDecay_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4RadioactiveDecay.hh
//
// Version:             0.b.4
// Date:                14/04/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
// 17 October 2011, L Desorgher - Add the method AddUserDecayDataFile
//
// 01 June 2011, M. Kelsey -- Add directional biasing interface to allow for
//		"collimation" of decay daughters.
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// 13 April 2000, F Lei, DERA UK
// 0.b.4 release. No change to this file     
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <map>
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4ios.hh"
#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4ParticleChangeForRadDecay.hh"

#include "G4NucleusLimits.hh"
#include "G4RadioactiveDecayRate.hh"
#include "G4RadioactiveDecayRateVector.hh"
#include "G4RadioactivityTable.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

class G4Fragment;
class G4RadioactiveDecaymessenger;
class G4PhotonEvaporation;

typedef std::vector<G4RadioactiveDecayRateVector> G4RadioactiveDecayRateTable;
typedef std::vector<G4RadioactiveDecayRate> G4RadioactiveDecayRates;
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

    // Return true if the specified isotope is
    //  1) defined as "nucleus" and
    //  2) it is within theNucleusLimit
    G4bool IsApplicable(const G4ParticleDefinition&);

    // Return decay table if it exists, if not, load it from file
    G4DecayTable* GetDecayTable(const G4ParticleDefinition*);

    // Select a logical volume in which RDM applies
    void SelectAVolume(const G4String aVolume);

    // Remove a logical volume from the RDM applied list
    void DeselectAVolume(const G4String aVolume);

    // Select all logical volumes for the application of RDM
    void SelectAllVolumes();

    // Remove all logical volumes from RDM applications
    void DeselectAllVolumes();

    // Set the decay biasing scheme using the data in "filename"
    void SetDecayBias(G4String filename);

    // Set the half-life threshold for isomer production
    void SetHLThreshold(G4double hl) {halflifethreshold = hl;}

    // Enable/disable ICM
    void SetICM(G4bool icm) {applyICM = icm;} 

    // Enable/disable ARM
    void SetARM(G4bool arm) {applyARM = arm;}

    // Set source exposure function using histograms in "filename"
    void SetSourceTimeProfile(G4String filename);

    G4bool IsRateTableReady(const G4ParticleDefinition &);
    // Returns true if the coefficient and decay time table for all the
    // descendants of the specified isotope are ready.
    // used in VR decay mode only

    void AddDecayRateTable(const G4ParticleDefinition&);
    // Calculates the coefficient and decay time table for all the descendents
    // of the specified isotope.  Adds the calculated table to the private data
    // member "theDecayRateTableVector".
    // used in VR decay mode only 

    void GetDecayRateTable(const G4ParticleDefinition&);
    // Used to retrieve the coefficient and decay time table for all the
    // descendants of the specified isotope from "theDecayRateTableVector"
    // and place it in "theDecayRateTable".
    // used in VR decay mode only 

    void SetDecayRate(G4int,G4int,G4double, G4int, std::vector<G4double>,
                      std::vector<G4double>);
    // Sets "theDecayRate" with data supplied in the arguements.
    // used in VR decay mode only 

    std::vector<G4RadioactivityTable*> GetTheRadioactivityTables()
       {return theRadioactivityTables;}
    // Return vector of G4Radioactivity map - should be used in VR mode only

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

    inline G4NucleusLimits GetNucleusLimits() const
      {return theNucleusLimits;}
    // Returns theNucleusLimits which specifies the range of isotopes
    // the G4RadioactiveDecay applies

    inline void SetAnalogueMonteCarlo (G4bool r ) { 
      AnalogueMC  = r; 
      if (!AnalogueMC) halflifethreshold = 1e-6*CLHEP::s;
    }
    // Controls whether G4RadioactiveDecay runs in analogue mode or
    // variance reduction mode.

    inline void SetFBeta (G4bool r ) { FBeta  = r; }
    // Controls whether G4RadioactiveDecay uses fast beta simulation mode

    inline G4bool IsAnalogueMonteCarlo () {return AnalogueMC;}
    // Returns true if the simulation is an analogue Monte Carlo, and false if
    // any of the biassing schemes have been selected.

    inline void SetBRBias (G4bool r) {
      BRBias = r;
      SetAnalogueMonteCarlo(0);
     }
     // Sets whether branching ration bias scheme applies.

    inline void SetSplitNuclei (G4int r) {
      NSplit = r;
      SetAnalogueMonteCarlo(0);
    }
    // Sets the N number for the Nuclei spliting bias scheme

    inline G4int GetSplitNuclei () {return NSplit;}
    //  Returns the N number used for the Nuclei spliting bias scheme

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

    inline void SetDecayCollimation(const G4ThreeVector& theDir,
                                    G4double halfAngle = 0.*CLHEP::deg) {
      SetDecayDirection(theDir);
      SetDecayHalfAngle(halfAngle);
    }

    // Force direction (random within half-angle) for "visible" daughters
    // (applies to electrons, positrons, gammas, neutrons, protons or alphas)

    void BuildPhysicsTable(const G4ParticleDefinition &);

    G4VParticleChange* DecayIt(const G4Track& theTrack,
                               const G4Step&  theStep);

  protected:

    G4DecayProducts* DoDecay(const G4ParticleDefinition& theParticleDef);

    // Apply directional bias for "visible" daughters (e+-, gamma, n, p, alpha)
    void CollimateDecay(G4DecayProducts* products);
    void CollimateDecayProduct(G4DynamicParticle* product);
    G4ThreeVector ChooseCollimationDirection() const;

    G4double GetMeanFreePath(const G4Track& theTrack, G4double previousStepSize,
                             G4ForceCondition* condition);

    G4double GetMeanLifeTime(const G4Track& theTrack,
                             G4ForceCondition* condition);

    G4double ConvolveSourceTimeProfile(const G4double, const G4double);

    G4double GetDecayTime();

    G4int GetDecayTimeBin(const G4double aDecayTime);

    //Add gamma,Xray,conversion,and auger electrons for bias mode
    void AddDeexcitationSpectrumForBiasMode(G4ParticleDefinition* apartDef,
                                            G4double weight,
                                            G4double currenTime,
                                            std::vector<double>& weights_v,
                                            std::vector<double>& times_v,
                                            std::vector<G4DynamicParticle*>& secondaries_v);

  private:

    G4RadioactiveDecay(const G4RadioactiveDecay &right);
    G4RadioactiveDecay & operator=(const G4RadioactiveDecay &right);

    G4RadioactiveDecaymessenger* theRadioactiveDecaymessenger;
    G4PhotonEvaporation* photonEvaporation;

    G4NucleusLimits theNucleusLimits;

    G4bool isInitialised;
    G4bool AnalogueMC;
    G4bool BRBias;
    G4bool FBeta;
    G4int NSplit;

    G4double halflifethreshold;
    G4bool applyICM;
    G4bool applyARM;

    // Parameters for pre-collimated (biased) decay products
    G4ThreeVector forceDecayDirection;
    G4double      forceDecayHalfAngle;
    static const G4ThreeVector origin;	// (0,0,0) for convenience

    G4int NSourceBin;
    G4double SBin[100];
    G4double SProfile[100];
    G4int NDecayBin;
    G4double DBin[100];
    G4double DProfile[100];

    std::vector<G4String> ValidVolumes;
    bool isAllVolumesMode;

    G4RadioactiveDecayRate theDecayRate;
    G4RadioactiveDecayRates theDecayRateVector;
    G4RadioactiveDecayRateVector theDecayRateTable;
    G4RadioactiveDecayRateTable theDecayRateTableVector;

    // for the radioactivity tables
    std::vector<G4RadioactivityTable*> theRadioactivityTables;
    G4int decayWindows[100];
    static const G4double levelTolerance;

    //User define radioactive decay data files replacing some files in the G4RADECAY database
    std::map<G4int, G4String> theUserRadioactiveDataFiles;

    // Library of decay tables
    DecayTableMap* dkmap;
#ifdef G4MULTITHREADED
    static DecayTableMap* master_dkmap;
#endif

    // Remainder of life time at rest
    G4double fRemainderLifeTime;
    G4int verboseLevel;


    // ParticleChange for decay process
    G4ParticleChangeForRadDecay fParticleChangeForRadDecay;

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
#endif
};

#endif

