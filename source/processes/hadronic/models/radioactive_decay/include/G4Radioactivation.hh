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
#ifndef G4Radioactivation_h
#define G4Radioactivation_h 1

#include "G4RadioactiveDecayBase.hh"

#include <vector>
#include <map>
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4ios.hh"
#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4ParticleChangeForRadDecay.hh"

#include "G4NucleusLimits.hh"
#include "G4RadioactiveDecayRatesToDaughter.hh"
#include "G4RadioactiveDecayChainsFromParent.hh"
#include "G4RadioactivityTable.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

class G4Fragment;
class G4RadioactivationMessenger;

typedef std::vector<G4RadioactiveDecayChainsFromParent> G4RadioactiveDecayParentChainTable;
typedef std::vector<G4RadioactiveDecayRatesToDaughter> G4RadioactiveDecayRates;
typedef std::map<G4String, G4DecayTable*> DecayTableMap;


class G4Radioactivation : public G4RadioactiveDecayBase
{
  public: // with description

    G4Radioactivation(const G4String& processName="Radioactivation");
    ~G4Radioactivation();

    virtual void ProcessDescription(std::ostream& outFile) const;

    // Set the decay biasing scheme using the data in "filename"
    void SetDecayBias(G4String filename);

    // Set the half-life threshold for isomer production
    void SetHLThreshold(G4double hl) {halflifethreshold = hl;}

    void SetSourceTimeProfile(G4String filename);
    // Set source exposure function using histograms in "filename"

    G4bool IsRateTableReady(const G4ParticleDefinition &);
    // Returns true if the coefficient and decay time table for all the
    // descendants of the specified isotope are ready.
    // used in VR decay mode only

    void CalculateChainsFromParent(const G4ParticleDefinition&);
    // Calculates the coefficient and decay time table for all the descendents
    // of the specified isotope.  Adds the calculated table to the private data
    // member "theParentChainTable".
    // used in VR decay mode only 

    void GetChainsFromParent(const G4ParticleDefinition&);
    // Used to retrieve the coefficient and decay time table for all the
    // descendants of the specified isotope from "theParentChainTable"
    // and place it in "chainsFromParent".
    // used in VR decay mode only 

    void SetDecayRate(G4int,G4int,G4double, G4int, std::vector<G4double>,
                      std::vector<G4double>);
    // Sets "theDecayRate" with data supplied in the arguements.
    // used in VR decay mode only 

    std::vector<G4RadioactivityTable*> GetTheRadioactivityTables()
       {return theRadioactivityTables;}
    // Return vector of G4Radioactivity map - should be used in VR mode only


    inline void  SetVerboseLevel(G4int value) {verboseLevel = value;}
    // Sets the VerboseLevel which controls duggering display

    inline G4int GetVerboseLevel() const {return verboseLevel;}
    // Returns the VerboseLevel which controls level of debugging output

     // Sets whether branching ration bias scheme applies.
    inline void SetBRBias(G4bool r) {
      BRBias = r;
     }

    // Sets the number of times a nucleus will decay when biased
    inline void SetSplitNuclei(G4int r) {
      NSplit = r;
    }

    //  Returns the nuclear splitting number
    inline G4int GetSplitNuclei () {return NSplit;}

    G4VParticleChange* DecayIt(const G4Track& theTrack,
                               const G4Step&  theStep);

  protected:

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

    G4RadioactivationMessenger* theRadioactivationMessenger;

  private:

    G4bool BRBias;
    G4int NSplit;

    G4double halflifethreshold;

    G4int NSourceBin;
    G4double SBin[100];
    G4double SProfile[100];
    G4int NDecayBin;
    G4double DBin[100];
    G4double DProfile[100];

    G4RadioactiveDecayRatesToDaughter ratesToDaughter;
    G4RadioactiveDecayRates theDecayRateVector;
    G4RadioactiveDecayChainsFromParent chainsFromParent;
    G4RadioactiveDecayParentChainTable theParentChainTable;

    // for the radioactivity tables
    std::vector<G4RadioactivityTable*> theRadioactivityTables;
    G4int decayWindows[100];

    // Remainder of life time at rest
//    G4double fRemainderLifeTime;
    G4int verboseLevel;

};

#endif

