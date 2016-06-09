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
//
 // This is the top level Hadronic Process class
 // The inelastic, elastic, capture, and fission processes
 // should derive from this class
 // This is an abstract base class, since the pure virtual function
 // PostStepDoIt has not been defined yet.
 // Note:  there is no .cc file
 //
 // original by H.P.Wellisch
 // J.L. Chuma, TRIUMF, 10-Mar-1997
 // Last modified: 04-Apr-1997
 
#ifndef G4HadronicProcess_h
#define G4HadronicProcess_h 1
 
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4EnergyRangeManager.hh"
#include "G4Nucleus.hh" 
#include "G4ReactionProduct.hh"
#include <vector>
#include "G4VIsotopeProduction.hh"
#include "G4IsoParticleChange.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4VLeadingParticleBiasing.hh"
#include "G4Delete.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4HadronicException.hh"
#include "G4Fuzzy.hh"

class G4Track;
class G4Step;
class G4Element;
class G4ParticleChange;

 class G4HadronicProcess : public G4VDiscreteProcess
 {
 public:
    
    G4HadronicProcess( const G4String &processName = "Hadronic", 
                       G4ProcessType   aType = fHadronic );    
    virtual ~G4HadronicProcess();

    void RegisterMe( G4HadronicInteraction *a );

    void AddDataSet(G4VCrossSectionDataSet * aDataSet)
    {
       theCrossSectionDataStore->AddDataSet(aDataSet);
    }
        
    virtual G4VParticleChange *PostStepDoIt( const G4Track &aTrack, 
                                            const G4Step &aStep ) = 0;
        
    virtual G4double GetMicroscopicCrossSection( const G4DynamicParticle *aParticle, 
                                                 const G4Element *anElement, 
						 G4double aTemp ) = 0;
    
    G4double GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *);

    // Set methods for isotope production
    
    static void EnableIsotopeProductionGlobally();
    static void DisableIsotopeProductionGlobally();
    
    void EnableIsotopeCounting()  {isoIsOnAnyway = 1;}
    void DisableIsotopeCounting() {isoIsOnAnyway = -1;}
    
    void RegisterIsotopeProductionModel(G4VIsotopeProduction * aModel)
    { theProductionModels.push_back(aModel); }

    static G4IsoParticleChange * GetIsotopeProductionInfo() 
    { 
      G4IsoParticleChange * anIsoResult = theIsoResult;
      if(theIsoResult) theOldIsoResult = theIsoResult;
      theIsoResult = 0;
      return anIsoResult;
    }

    // use this with bool only
    static G4bool AlwaysKillLeadingHadron(G4Fuzzy aB=G4Fuzzy())
    {
       static G4Fuzzy state = G4Fuzzy();
       if(getenv("AlwaysKillLeadingHadron")) return true;
       if (!aB.first) state = aB;
       return state.second;
    }

    void BiasCrossSectionByFactor(G4double aScale) 
    {
      xBiasOn = true;
      aScaleFactor = aScale;
      G4String it = GetProcessName(); 
      if( (it != "PhotonInelastic") && 
          (it != "ElectroNuclear") && 
	  (it != "PositronNuclear") )
      {
        G4Exception("G4HadronicProcess", "007", FatalException,
                    "Cross-section biasing available only for gamma and electro nuclear reactions.");
      }
      if(aScale<100)
      {
        G4Exception("G4HadronicProcess", "001", JustWarning,
                    "Cross-section bias readjusted to be above safe limit. New value is 100");        
        aScaleFactor = 100.;
      }
    }
    
 protected:
    
    virtual void ResetNumberOfInteractionLengthLeft()
    {
      G4VProcess::theNumberOfInteractionLengthLeft =  -std::log( G4UniformRand() );
      theInitialNumberOfInteractionLength = G4VProcess::theNumberOfInteractionLengthLeft;
      // hpw ReStarted = true;
    }

    G4VParticleChange *GeneralPostStepDoIt( const G4Track &aTrack, 
                                           const G4Step &aStep );
    
    void SetDispatch( G4HadronicProcess *value )
    { dispatch=value; }
    
    G4Element * ChooseAandZ( const G4DynamicParticle *aParticle,
                             const G4Material *aMaterial );

    inline const G4EnergyRangeManager &GetEnergyRangeManager() const
    { return theEnergyRangeManager; }
    
    inline void SetEnergyRangeManager( const G4EnergyRangeManager &value )
    { theEnergyRangeManager = value; }

    inline G4HadronicInteraction *ChooseHadronicInteraction(
     G4double kineticEnergy, G4Material *aMaterial, G4Element *anElement )
    {
      G4EnergyRangeManager* ERMan = GetManagerPointer();
      if(!ERMan->GetHadronicInteractionCounter())
        G4cout<<"*G4HadronicProcess::ChooseHadronicInteraction: process = "
              <<GetProcessName()<<", nM="<<ERMan->GetHadronicInteractionCounter()<<G4endl;
      return ERMan->GetHadronicInteraction(kineticEnergy, aMaterial, anElement);
    }

    inline G4HadronicInteraction *GetHadronicInteraction()
    { return theInteraction; }
    
 public:
    inline G4EnergyRangeManager *GetManagerPointer()
    { return &theEnergyRangeManager; }
 protected:

    G4double GetCurrentZ()
    { return currentZ; }
    
    G4double GetCurrentN()
    { return currentN; }
    
    G4CrossSectionDataStore* GetCrossSectionDataStore()
    {
       return theCrossSectionDataStore;
    }
   
    G4double GetLastCrossSection() {return theLastCrossSection;}
 private:
    
    G4HadFinalState * DoIsotopeCounting(G4HadFinalState * aResult,
                                          const G4Track & aTrack,
                                          const G4Nucleus & aNucleus);
                                          
    G4IsoResult * ExtractResidualNucleus(const G4Track & aTrack,
                                         const G4Nucleus & aNucleus,
                                         G4HadFinalState * aResult);

    G4double GetTotalNumberOfInteractionLengthTraversed()
    {
      return theInitialNumberOfInteractionLength
            -G4VProcess::theNumberOfInteractionLengthLeft;
    }
    
    G4double GetDistanceToBoundary(const G4Track & aT);

    void FillTotalResult(G4HadFinalState * aR, const G4Track & aT);
    
    void SetCrossSectionDataStore(G4CrossSectionDataStore* aDataStore)
    {
       theCrossSectionDataStore = aDataStore;
    }
    
    G4double XBiasSurvivalProbability();
    G4double XBiasSecondaryWeight();
    
 private:
    
    G4EnergyRangeManager theEnergyRangeManager;
    
    G4HadronicInteraction *theInteraction;

    G4CrossSectionDataStore* theCrossSectionDataStore;
 
    G4Nucleus targetNucleus;
    
    G4double currentZ;
    G4double currentN;
    G4HadronicProcess *dispatch;

// swiches for isotope production
    
    static G4bool isoIsEnabled; // true or false; local swich overrides
    G4int isoIsOnAnyway; // true(1), false(-1) or default(0)
    
    G4IsoParticleChange theIsoPC;
    std::vector<G4VIsotopeProduction *> theProductionModels;
    
    std::vector<G4VLeadingParticleBiasing *> theBias;

    static G4IsoParticleChange * theIsoResult;
    static G4IsoParticleChange * theOldIsoResult;
    
    G4ParticleChange * theTotalResult; 
    
    G4double theInitialNumberOfInteractionLength;   

    G4double aScaleFactor;
    G4bool xBiasOn;
    G4double theLastCrossSection;

    G4int ModelingState;
 };
 
#endif
 
