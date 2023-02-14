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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// G4HadronicProcess
//
// This is the top level Hadronic Process class
// The inelastic, elastic, capture, and fission processes
// should derive from this class
//
// original by H.P.Wellisch
// J.L. Chuma, TRIUMF, 10-Mar-1997
// Last modified: 04-Apr-1997
// 19-May-2008 V.Ivanchenko cleanup and added comments
// 05-Jul-2010 V.Ivanchenko cleanup commented lines
// 28-Jul-2012 M.Maire add function GetTargetDefinition() 
// 14-Sep-2012 Inherit from RestDiscrete, use subtype code (now in ctor) to
//		configure base-class
// 28-Sep-2012 M. Kelsey -- Undo inheritance change, keep new ctor

#ifndef G4HadronicProcess_h
#define G4HadronicProcess_h 1
 
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4EnergyRangeManager.hh"
#include "G4Nucleus.hh" 
#include "G4ReactionProduct.hh"
#include "G4HadronicProcessType.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4HadXSTypes.hh"
#include <vector>

class G4Track;
class G4Step;
class G4Element;
class G4ParticleChange;
class G4HadronicInteraction;
class G4HadronicProcessStore;
class G4VCrossSectionDataSet;
class G4VLeadingParticleBiasing;
class G4ParticleDefinition;

class G4HadronicProcess : public G4VDiscreteProcess
{
public:
  G4HadronicProcess(const G4String& processName="Hadronic",
		    G4ProcessType procType=fHadronic);    

  // Preferred signature for subclasses, specifying their subtype here
  G4HadronicProcess(const G4String& processName, 
		    G4HadronicProcessType subType);    

  ~G4HadronicProcess() override;

  // register generator of secondaries
  void RegisterMe(G4HadronicInteraction* a);

  // get cross section per element
  G4double GetElementCrossSection(const G4DynamicParticle * part, 
				  const G4Element * elm, 
				  const G4Material* mat = nullptr);

  // obsolete method to get cross section per element
  inline
  G4double GetMicroscopicCrossSection(const G4DynamicParticle * part, 
				      const G4Element * elm, 
				      const G4Material* mat = nullptr);

  // initialisation for a new track
  void StartTracking(G4Track* track) override;

  // compute step limit
  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
           G4double, G4ForceCondition*) override;

  // generic PostStepDoIt recommended for all derived classes
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
				  const G4Step& aStep) override;

  // initialisation of physics tables and G4HadronicProcessStore
  void PreparePhysicsTable(const G4ParticleDefinition&) override;

  // build physics tables and print out the configuration of the process
  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  // dump physics tables 
  void DumpPhysicsTable(const G4ParticleDefinition& p);

  // add cross section data set
  void AddDataSet(G4VCrossSectionDataSet * aDataSet);

  // access to the list of hadronic interactions
  std::vector<G4HadronicInteraction*>& GetHadronicInteractionList();

  // access to an hadronic interaction by name
  G4HadronicInteraction* GetHadronicModel(const G4String&);

  // access to the chosen generator
  inline G4HadronicInteraction* GetHadronicInteraction() const;
  
  // get inverse cross section per volume
  G4double GetMeanFreePath(const G4Track &aTrack, G4double, 
			   G4ForceCondition *) override;

  // access to the target nucleus
  inline const G4Nucleus* GetTargetNucleus() const;

  inline G4Nucleus* GetTargetNucleusPointer();
  
  inline const G4Isotope* GetTargetIsotope();

  // methods needed for implementation of integral XS
  G4double ComputeCrossSection(const G4ParticleDefinition*,
                               const G4Material*,
                               const G4double kinEnergy);

  inline G4HadXSType CrossSectionType() const;
  inline void SetCrossSectionType(G4HadXSType val);

  void ProcessDescription(std::ostream& outFile) const override;
 
  // scale cross section
  void BiasCrossSectionByFactor(G4double aScale);
  void MultiplyCrossSectionBy(G4double factor);
  inline G4double CrossSectionFactor() const;

  // Integral option 
  inline void SetIntegral(G4bool val);

  // Energy-momentum non-conservation limits and reporting
  inline void SetEpReportLevel(G4int level);
  inline void SetEnergyMomentumCheckLevels(G4double relativeLevel,
                                           G4double absoluteLevel);
  inline std::pair<G4double, G4double> GetEnergyMomentumCheckLevels() const;

  // access to the cross section data store
  inline G4CrossSectionDataStore* GetCrossSectionDataStore();

  // access to the data for integral XS method
  inline std::vector<G4TwoPeaksHadXS*>* TwoPeaksXS() const;
  inline std::vector<G4double>* EnergyOfCrossSectionMax() const;

  // hide assignment operator as private 
  G4HadronicProcess& operator=(const G4HadronicProcess& right) = delete;
  G4HadronicProcess(const G4HadronicProcess&) = delete;

protected:

  // generic method to choose secondary generator 
  // recommended for all derived classes
  inline G4HadronicInteraction* ChooseHadronicInteraction(
      const G4HadProjectile & aHadProjectile, G4Nucleus& aTargetNucleus,
      const G4Material* aMaterial, const G4Element* anElement);
              
  // access to the cross section data set
  inline G4double GetLastCrossSection();

  // fill result
  void FillResult(G4HadFinalState* aR, const G4Track& aT);

  void DumpState(const G4Track&, const G4String&, G4ExceptionDescription&);

  // Check the result for catastrophic energy non-conservation
  G4HadFinalState* CheckResult(const G4HadProjectile& thePro,
			       const G4Nucleus& targetNucleus, 
			       G4HadFinalState* result);

  // Check 4-momentum balance
  void CheckEnergyMomentumConservation(const G4Track&, const G4Nucleus&);

private:

  void InitialiseLocal();
  void UpdateCrossSectionAndMFP(const G4double kinEnergy);
  void RecomputeXSandMFP(const G4double kinEnergy);

  inline void DefineXSandMFP();
  inline void ComputeXSandMFP();

  G4double XBiasSurvivalProbability();
  G4double XBiasSecondaryWeight();

  // Set E/p conservation check levels from environment variables
  void GetEnergyMomentumCheckEnvvars();

protected:

  G4HadProjectile thePro;

  G4ParticleChange* theTotalResult;
  G4CrossSectionDataStore* theCrossSectionDataStore;

  G4double fWeight = 1.0;
  G4double aScaleFactor = 1.0;
  G4double theLastCrossSection = 0.0;
  G4double mfpKinEnergy = DBL_MAX;
  G4long epReportLevel = 0;

  G4HadXSType fXSType = fHadNoIntegral;

private:
    
  G4EnergyRangeManager theEnergyRangeManager;
  G4Nucleus targetNucleus;
    
  G4HadronicInteraction* theInteraction = nullptr;
  G4HadronicProcessStore* theProcessStore;
  const G4HadronicProcess* masterProcess = nullptr;
  const G4ParticleDefinition* firstParticle = nullptr;
  const G4ParticleDefinition* currentParticle = nullptr;
  const G4Material* currentMat = nullptr;
  const G4DynamicParticle* fDynParticle = nullptr;

  std::vector<G4double>* theEnergyOfCrossSectionMax = nullptr;
  std::vector<G4TwoPeaksHadXS*>* fXSpeaks = nullptr;
     
  G4double theMFP = DBL_MAX;
  G4double minKinEnergy;

  // counters
  G4int nMatWarn = 0;
  G4int nKaonWarn = 0;
  G4int nICelectrons = 0;
  G4int matIdx = 0;

  // flags
  G4bool levelsSetByProcess = false;
  G4bool G4HadronicProcess_debug_flag = false;
  G4bool useIntegralXS = true;
  G4bool isMaster = true;

  G4ThreeVector unitVector;

  // Energy-momentum checking
  std::pair<G4double, G4double> epCheckLevels;
  std::vector<G4VLeadingParticleBiasing*> theBias;
};

inline G4double G4HadronicProcess::
GetMicroscopicCrossSection(const G4DynamicParticle * part, 
			   const G4Element * elm, 
			   const G4Material* mat)
{ 
  return GetElementCrossSection(part, elm, mat);
}

inline G4HadronicInteraction* 
G4HadronicProcess::GetHadronicInteraction() const
{ 
  return theInteraction;
}

inline const G4Nucleus*
G4HadronicProcess::GetTargetNucleus() const
{
  return &targetNucleus;
}
  
inline const G4Isotope* G4HadronicProcess::GetTargetIsotope()
{ 
  return targetNucleus.GetIsotope();
}

inline G4HadXSType 
G4HadronicProcess::CrossSectionType() const
{ 
  return fXSType;
}

inline void 
G4HadronicProcess::SetCrossSectionType(G4HadXSType val)
{ 
  fXSType = val;
}

inline G4double G4HadronicProcess::CrossSectionFactor() const
{ 
  return aScaleFactor;
}

inline void G4HadronicProcess::SetIntegral(G4bool val)
{ 
  useIntegralXS = val;
}

inline void G4HadronicProcess::SetEpReportLevel(G4int level)
{ 
  epReportLevel = level;
}

inline void 
G4HadronicProcess::SetEnergyMomentumCheckLevels(G4double relativeLevel, 
                                                G4double absoluteLevel)
{ 
  epCheckLevels.first = relativeLevel;
  epCheckLevels.second = absoluteLevel;
  levelsSetByProcess = true;
}

inline std::pair<G4double, G4double>
G4HadronicProcess::GetEnergyMomentumCheckLevels() const
{ 
  return epCheckLevels;
}

inline G4CrossSectionDataStore*
G4HadronicProcess::GetCrossSectionDataStore()
{
  return theCrossSectionDataStore;
}

inline std::vector<G4TwoPeaksHadXS*>* 
G4HadronicProcess::TwoPeaksXS() const
{ 
  return fXSpeaks;
}

inline std::vector<G4double>*
G4HadronicProcess::EnergyOfCrossSectionMax() const
{
  return  theEnergyOfCrossSectionMax;
}

inline G4HadronicInteraction* G4HadronicProcess::
ChooseHadronicInteraction(const G4HadProjectile& aHadProjectile,
                          G4Nucleus& aTargetNucleus,
                          const G4Material* aMaterial,
                          const G4Element* anElement)
{ 
  return theEnergyRangeManager.GetHadronicInteraction(aHadProjectile, 
                                                      aTargetNucleus,
                                                      aMaterial,anElement);
}

inline G4Nucleus* G4HadronicProcess::GetTargetNucleusPointer() 
{ 
  return &targetNucleus;
}

inline G4double G4HadronicProcess::GetLastCrossSection() 
{ 
  return theLastCrossSection;
}

inline void G4HadronicProcess::DefineXSandMFP()
{
  theLastCrossSection = aScaleFactor*
    theCrossSectionDataStore->GetCrossSection(fDynParticle, currentMat);
  theMFP = (theLastCrossSection > 0.0) ? 1.0/theLastCrossSection : DBL_MAX;
}

inline void G4HadronicProcess::ComputeXSandMFP()
{
  theLastCrossSection = aScaleFactor*
    theCrossSectionDataStore->ComputeCrossSection(fDynParticle, currentMat);
  theMFP = (theLastCrossSection > 0.0) ? 1.0/theLastCrossSection : DBL_MAX;
}
 
#endif
 
