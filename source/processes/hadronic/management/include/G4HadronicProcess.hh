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
// $Id: G4HadronicProcess.hh,v 1.42 2010-07-05 14:50:15 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//

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

#include "G4CrossSectionDataStore.hh"
#include "G4HadronicProcessType.hh"

class G4Track;
class G4Step;
class G4Element;
class G4ParticleChange;


class G4HadronicProcess : public G4VDiscreteProcess
{
public:
    
  G4HadronicProcess(const G4String& processName = "Hadronic", 
		    G4ProcessType aType = fHadronic);    

  virtual ~G4HadronicProcess();

  // register generator of secondaries
  void RegisterMe(G4HadronicInteraction* a);

  // get cross section per element
  virtual 
  G4double GetMicroscopicCrossSection(const G4DynamicParticle *aParticle, 
				      const G4Element *anElement, 
				      G4double aTemp );

  // generic PostStepDoIt recommended for all derived classes
  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
					  const G4Step& aStep);

  // initialisation of physics tables and G4HadronicProcessStore
  virtual void PreparePhysicsTable(const G4ParticleDefinition&);

  // build physics tables and print out the configuration of the process
  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  // dump physics tables 
  inline void DumpPhysicsTable(const G4ParticleDefinition& p)
  { theCrossSectionDataStore->DumpPhysicsTable(p); }

  // add cross section data set
  inline void AddDataSet(G4VCrossSectionDataSet * aDataSet)
  { theCrossSectionDataStore->AddDataSet(aDataSet);}

  // access to the manager
  inline G4EnergyRangeManager *GetManagerPointer()
  { return &theEnergyRangeManager; }
          
  // get inverse cross section per volume
  G4double GetMeanFreePath(const G4Track &aTrack, G4double, 
			   G4ForceCondition *);

protected:    

  // reset number of interaction length and save  
  virtual void ResetNumberOfInteractionLengthLeft()
  { G4VProcess::ResetNumberOfInteractionLengthLeft(); 
    theInitialNumberOfInteractionLength = 
      G4VProcess::theNumberOfInteractionLengthLeft;
  }

  // generic method to choose secondary generator 
  // recommended for all derived classes
  inline G4HadronicInteraction *ChooseHadronicInteraction(
      G4double kineticEnergy, G4Material *aMaterial, G4Element *anElement )
  { return theEnergyRangeManager.GetHadronicInteraction(kineticEnergy,
							aMaterial,anElement);
  }

public:

  // Methods for isotope production    
  static void EnableIsotopeProductionGlobally();
  static void DisableIsotopeProductionGlobally();
    
  void EnableIsotopeCounting()  {isoIsOnAnyway = 1;}
  void DisableIsotopeCounting() {isoIsOnAnyway = -1;}
    
  void RegisterIsotopeProductionModel(G4VIsotopeProduction * aModel)
  { theProductionModels.push_back(aModel); }

  static G4IsoParticleChange * GetIsotopeProductionInfo();

  void BiasCrossSectionByFactor(G4double aScale);

  // Energy-momentum non-conservation limits and reporting
  inline void SetEpReportLevel(G4int level)
  { epReportLevel = level; }

  inline void SetEnergyMomentumCheckLevels(G4double relativeLevel, G4double absoluteLevel)
  { epCheckLevels.first = relativeLevel;
    epCheckLevels.second = absoluteLevel;
    levelsSetByProcess = true;
  }

  inline std::pair<G4double, G4double> GetEnergyMomentumCheckLevels() const
  { return epCheckLevels; }

protected:
            
  // obsolete method will be removed
  inline const G4EnergyRangeManager &GetEnergyRangeManager() const
  { return theEnergyRangeManager; }
    
  // obsolete method will be removed
  inline void SetEnergyRangeManager( const G4EnergyRangeManager &value )
  { theEnergyRangeManager = value; }

  // access to the chosen generator
  inline G4HadronicInteraction *GetHadronicInteraction()
  { return theInteraction; }
    
  // access to the cross section data store
  inline G4CrossSectionDataStore* GetCrossSectionDataStore()
  { return theCrossSectionDataStore; }
   
  // access to the cross section data set
  inline G4double GetLastCrossSection() 
  { return theLastCrossSection; }

private:

  void DumpState(const G4Track&, const G4String&);
    
  void FillTotalResult(G4HadFinalState * aR, const G4Track & aT);

  G4HadFinalState * DoIsotopeCounting(G4HadFinalState * aResult,
				      const G4Track & aTrack,
				      const G4Nucleus & aNucleus);
                                          
  G4IsoResult * ExtractResidualNucleus(const G4Track & aTrack,
				       const G4Nucleus & aNucleus,
				       G4HadFinalState * aResult);

  inline G4double GetTotalNumberOfInteractionLengthTraversed()
  { return theInitialNumberOfInteractionLength
      -G4VProcess::theNumberOfInteractionLengthLeft;
  }
            
  G4double XBiasSurvivalProbability();
  G4double XBiasSecondaryWeight();

  void CheckEnergyMomentumConservation(const G4Track&, const G4Nucleus&);
    
private:
    
  G4EnergyRangeManager theEnergyRangeManager;
    
  G4HadronicInteraction *theInteraction;

  G4CrossSectionDataStore* theCrossSectionDataStore;
 
  G4Nucleus targetNucleus;
    
  G4HadronicProcess *dispatch;

  bool G4HadronicProcess_debug_flag;

  // Energy-momentum checking
  G4int epReportLevel;
  std::pair<G4double, G4double> epCheckLevels;
  G4bool levelsSetByProcess;

  // swiches for isotope production    
  static G4bool isoIsEnabled; // true or false; local swich overrides
  G4int isoIsOnAnyway; // true(1), false(-1) or default(0)
    
  G4IsoParticleChange theIsoPC;
  std::vector<G4VIsotopeProduction *> theProductionModels;
    
  std::vector<G4VLeadingParticleBiasing *> theBias;

  static G4IsoParticleChange* theIsoResult;
  static G4IsoParticleChange* theOldIsoResult;
    
  G4ParticleChange* theTotalResult; 
    
  G4double theInitialNumberOfInteractionLength;   

  G4double aScaleFactor;
  G4bool xBiasOn;
  G4double theLastCrossSection;

  G4int ModelingState;
};
 
#endif
 
