// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HadronElasticProcess.hh,v 1.4 2000-12-14 08:47:20 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Hadron Elastic Scattering Process -- header file
// F.W. Jones, TRIUMF, 04-JUN-96
//  
// For further comments see G4HadronElasticProcess.cc
//
// Modified by FWJ 03-DEC-96: uses G4LCrossSectionData
// Modified by FWJ    FEB-97: adapted to new tracking design
// Modified by FWJ    MAR-97: newer new tracking design
//
// Modified by FWJ 27-MAR-97: first version for Alpha release
// 14-APR-97 FWJ: cross section data class name changed
//
// changing design to use Model-scheme: HPW, 20th June 1997
//
// 14-APR-98 FWJ: variant G4HadronElastic process for
// G4CrossSectionDataSet/DataStore class design.
// 29-JUN-98 FWJ: default data set G4HadronCrossSections
//

// Class Description
// Process for hadron nuclear elastic scattering; 
// to be used in your physics list in case you need this physics.
// Class Description - End


#ifndef G4HadronElasticProcess_h
#define G4HadronElasticProcess_h 1
 
//#include <rw/tphdict.h>
#include "globals.hh"
#include "G4HadronicProcess.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4HadronElasticDataSet.hh"
//#include "Randomize.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
//#include "G4ElementTable.hh"
#include "G4VDiscreteProcess.hh"
//#include "G4PhysicsTable.hh"
//#include "G4PhysicsVector.hh"
//#include "G4LPhysicsFreeVector.hh"
//#include "G4HadronCrossSections.hh"
#include "G4LightMedia.hh"
#include "G4Step.hh"
#include "G4TrackStatus.hh"


class G4HadronElasticProcess : public G4HadronicProcess
{
public:

   G4HadronElasticProcess(const G4String& processName = "LElastic");

   ~G4HadronElasticProcess();
 
   inline
   G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize,
                            G4ForceCondition* condition)
   {
      return GetMeanFreePathBasic(aTrack.GetDynamicParticle(),
                                  aTrack.GetMaterial());
   }

   G4double GetMeanFreePathBasic(const G4DynamicParticle* aParticle,
                                 const G4Material* aMaterial);

   G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
   {
      SetDispatch(this);
      return G4HadronicProcess::GeneralPostStepDoIt(aTrack, aStep);
   }

   void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);

   void SetCrossSectionDataStore(G4CrossSectionDataStore* aDataStore)
   {
      theCrossSectionDataStore = aDataStore;
   }

   G4CrossSectionDataStore* GetCrossSectionDataStore()
   {
      return theCrossSectionDataStore;
   }

  //   inline
  //   void setCrossSectionDataVerboseLevel(G4int value)
  //   {
  //      theCrossSectionData.SetVerboseLevel(value);
  //   }

   G4bool IsApplicable(const G4ParticleDefinition& aParticleType);

   void DumpPhysicsTable(const G4ParticleDefinition& aParticleType);

private:

// Note: DEC compiler won't accept "const G4ParticleDefinition" here
// because it instantiates the ClearAndDestroy method that deletes the 
// pointers in the dictionary.
  //   G4RWTPtrHashDictionary<G4ParticleDefinition, G4PhysicsTable> 
  //     thePhysicsDictionary;

  //   G4HadronCrossSections theCrossSectionData;

   G4CrossSectionDataStore* theCrossSectionDataStore;

   G4double GetMicroscopicCrossSection(const G4DynamicParticle* aParticle,
                                       const G4Element* anElement);

};
#endif
