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
//
 // Hadronic Inelastic Process class
 // The specific particle inelastic processes derive from this class
 // This is an abstract base class, since the pure virtual function
 // PostStepDoIt has not been defined yet.
 //
 // J.L. Chuma, TRIUMF, 10-Mar-1997
 // Last modified: 27-Mar-1997
//
// 14-APR-98 F.W.Jones: variant G4HadronInelastic process for
// G4CrossSectionDataSet/DataStore class design.
// 29-JUN-98 F.W.Jones: default data set G4HadronCrossSections
//

#ifndef G4HadronInelasticProcess_h
#define G4HadronInelasticProcess_h 1

#include "G4HadronicProcess.hh"
//#include "G4LPhysicsFreeVector.hh"
#include "G4HadronCrossSections.hh" 
#include "G4CrossSectionDataStore.hh"
#include "G4HadronInelasticDataSet.hh"
#include "G4ParticleChange.hh"
 

 class G4HadronInelasticProcess : public G4HadronicProcess
 {
 public:
    
    G4HadronInelasticProcess(
     const G4String &processName,
     G4ParticleDefinition *aParticle ) :
      G4HadronicProcess( processName ),
      theCrossSectionDataStore(new G4CrossSectionDataStore)
    {
      theCrossSectionDataStore->AddDataSet(new G4HadronInelasticDataSet);
      theParticle = aParticle;
      aScaleFactor = 1;
      //      BuildThePhysicsTable();
    }
    
    virtual ~G4HadronInelasticProcess()
    { }
    
    G4double GetMeanFreePath(
     const G4Track &aTrack,
     G4double previousStepSize,
     G4ForceCondition *condition );
    
    void BuildThePhysicsTable();
    
    void BiasCrossSectionByFactor(G4double aScale) {aScaleFactor = aScale;}
    
    G4double GetMicroscopicCrossSection(
     const G4DynamicParticle *aParticle,
     const G4Element *anElement,
     G4double aTemp);

    G4VParticleChange *PostStepDoIt(
     const G4Track &aTrack, const G4Step &aStep )
    {
      if(0==theLastCrossSection&&!getenv("DebugNeutronHP"))
      {
        G4cerr << "G4HadronInelasticProcess: called for final state, while cross-section was zero"<<G4endl;
	G4cerr << "                          Returning empty particle change...."<<G4endl;
	G4double dummy=0;
	G4ForceCondition condition;
	G4double it = GetMeanFreePath(aTrack, dummy, &condition);
	G4cerr << "                          current MeanFreePath is "<<it<<G4endl;
	theParticleChange.Initialize(aTrack);
	return &theParticleChange;
      }
      SetDispatch( this );
      return G4HadronicProcess::GeneralPostStepDoIt( aTrack, aStep );
    }
    
   void SetCrossSectionDataStore(G4CrossSectionDataStore* aDataStore)
   {
      theCrossSectionDataStore = aDataStore;
   }

   G4CrossSectionDataStore* GetCrossSectionDataStore()
   {
      return theCrossSectionDataStore;
   }

 protected:

   //    G4HadronicCrossSections theCrossSectionData;
   G4CrossSectionDataStore* theCrossSectionDataStore;

    G4ParticleDefinition *theParticle;
    
 private:
   G4double aScaleFactor;
   G4double theLastCrossSection;
   G4ParticleChange theParticleChange;
 };
 
#endif
 
