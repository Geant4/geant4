// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HadronInelasticProcess.hh,v 1.3 2000-09-08 08:41:46 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
      aScaleFactor = 0;
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
     const G4Element *anElement);

    G4VParticleChange *PostStepDoIt(
     const G4Track &aTrack, const G4Step &aStep )
    {
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
 };
 
#endif
 
