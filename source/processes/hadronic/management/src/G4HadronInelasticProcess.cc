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
 // Hadronic Inelastic Process Class
 // J.L. Chuma, TRIUMF, 24-Mar-1997
 // Last modified: 27-Mar-1997
 // J.P. Wellisch: Bug hunting, 23-Apr-97
 // Modified by J.L.Chuma 8-Jul-97 to eliminate possible division by zero for sigma
//
// 14-APR-98 F.W.Jones: variant G4HadronInelastic process for
// G4CrossSectionDataSet/DataStore class design.
//
// 17-JUN-98 F.W.Jones: removed extraneous code causing core dump.
//
 
#include "G4HadronInelasticProcess.hh"
#include "G4GenericIon.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4HadronicException.hh"
  
 void G4HadronInelasticProcess::BuildThePhysicsTable()
  {
    if (!G4HadronicProcess::GetCrossSectionDataStore()) {
      return;
    }
    G4HadronicProcess::GetCrossSectionDataStore()->BuildPhysicsTable(*theParticle);
  }
 
 G4HadronInelasticProcess::G4HadronInelasticProcess(
  const G4String &processName,
  G4ParticleDefinition *aParticle ) :
   G4HadronicProcess( processName )
 {
   G4HadronicProcess::AddDataSet(new G4HadronInelasticDataSet);
   theParticle = aParticle;
 }

 G4HadronInelasticProcess::~G4HadronInelasticProcess() { }

 G4VParticleChange *G4HadronInelasticProcess::
 PostStepDoIt(const G4Track &aTrack, const G4Step &aStep)
 {
   if(0==GetLastCrossSection()&&!getenv("DebugNeutronHP"))
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

 G4bool G4HadronInelasticProcess::
 IsApplicable(const G4ParticleDefinition& aP)
 {
    return  theParticle == &aP || theParticle == G4GenericIon::GenericIon();
 }

 G4double G4HadronInelasticProcess::GetMicroscopicCrossSection(
  const G4DynamicParticle *aParticle,
  const G4Element *anElement,
  G4double aTemp)
  {
    // returns the microscopic cross section in GEANT4 internal units
    
   if (!G4HadronicProcess::GetCrossSectionDataStore()) 
   {
      throw G4HadronicException(__FILE__, __LINE__, 
      "G4HadronInelasticProcess::GetMicroscopicCrossSection: "
                  "no CrossSectionDataStore");
      return DBL_MIN;
   }
   return G4HadronicProcess::GetCrossSectionDataStore()->GetCrossSection(aParticle, anElement, aTemp);

  }
 
 /* end of file */
