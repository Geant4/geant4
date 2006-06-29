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
// G4 Fission Process -- header file
// F.W. Jones, TRIUMF, 03-DEC-96
//  
// For further comments see G4HadronFissionProcess.cc.
//
// 27-MAR-97 F.W.Jones: first version for Alpha release
// 14-APR-97 FWJ: cross section data class name changed
//
// 19-MAY-98 FWJ: variant G4HadronFission process for
// G4CrossSectionDataSet/DataStore class design.
// 29-JUN-98 FWJ: default data set G4HadronCrossSections
//


// Class Description
// Process for induced fission; 
// to be used in your physics list in case you need this physics.
// Class Description - End

#ifndef G4HadronFissionProcess_h
#define G4HadronFissionProcess_h 1
 
//#include <rw/tphdict.h>
#include "globals.hh"
//#include "Randomize.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
//#include "G4ElementTable.hh"
#include "G4VDiscreteProcess.hh"
//#include "G4PhysicsTable.hh"
//#include "G4PhysicsVector.hh"
//#include "G4LPhysicsFreeVector.hh"
//#include "G4HadronCrossSections.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Step.hh"
#include "G4TrackStatus.hh"
#include "G4HadronicProcess.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4HadronFissionDataSet.hh"


class G4HadronFissionProcess : public G4HadronicProcess
{
public:

   G4HadronFissionProcess(const G4String& processName ="LFission");

   ~G4HadronFissionProcess();
 
   G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
   {
      SetDispatch(this);
      return G4HadronicProcess::GeneralPostStepDoIt(aTrack, aStep);
   }

   G4bool IsApplicable(const G4ParticleDefinition& )
   {
      return true;
   }

   void BuildThePhysicsTable(G4ParticleDefinition& aParticleType);

   void DumpPhysicsTable(const G4ParticleDefinition& aParticleType);

private:

   G4double GetMicroscopicCrossSection(const G4DynamicParticle* aParticle,
                                       const G4Element* anElement, 
				       G4double aTemp);
};
#endif
